"""Bounded byte relay inside one private netns, only to the fixed model UDS."""

import selectors
import socket
import subprocess
import sys
import threading
import time

CONNECT_SECONDS = 1.0
IDLE_SECONDS = 2.0
TOTAL_SECONDS = 8.0
MAX_BYTES = 1048576
MAX_BUFFER = 65536
MAX_CONNECTIONS = 2


def forward(
    client, upstream_path="/model.sock", *, idle=IDLE_SECONDS, total=TOTAL_SECONDS
):
    started = time.monotonic()
    transferred = 0
    with client, socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as upstream:
        upstream.settimeout(min(CONNECT_SECONDS, total))
        upstream.connect(upstream_path)
        client.setblocking(False)
        upstream.setblocking(False)
        buffers = {client: bytearray(), upstream: bytearray()}
        peers = {client: upstream, upstream: client}
        last = time.monotonic()
        with selectors.DefaultSelector() as selector:
            while True:
                now = time.monotonic()
                remaining = min(started + total - now, last + idle - now)
                if remaining <= 0:
                    return "deadline"
                for sock, peer in peers.items():
                    events = (
                        selectors.EVENT_READ if len(buffers[peer]) < MAX_BUFFER else 0
                    )
                    if buffers[sock]:
                        events |= selectors.EVENT_WRITE
                    try:
                        selector.unregister(sock)
                    except KeyError:
                        pass
                    if events:
                        selector.register(sock, events)
                for key, events in selector.select(remaining):
                    sock = key.fileobj
                    if events & selectors.EVENT_WRITE:
                        try:
                            count = sock.send(buffers[sock])
                        except BlockingIOError:
                            count = 0
                        if count:
                            del buffers[sock][:count]
                            last = time.monotonic()
                    if events & selectors.EVENT_READ:
                        try:
                            data = sock.recv(MAX_BUFFER - len(buffers[peers[sock]]))
                        except BlockingIOError:
                            continue
                        if not data:
                            return "disconnect"
                        transferred += len(data)
                        if transferred > MAX_BYTES:
                            return "byte_budget"
                        buffers[peers[sock]].extend(data)
                        last = time.monotonic()


class Relay:
    def __init__(self, listener, upstream="/model.sock"):
        self.listener, self.upstream = listener, upstream
        self.stop = threading.Event()
        self.slots = threading.BoundedSemaphore(MAX_CONNECTIONS)
        self.lock = threading.Lock()
        self.clients = set()
        self.threads = []
        listener.settimeout(0.1)
        self.thread = threading.Thread(target=self.accept, daemon=True)
        self.thread.start()

    def accept(self):
        while not self.stop.is_set():
            try:
                client, _ = self.listener.accept()
            except socket.timeout:
                continue
            except OSError:
                return
            if not self.slots.acquire(blocking=False):
                client.close()
                continue
            with self.lock:
                self.clients.add(client)
                self.threads = [t for t in self.threads if t.is_alive()]
                worker = threading.Thread(
                    target=self.handle, args=(client,), daemon=True
                )
                self.threads.append(worker)
            worker.start()

    def handle(self, client):
        try:
            forward(client, self.upstream)
        except OSError:
            pass
        finally:
            client.close()
            with self.lock:
                self.clients.discard(client)
            self.slots.release()

    def close(self):
        self.stop.set()
        self.listener.close()
        with self.lock:
            for client in self.clients:
                try:
                    client.shutdown(socket.SHUT_RDWR)
                except OSError:
                    pass
        self.thread.join(0.5)
        for thread in self.threads:
            thread.join(TOTAL_SECONDS + 0.2)
        if self.thread.is_alive() or any(t.is_alive() for t in self.threads):
            raise RuntimeError("relay close deadline")


def main():
    listener = socket.socket()
    listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listener.bind(("127.0.0.1", 18765))
    listener.listen(MAX_CONNECTIONS)
    relay = Relay(listener)
    process = subprocess.Popen(
        [
            "/runtime/peri",
            "--db-path",
            "/state/threads.db",
            "acp",
            "--cwd",
            "/work",
            "--settings-stdin",
        ],
        close_fds=True,
    )
    try:
        return process.wait()
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=2)
        relay.close()


if __name__ == "__main__":
    sys.exit(main())
