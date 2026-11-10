"""Inside one private netns: fixed loopback HTTP byte relay to /model.sock."""

import selectors
import socket
import subprocess
import sys
import threading


def forward(client):
    with client, socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as upstream:
        upstream.connect("/model.sock")
        with selectors.DefaultSelector() as selector:
            selector.register(client, selectors.EVENT_READ, upstream)
            selector.register(upstream, selectors.EVENT_READ, client)
            while True:
                for key, _events in selector.select(timeout=30):
                    data = key.fileobj.recv(65536)
                    if not data:
                        return
                    key.data.sendall(data)


def main():
    listener = socket.socket()
    listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listener.bind(("127.0.0.1", 18765))
    listener.listen(8)

    def accept():
        while True:
            try:
                client, _ = listener.accept()
            except OSError:
                return
            threading.Thread(target=forward, args=(client,), daemon=True).start()

    threading.Thread(target=accept, daemon=True).start()
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
        listener.close()
        if process.poll() is None:
            process.kill()
            process.wait()


if __name__ == "__main__":
    sys.exit(main())
