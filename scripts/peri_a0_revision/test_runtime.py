"""Bounded transport behavior; real sockets, no model accounts or platform services."""

from pathlib import Path
import socket
import tempfile
import threading
import time
import unittest
from relay import forward, Relay, MAX_CONNECTIONS
from stub import ModelStub


class RuntimeTests(unittest.TestCase):
    def forward_case(self, mode):
        with tempfile.TemporaryDirectory(prefix="a0r-") as directory:
            path = str(Path(directory) / "up.sock")
            listener = socket.socket(socket.AF_UNIX)
            listener.bind(path)
            listener.listen(1)
            left, right = socket.socketpair()
            outcomes = []

            def worker():
                try:
                    outcomes.append(forward(left, path, idle=0.15, total=0.4))
                except OSError:
                    outcomes.append("disconnect")

            thread = threading.Thread(target=worker)
            thread.start()
            upstream, _ = listener.accept()
            start = time.monotonic()
            try:
                if mode == "disconnect":
                    right.close()
                elif mode == "blocked_write":
                    right.setblocking(False)
                    while thread.is_alive() and time.monotonic() - start < 0.5:
                        try:
                            right.send(b"x" * 65536)
                        except BlockingIOError:
                            time.sleep(0.005)
                thread.join(0.8)
                self.assertFalse(thread.is_alive())
                self.assertLess(time.monotonic() - start, 0.9)
                self.assertIn(outcomes[0], ["deadline", "disconnect", "byte_budget"])
            finally:
                right.close()
                upstream.close()
                listener.close()

    def test_idle_absolute_and_write_backpressure_deadlines(self):
        self.forward_case("idle")
        self.forward_case("blocked_write")

    def test_disconnect_closes_forwarder(self):
        self.forward_case("disconnect")

    def test_connection_cap_and_close(self):
        with tempfile.TemporaryDirectory(prefix="a0r-") as directory:
            upstream = socket.socket(socket.AF_UNIX)
            path = str(Path(directory) / "up.sock")
            upstream.bind(path)
            upstream.listen(4)
            listener = socket.socket()
            listener.bind(("127.0.0.1", 0))
            listener.listen(4)
            relay = Relay(listener, path)
            clients = []
            accepted = []
            try:
                for _ in range(MAX_CONNECTIONS):
                    clients.append(
                        socket.create_connection(listener.getsockname(), timeout=1)
                    )
                    accepted.append(upstream.accept()[0])
                extra = socket.create_connection(listener.getsockname(), timeout=1)
                self.assertEqual(extra.recv(1), b"")
                extra.close()
                self.assertLessEqual(len(relay.threads), MAX_CONNECTIONS)
                relay.close()
                self.assertTrue(all(not t.is_alive() for t in relay.threads))
            finally:
                relay.close()
                for connection in clients + accepted:
                    connection.close()
                upstream.close()

    def test_stub_rejects_oversized_body_and_closes(self):
        with tempfile.TemporaryDirectory(prefix="a0r-") as directory:
            stub = ModelStub(Path(directory) / "stub.sock")
            client = socket.socket(socket.AF_UNIX)
            try:
                client.connect(str(stub.path))
                client.settimeout(2)
                client.sendall(
                    b"POST /v1/chat/completions HTTP/1.0\r\nAuthorization: Bearer synthetic-a0-key\r\nContent-Length: 262145\r\n\r\n"
                )
                self.assertIn(b"400", client.recv(1024))
                self.assertIn("request_size", stub.failures)
            finally:
                client.close()
                stub.close()
            self.assertFalse(stub.path.exists())
            self.assertFalse(stub.thread.is_alive())


if __name__ == "__main__":
    unittest.main(verbosity=2)
