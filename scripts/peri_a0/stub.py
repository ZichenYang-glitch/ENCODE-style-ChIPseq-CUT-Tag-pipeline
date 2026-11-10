"""Non-forwarding Unix HTTP stub. Every parsed request is counted before checks."""

from collections import deque
import http.server
import json
import socketserver
import threading


class ModelStub:
    def __init__(self, path):
        self.path = path
        self.requests, self.failures = [], []
        self.responses = deque()
        self.lock = threading.Lock()
        self.stopping = threading.Event()
        self.on_catalog_failure = None
        owner = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def setup(self):
                super().setup()
                self.connection.settimeout(5)

            def log_message(self, *_args):
                pass

            def parse_request(self):
                valid = super().parse_request()
                with owner.lock:
                    self.record = {
                        "method": getattr(self, "command", None),
                        "path": getattr(self, "path", None),
                    }
                    if len(owner.requests) >= 256:
                        owner.failures.append("request_budget")
                        return False
                    owner.requests.append(self.record)
                    if not valid:
                        owner.failures.append("malformed_http")
                return valid

            def deny_method(self):
                owner.failures.append("unexpected_method")
                self.send_error(405)

            do_GET = do_PUT = do_DELETE = do_CONNECT = do_OPTIONS = do_PATCH = (
                do_HEAD
            ) = deny_method

            def do_POST(self):
                try:
                    if self.path != "/v1/chat/completions":
                        raise ValueError("unexpected_path")
                    if self.headers.get("Authorization") != "Bearer synthetic-a0-key":
                        raise ValueError("invalid_synthetic_credential")
                    length = int(self.headers.get("Content-Length", "0"))
                    if not 0 < length <= 262144:
                        raise ValueError("request_size")
                    body = json.loads(self.rfile.read(length))
                    self.record["body"] = body
                    if not isinstance(body, dict) or not isinstance(
                        body.get("tools"), list
                    ):
                        raise ValueError("malformed_provider_body")
                    names = {item["function"]["name"] for item in body.get("tools", [])}
                    if names != {"SearchExtraTools", "ExecuteExtraTool"}:
                        owner.failures.append("PROVIDER_TOOL_CATALOG_MISMATCH")
                        if owner.on_catalog_failure:
                            owner.on_catalog_failure("PROVIDER_TOOL_CATALOG_MISMATCH")
                        raise ValueError("provider_tool_catalog_mismatch")
                    with owner.lock:
                        if not owner.responses:
                            raise ValueError("unscripted_model_request")
                        response = owner.responses.popleft()
                    if response is None:
                        owner.stopping.wait(10)
                        return
                    self.send_response(200)
                    self.send_header("Content-Type", "text/event-stream")
                    self.end_headers()
                    self.wfile.write(
                        (
                            "data: " + json.dumps(response) + "\n\ndata: [DONE]\n\n"
                        ).encode()
                    )
                    self.wfile.flush()
                except (
                    ValueError,
                    OSError,
                    TypeError,
                    KeyError,
                    AttributeError,
                ) as error:
                    code = (
                        str(error)
                        if isinstance(error, ValueError)
                        else "malformed_provider_or_io"
                    )
                    owner.failures.append(code)
                    if owner.on_catalog_failure:
                        owner.on_catalog_failure("PROVIDER_REQUEST_REJECTED")
                    try:
                        self.send_error(400)
                    except OSError:
                        pass

        class Server(socketserver.ThreadingMixIn, socketserver.UnixStreamServer):
            daemon_threads = False
            block_on_close = True

        self.server = Server(str(path), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def script(self, *actions):
        with self.lock:
            if self.responses:
                raise RuntimeError("previous script was not fully consumed")
            for index, action in enumerate(actions):
                if action is None:
                    self.responses.append(None)
                    continue
                if isinstance(action, str):
                    delta = {"role": "assistant", "content": action}
                    finish = "stop"
                else:
                    name, arguments = action
                    delta = {
                        "role": "assistant",
                        "tool_calls": [
                            {
                                "index": 0,
                                "id": f"call-a0-{len(self.requests)}-{index}",
                                "type": "function",
                                "function": {
                                    "name": name,
                                    "arguments": json.dumps(arguments),
                                },
                            }
                        ],
                    }
                    finish = "tool_calls"
                self.responses.append(
                    {
                        "id": f"a0-{index}",
                        "object": "chat.completion.chunk",
                        "model": "a0-stub",
                        "choices": [
                            {"index": 0, "delta": delta, "finish_reason": finish}
                        ],
                        "usage": {"prompt_tokens": 5, "completion_tokens": 3},
                    }
                )

    def close(self):
        self.stopping.set()
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)
        self.path.unlink(missing_ok=True)
