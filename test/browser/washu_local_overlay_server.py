"""Serve the installed HelixWeave page and its real API for an anonymous file viewer."""

from __future__ import annotations

import argparse
from pathlib import Path
import socket

import uvicorn

from encode_pipeline.api.main import create_app
from encode_pipeline.deployment.frontend import ManifestStaticApplication
from encode_pipeline.frontend_assets import load_packaged_frontend_assets


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    args = parser.parse_args()
    app = create_app(
        database_url=f"sqlite:///{args.work / 'anonymous.db'}",
        workspace_root=args.work / "workspace",
    )
    app.mount(
        "/", ManifestStaticApplication(load_packaged_frontend_assets()), name="frontend"
    )
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        listener.listen()
        print(f"WASHU_VIEWER_PORT={listener.getsockname()[1]}", flush=True)
        uvicorn.Server(uvicorn.Config(app, log_level="warning")).run(sockets=[listener])


if __name__ == "__main__":
    main()
