"""Build Explorer data, start a local server, and open the site in a browser."""

from __future__ import annotations

import argparse
import http.server
import socketserver
import subprocess
import sys
import threading
import webbrowser
from pathlib import Path



class ReuseTCPServer(socketserver.TCPServer):
    allow_reuse_address = True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-build", action="store_true", help="Serve the current static build without regenerating its data.")
    args = parser.parse_args()
    if not args.no_build:
        subprocess.run([sys.executable, "sites/pathogen-atlas-explorer/scripts/build_data.py"], check=True)
    site = Path("sites/pathogen-atlas-explorer/dist").resolve()
    handler = lambda *handler_args, **handler_kwargs: http.server.SimpleHTTPRequestHandler(*handler_args, directory=str(site), **handler_kwargs)
    with ReuseTCPServer(("127.0.0.1", args.port), handler) as server:
        url = f"http://127.0.0.1:{args.port}"
        threading.Timer(.5, lambda: webbrowser.open(url)).start()
        print(f"Pathogen Atlas Explorer: {url} (Ctrl-C to stop)", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("Explorer server stopped.")


if __name__ == "__main__":
    main()
