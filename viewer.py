from __future__ import annotations

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Embedded 3D protein viewer (Mol* via WebView).")
    parser.add_argument("--url", default="https://molstar.org/viewer/", help="URL to open inside embedded viewer")
    parser.add_argument("--pdb", default="", help="Optional PDB ID to load in Mol*")
    args = parser.parse_args()

    url = args.url
    if args.pdb:
        url = f"https://molstar.org/viewer/?pdb={args.pdb.strip()}"

    try:
        import webview  # type: ignore
    except Exception:
        print(
            "Missing dependency: pywebview\n\n"
            "Install it with:\n"
            "  pip install pywebview\n\n"
            "Then re-run viewer.py.",
            file=sys.stderr,
        )
        return 2

    webview.create_window("3D Viewer (Mol*)", url, width=1100, height=750)
    webview.start()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
