#!/usr/bin/env python3
"""Export the project's LyX document to PDF using an installed LyX.

Put this script beside 'riemanni hüpoteesist.lyx', then run:
    python3 build_pdf.py
    python3 build_pdf.py --watch

Requirements: Python 3.8+, LyX, and a working LaTeX installation.
No Python packages are needed. The .lyx source is never modified.
The PDF is replaced only after a successful export. Nothing is committed.

Other examples:
    python3 build_pdf.py "another document.lyx"
    python3 build_pdf.py --lyx "/path/to/lyx"
    python3 build_pdf.py --format pdf4   # XeLaTeX, if configured in LyX

Watch mode monitors the main .lyx file only, not images or included files.
Save the document in LyX before building or committing.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


DEFAULT_SOURCE = Path(__file__).resolve().parent / "riemanni hüpoteesist.lyx"


def build_pdf(source, output, lyx, export_format):
    """Stage a fresh export; leave an existing output untouched on failure."""
    output.parent.mkdir(parents=True, exist_ok=True)
    # Stage on the same filesystem so os.replace() can publish atomically.
    with tempfile.TemporaryDirectory(prefix=".lyx-export-", dir=output.parent) as tmp:
        staged = Path(tmp) / output.name
        command = [
            lyx, "--no-remote", "-batch", "--force-overwrite", "all",
            "--export-to", export_format, str(staged), str(source),
        ]
        env = os.environ.copy()
        if sys.platform.startswith("linux") and not (
            env.get("DISPLAY") or env.get("WAYLAND_DISPLAY")
        ):
            env.setdefault("QT_QPA_PLATFORM", "offscreen")
        print("Building: {}".format(source), flush=True)
        result = subprocess.run(command, cwd=source.parent, env=env)
        if result.returncode != 0:
            raise RuntimeError(
                "LyX export failed (exit {}). Previous PDF kept. "
                "See the LyX/LaTeX messages above.".format(result.returncode)
            )
        if not staged.is_file():
            raise RuntimeError("LyX produced no PDF. Previous PDF kept.")
        with staged.open("rb") as handle:
            if handle.read(5) != b"%PDF-":
                raise RuntimeError("LyX output is not a PDF. Previous PDF kept.")
        os.replace(staged, output)
    print("Updated: {}".format(output), flush=True)


def fingerprint(source):
    stat = source.stat()
    return stat.st_mtime_ns, stat.st_size


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, help="PDF path; default: beside the source")
    parser.add_argument("--lyx", default=os.environ.get("LYX", "lyx"),
                        help="LyX executable name or full path (or set LYX)")
    parser.add_argument("--format", default="pdf2", dest="export_format",
                        choices=("pdf", "pdf2", "pdf3", "pdf4", "pdf5"),
                        help="LyX export format; default pdf2 (pdfLaTeX)")
    parser.add_argument("--watch", action="store_true",
                        help="Build now, then rebuild after the main .lyx file changes")
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    output = (args.output.expanduser().resolve() if args.output
              else source.with_suffix(".pdf"))
    if not source.is_file() or source.suffix.lower() != ".lyx":
        parser.error("Source must be an existing .lyx file: {}".format(source))
    if output.suffix.lower() != ".pdf":
        parser.error("Output must have a .pdf extension")
    lyx = shutil.which(args.lyx)
    if lyx is None:
        parser.error("LyX was not found. Install LyX and LaTeX, then put LyX on PATH "
                     "or pass --lyx with its executable's full path.")

    if not args.watch:
        try:
            build_pdf(source, output, lyx, args.export_format)
        except (OSError, RuntimeError) as error:
            print("Error: {}".format(error), file=sys.stderr)
            return 1
        return 0

    print("Watching the main LyX file; press Ctrl+C to stop.", flush=True)
    previous = None
    while True:
        try:
            current = fingerprint(source)
            if current != previous:
                # A change made during export remains detectable on the next pass.
                previous = current
                # Wait briefly for an editor's save operation to finish.
                time.sleep(0.5)
                if fingerprint(source) != current:
                    previous = None
                    continue
                build_pdf(source, output, lyx, args.export_format)
        except (OSError, RuntimeError) as error:
            print("Error: {}".format(error), file=sys.stderr)
            # Failed exports are retried after the next source change.
        time.sleep(1)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        sys.exit(130)
