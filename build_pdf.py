#!/usr/bin/env python3
"""Export the project's LyX document to PDF using an installed LyX.

Put this script beside 'riemanni hüpoteesist.lyx', then run:
    python3 build_pdf.py
    python3 build_pdf.py --watch

Requirements: Python 3.8+, LyX, and a working LaTeX installation.
No Python packages are needed. The .lyx source is never modified.
The PDF is replaced only after a successful export. Nothing is committed.
If a pdfLaTeX export fails, diagnostic files are saved beside the PDF and
the first LaTeX error is printed. These files do not replace the PDF.

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


def diagnose_pdflatex(source, output, lyx, env):
    """Export TeX separately and retain a noninteractive compiler log."""
    debug = Path(tempfile.mkdtemp(prefix="lyx-diagnostics-", dir=output.parent))
    print("Saving diagnostics: {}".format(debug), flush=True)
    tex = debug / "document.tex"
    export_log = debug / "lyx-export.log"
    with export_log.open("wb") as log:
        result = subprocess.run(
            [lyx, "--no-remote", "-batch", "--force-overwrite", "all",
             "--export-to", "pdflatex", str(tex), str(source)],
            cwd=source.parent, env=env, stdin=subprocess.DEVNULL,
            stdout=log, stderr=subprocess.STDOUT,
        )
    if result.returncode != 0 or not tex.is_file():
        print("Could not export diagnostic TeX; see {}".format(export_log), flush=True)
        print(export_log.read_text(encoding="utf-8", errors="replace")[-6000:])
        return
    engine = shutil.which("pdflatex")
    if engine is None:
        print("pdflatex is not on PATH. Diagnostic TeX saved at {}".format(tex))
        return
    tex_env = env.copy()
    # Keep source-relative inputs available, and preserve TeX's default paths.
    tex_env["TEXINPUTS"] = str(source.parent) + os.pathsep + env.get("TEXINPUTS", "")
    console_log = debug / "pdflatex-console.log"
    with console_log.open("wb") as log:
        subprocess.run(
            [engine, "-interaction=nonstopmode", "-halt-on-error", tex.name],
            cwd=debug, env=tex_env, stdin=subprocess.DEVNULL,
            stdout=log, stderr=subprocess.STDOUT,
        )
    lines = console_log.read_text(encoding="utf-8", errors="replace").splitlines()
    first_error = next((i for i, line in enumerate(lines) if line.startswith("!")), None)
    if first_error is not None:
        print("\nFirst LaTeX error:\n" + "\n".join(lines[first_error:first_error + 15]))
    else:
        print("\nDiagnostic compiler output:\n" + "\n".join(lines[-25:]))
    print("Full diagnostic output: {}".format(console_log), flush=True)
    if (debug / "document.log").is_file():
        print("Full LaTeX log: {}".format(debug / "document.log"), flush=True)


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
            if export_format == "pdf2":
                try:
                    diagnose_pdflatex(source, output, lyx, env)
                except OSError as error:
                    print("Could not complete diagnostics: {}".format(error), file=sys.stderr)
            raise RuntimeError(
                "LyX export failed (exit {}). Previous PDF kept. "
                "See the diagnostic messages above.".format(result.returncode)
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