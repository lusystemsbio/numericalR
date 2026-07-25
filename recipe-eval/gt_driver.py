#!/usr/bin/env python3
"""Run under the capture shim. Reads a JSON list of [section_code, python_code]
blocks (one chapter's Python, in order) and execs each in a shared namespace,
tagging its plots with the section via sitecustomize.mark(). Per-block failures
are logged, not fatal, so one bad chunk doesn't lose the rest of the chapter."""
import sys, json
import sitecustomize as cap

blocks = json.load(open(sys.argv[1]))
ns = {"__name__": "__main__"}
for sec, code in blocks:
    cap.mark(sec)
    try:
        exec(compile(code, f"<{sec}>", "exec"), ns)
    except Exception as e:
        print(f"BLOCKFAIL {sec}: {type(e).__name__}: {e}", file=sys.stderr)
