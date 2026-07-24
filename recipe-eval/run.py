#!/usr/bin/env python3
"""Execute every generated script in a clean sandbox and capture the result:
exit status, stdout, stderr tail, and whether a figure was written. This is the
'does it run + what did it output' layer; scoring reads run_results.json."""
import glob, json, os, subprocess

BASE = os.environ.get("RECIPE_EVAL_BASE", "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval")
GEN = os.path.join(BASE, "gen")
PY = "/usr/local/bin/python3"

def main():
    res = {}
    for f in sorted(glob.glob(os.path.join(GEN, "*.py"))):
        name = os.path.basename(f)[:-3]
        fig = os.path.join(GEN, name + ".png")
        had_fig = os.path.exists(fig)
        if had_fig: os.remove(fig)  # so we detect if THIS run makes it
        try:
            r = subprocess.run([PY, f], cwd=GEN, capture_output=True, text=True, timeout=180)
            res[name] = {"ran": r.returncode == 0, "exit": r.returncode,
                         "stdout": r.stdout[-4000:], "stderr": r.stderr[-1500:],
                         "fig": os.path.exists(fig)}
        except subprocess.TimeoutExpired:
            res[name] = {"ran": False, "exit": "timeout", "stdout": "", "stderr": "TIMEOUT", "fig": False}
        print(f"  {name}: ran={res[name]['ran']} exit={res[name]['exit']} fig={res[name]['fig']}")
    json.dump(res, open(os.path.join(BASE, "run_results.json"), "w"), indent=1)
    ran = sum(1 for v in res.values() if v["ran"])
    print(f"executed {len(res)} scripts; {ran} ran clean ({100*ran/max(len(res),1):.0f}%)")

if __name__ == "__main__":
    main()
