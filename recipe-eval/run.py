#!/usr/bin/env python3
"""Execute every generated script in a clean sandbox and capture the result:
exit status, stdout, stderr tail, and whether a figure was written. This is the
'does it run + what did it output' layer; scoring reads run_results.json.
Scripts run in parallel (each is its own process, so workers use many cores)."""
import glob, json, os, subprocess, concurrent.futures as cf

BASE = os.environ.get("RECIPE_EVAL_BASE", os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(BASE, "gen")
PY = os.environ.get("RECIPE_EVAL_PY", "/usr/local/bin/python3")
WORKERS = 6

def run_one(f):
    name = os.path.basename(f)[:-3]
    fig = os.path.join(GEN, name + ".png")
    if os.path.exists(fig): os.remove(fig)   # detect if THIS run writes it
    try:
        r = subprocess.run([PY, f], cwd=GEN, capture_output=True, text=True, timeout=180)
        v = {"ran": r.returncode == 0, "exit": r.returncode,
             "stdout": r.stdout[-4000:], "stderr": r.stderr[-1500:], "fig": os.path.exists(fig)}
    except subprocess.TimeoutExpired:
        v = {"ran": False, "exit": "timeout", "stdout": "", "stderr": "TIMEOUT", "fig": False}
    except Exception as e:
        v = {"ran": False, "exit": f"err:{e}", "stdout": "", "stderr": str(e), "fig": False}
    return name, v

def main():
    files = sorted(glob.glob(os.path.join(GEN, "*.py")))
    res = {}
    done = 0
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for name, v in ex.map(run_one, files):
            res[name] = v; done += 1
            if done % 25 == 0:
                print(f"  {done}/{len(files)} executed", flush=True)
    json.dump(res, open(os.path.join(BASE, "run_results.json"), "w"), indent=1)
    ran = sum(1 for v in res.values() if v["ran"])
    print(f"executed {len(res)} scripts; {ran} ran clean ({100*ran/max(len(res),1):.0f}%)")

if __name__ == "__main__":
    main()
