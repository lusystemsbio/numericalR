#!/usr/bin/env python3
"""Re-run every generation under the capture shim (in the stable runner venv) and
save its plotted series to gen_cap/<name>.json. Offline; no API usage."""
import glob, json, os, subprocess, concurrent.futures as cf

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "gen")
CAPDIR = os.path.join(HERE, "gen_cap"); os.makedirs(CAPDIR, exist_ok=True)
CAPTURE = os.path.join(HERE, "capture")
PY = os.environ.get("RECIPE_EVAL_PY", "/usr/local/bin/python3")   # stable runner venv
WORKERS = 6

def one(f):
    name = os.path.basename(f)[:-3]
    out = os.path.join(CAPDIR, name + ".json")
    env = dict(os.environ, PYTHONPATH=CAPTURE, CAPTURE_OUT=out, MPLBACKEND="Agg")
    try:
        subprocess.run([PY, f], cwd=GEN, env=env, capture_output=True, text=True, timeout=180)
    except Exception:
        pass
    return name, os.path.exists(out) and len(json.load(open(out))) or 0

def main():
    files = sorted(glob.glob(os.path.join(GEN, "*.py")))
    done = 0
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for name, n in ex.map(one, files):
            done += 1
            if done % 50 == 0: print(f"  {done}/{len(files)}", flush=True)
    got = len(glob.glob(os.path.join(CAPDIR, "*.json")))
    print(f"captured plotted data for {got}/{len(files)} generations -> gen_cap/")

if __name__ == "__main__":
    main()
