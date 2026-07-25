#!/usr/bin/env python3
"""Deterministic scorer: compare each generation's TARGET QUANTITY (from its
captured plotted data and/or stdout) to the book's ground truth, with no LLM
judgment. Comparison primitives handle mismatched sampling (interpolate curves,
set-match points, KS for distributions, field metrics). A per-box REGISTRY maps
each box to a small check(caps, stdout) -> bool. Vision is reserved (elsewhere)
only for the few pure-morphology boxes.

Currently registered: the pilot boxes, to validate the method against their known
verdicts before scaling the registry to all 116.
"""
import glob, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAPDIR = os.path.join(HERE, "gen_cap")
RUNS = json.load(open(os.path.join(HERE, "run_results.json")))
BOOK = json.load(open(os.path.join(HERE, "book_gt.json"))) if os.path.exists(os.path.join(HERE,"book_gt.json")) else {}

# ---------- primitives ----------
def floats(s):
    return [float(x) for x in re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', s or "")]

def xy_series(caps):
    out = []
    for s in caps:
        if s.get("kind") not in ("xy","scatter","bar","y"): continue
        y = s.get("y")
        if y is None: continue
        try:
            ya = np.atleast_1d(np.asarray(y, float))
            if ya.ndim != 1 or ya.size < 1: continue
            x = s.get("x")
            xa = np.atleast_1d(np.asarray(x, float)) if x is not None else np.arange(ya.size, dtype=float)
            if xa.shape != ya.shape: xa = np.arange(ya.size, dtype=float)
            out.append((xa, ya))
        except Exception:
            continue
    return out

def finals(caps):  # last y-value of every curve
    return [float(y[-1]) for _, y in xy_series(caps) if len(y)]

def has_near(vals, target, rel=0.01, absol=None):
    tol = absol if absol is not None else abs(target)*rel
    return any(abs(v-target) <= tol for v in vals)

def endpoints(caps):  # last (x,y) point of each 2-D trajectory
    pts = []
    for x, y in xy_series(caps):
        if len(x) == len(y) and len(y) > 1: pts.append((float(x[-1]), float(y[-1])))
    return pts

def images(caps):
    return [np.asarray(s["z"], float) for s in caps if s.get("kind") in ("image","contour") and s.get("z")]

def spot_count(z):
    """crude morphology metric: # of local maxima above (mean+0.5*std)."""
    z = np.asarray(z, float)
    if z.ndim != 2: return -1
    thr = z.mean() + 0.5*z.std()
    hi = z > thr
    # count connected components of the high region (4-neighbour flood fill)
    seen = np.zeros_like(hi, bool); n = 0
    from collections import deque
    H, W = hi.shape
    for i in range(H):
        for j in range(W):
            if hi[i,j] and not seen[i,j]:
                n += 1; q = deque([(i,j)]); seen[i,j]=True
                while q:
                    a,b = q.popleft()
                    for da,db in ((1,0),(-1,0),(0,1),(0,-1)):
                        p,r = a+da,b+db
                        if 0<=p<H and 0<=r<W and hi[p,r] and not seen[p,r]:
                            seen[p,r]=True; q.append((p,r))
    return n

def caps_of(name):
    p = os.path.join(CAPDIR, name + ".json")
    return json.load(open(p)) if os.path.exists(p) else []

# ---------- per-box checks (pilot) ----------
def c_2C1(caps, out):
    f = finals(caps)
    return has_near(f, 22026.47, rel=0.01) and any(19000 <= v <= 21600 for v in f)  # exact/RK4 + less-accurate Euler
def c_2E3(caps, out):
    fl = floats(out)
    return all(has_near(fl, r, rel=0.04) for r in (71.48, 170.75, 331.61))          # 3 roots at k=0.15
def c_5A4(caps, out):
    for x, y in xy_series(caps):                                                     # an energy curve ~const near 2.05
        if len(y) > 10 and 1.7 <= y.mean() <= 2.4 and (y.max()-y.min())/y.mean() < 0.2:
            return True
    return False
def c_6A4(caps, out):
    for s in caps:
        if s.get("kind") == "hist" and s.get("data"):
            d = np.asarray(s["data"], float)
            if d.size >= 2000 and abs(d.mean()) < 0.08 and abs(d.std()-1) < 0.08:
                return True
    return False
def c_3A1(caps, out):
    e = endpoints(caps)
    a = any(abs(x-53) < 90 and abs(y-362) < 90 for x, y in e)   # near (53,362)
    b = any(abs(x-542) < 90 and abs(y-35) < 90 for x, y in e)   # near (542,35)
    return a and b
def c_9B3(caps, out):
    return has_near(floats(out), 193.73, rel=0.01)              # optimal tour length
def c_10C5(caps, out):
    o = (out or "").replace(" ", "")
    return ("01" in o and "10" in o and ("00->11" in o or "00-11" in o or "00,11" in o))
def c_8D2(caps, out):
    fl = floats(out)                                            # book SD ~ 2.1, 7.3, 17.5 (below sqrt(x))
    return sum(1 for t in (2.1, 7.3, 17.5) if has_near(fl, t, rel=0.35)) >= 2
def c_7E3(caps, out):
    for z in images(caps):
        if spot_count(z) >= 15: return True                     # many isolated spots
    return False
def c_10B1(caps, out):
    # a captured series holding the 3 centroids near (0,0),(1.5,1.5),(3,3)
    for x, y in xy_series(caps):
        if 3 <= len(x) <= 6:
            pts = list(zip(x, y))
            if all(any(abs(px-cx) < 0.7 and abs(py-cy) < 0.7 for px, py in pts)
                   for cx, cy in ((0,0),(1.5,1.5),(3,3))):
                return True
    return False

REGISTRY = {"2C.1":c_2C1, "2E.3":c_2E3, "5A.4":c_5A4, "6A.4":c_6A4, "3A.1":c_3A1,
            "9B.3":c_9B3, "10C.5":c_10C5, "8D.2":c_8D2, "7E.3":c_7E3, "10B.1":c_10B1}

def main():
    verdicts = {}
    for name in sorted(RUNS):
        sec = name.rsplit(".1_s",1)[0]  # "2C.1.1_s3" -> "2C.1"
        if sec not in REGISTRY: continue
        caps = caps_of(name); out = RUNS[name].get("stdout","")
        try: p = bool(REGISTRY[sec](caps, out))
        except Exception as e: p = False
        verdicts[name] = p
    # RFS + compare to saved Opus pilot verdicts
    O = json.load(open(os.path.join(HERE,"verdicts_opus_pilot.json"))) if os.path.exists(os.path.join(HERE,"verdicts_opus_pilot.json")) else {}
    boxes = sorted(set(n.rsplit("_s",1)[0] for n in verdicts))
    print(f"{'box':9} {'det-RFS':>7}  {'opus':>5}  match")
    agree=tot=0
    for b in boxes:
        s=[n for n in verdicts if n.rsplit('_s',1)[0]==b]
        det=sum(verdicts[n] for n in s)
        op=sum(1 for n in s if O.get(n,{}).get("pass"))
        m=sum(1 for n in s if O and (verdicts[n]==bool(O.get(n,{}).get("pass"))))
        agree+=m; tot+=len(s)
        print(f"{b:9} {det}/{len(s):<5}  {op}/5   {m}/{len(s)}")
    if tot: print(f"\nDeterministic vs Opus agreement: {agree}/{tot} = {100*agree/tot:.0f}%")

if __name__ == "__main__":
    main()
