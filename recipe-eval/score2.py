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

# ---------- more primitives ----------
def nf(caps, target, rel=0.03, n=1):           # >=n curves ending near target
    return sum(1 for f in finals(caps) if abs(f-target) <= abs(target)*rel) >= n
def ally(caps):                                # every plotted y value
    v = []
    for _, y in xy_series(caps): v += list(y)
    return v
def yspan_covers(caps, lo, hi):                # plotted y reaches both a low and a high band
    v = ally(caps)
    return bool(v) and min(v) < lo and max(v) > hi
def local_minima(y):
    y = np.asarray(y, float); return int(((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:])).sum())

# ---------- Part 2 checks ----------
def c_2A1(caps, out):  return nf(caps, 500, 0.03, 3)                 # ICs converge to g/k=500
def c_2B1(caps, out):  return nf(caps, 500, 0.03, 1)                 # Euler reaches 500
def c_2B4(caps, out):  return nf(caps, 500, 0.02, 1)                 # Heun vs exact -> 500
def c_2B5(caps, out):  return nf(caps, 500, 0.02, 1)                 # RK2 -> 500
def c_2B6(caps, out):  return nf(caps, 500, 0.02, 1)                 # RK4 -> 500
def c_2B8(caps, out):  return nf(caps, 5, 0.1, 1)                    # backward Euler stable -> g/k=5
def c_2C2(caps, out):  return nf(caps, 100, 0.05, 1)                 # logistic -> B=100
def c_2D1(caps, out):  return has_near(floats(out), 250, rel=0.08)  # steady state ~250 (stable)
def c_2D3(caps, out):                                               # potential has two wells (k=0.15)
    return any(local_minima(y) >= 2 for _, y in xy_series(caps) if len(y) > 20)
def c_2E4(caps, out):  return all(has_near(floats(out), r, rel=0.05) for r in (71.5,170.8,331.6))
def c_scurve(caps, out): return yspan_covers(caps, 120, 320)        # bifurcation S-curve spans both branches
# ---------- primitives for Parts 3-10 ----------
def pts(caps):
    P = []
    for x, y in xy_series(caps):
        if len(x) == len(y): P += list(zip(x.tolist(), y.tolist()))
    return P
def produced(c, o):            # setup boxes: at least ran and emitted a plot or output
    return bool(c) or len((o or '').strip()) > 30
def reaches(caps, targets, tol):               # endpoints reach every target state
    e = endpoints(caps)
    return all(any(abs(px-tx) < tol and abs(py-ty) < tol for px, py in e) for tx, ty in targets)
def near_pt(caps, tx, ty, tol):                # any plotted point near (tx,ty)
    return any(abs(px-tx) < tol and abs(py-ty) < tol for px, py in pts(caps))
def hv(out, t, rel=0.05, absol=None):          # stdout has a value near t
    return has_near(floats(out), t, rel=rel, absol=absol)
def local_maxima(y):
    y = np.asarray(y, float)
    return int(((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])).sum()) if y.size > 2 else 0
def osc(caps, k=3):                             # some curve oscillates (>=k local maxima)
    return any(local_maxima(y) >= k for _, y in xy_series(caps) if len(y) > 20)
def periodic_peaks(caps, k=4):                  # 1D Turing profile: regular spatial peaks
    for _, y in xy_series(caps):
        if len(y) > 20 and local_maxima(y) >= k and np.ptp(y) > 0.3*(abs(np.mean(y))+1e-9):
            return True
    return False
def c_7C2(c,o): return periodic_peaks(c)
def c_7C3(c,o): return periodic_peaks(c)
def sustained_osc(caps, k=3):                   # oscillates AND is still oscillating at the end
    for _, y in xy_series(caps):
        if len(y) > 30 and local_maxima(y) >= k:
            y = np.asarray(y, float); n = len(y)
            overall = np.ptp(y)
            if overall > 1e-6 and np.ptp(y[-n//4:]) > 0.3*overall:   # last quarter still swings
                return True
    return False
def n_attractors(caps, tol):
    reps = []
    for p in endpoints(caps):
        if not any(abs(p[0]-r[0]) < tol and abs(p[1]-r[1]) < tol for r in reps): reps.append(p)
    return len(reps)
def hist_ok(caps, mean, std, mtol=0.1, stol=0.15, nmin=2000):
    for s in caps:
        if s.get("kind") == "hist" and s.get("data"):
            d = np.asarray(s["data"], float)
            if d.size >= nmin and abs(d.mean()-mean) < mtol and abs(d.std()-std) < stol*max(std,1):
                return True
    return False
def curve_max_at(caps, x0, xtol, ymin_below=None):   # a curve peaks near x=x0
    for x, y in xy_series(caps):
        if len(y) > 5 and abs(float(x[np.argmax(y)]) - x0) < xtol:
            return True
    return False

TOG = [(53,362),(542,35)]   # toggle stable steady states; saddle ~ (190,127)
def T(t): return lambda c,o: t                 # constant (lenient placeholder)

# ---------- Part 3 ----------
def c_3A4(c,o): return all(near_pt(c,*s,tol=45) for s in [(53,362),(190,127),(542,35)])  # nullclines through 3 states
def c_3B2(c,o): return reaches(c,TOG,90) or all(near_pt(c,*s,25) for s in TOG)           # find the states
def c_3B3(c,o): return ("stable" in o.lower()) and (("saddle" in o.lower()) or ("unstable" in o.lower()))
def c_3C2(c,o): return reaches(c,[(0,5),(8,1)],1.2)                                       # washout & coexistence
def c_3D2(c,o): return osc(c,2)                                                            # LV oscillations/closed orbits
def c_3D6(c,o): return len([1 for v in floats(o) if 0.005<abs(v)<3]) >= 3                 # fitted params printed
def c_3D7(c,o): return len([1 for v in floats(o) if 0.005<abs(v)<3]) >= 3
def c_3E1(c,o): return yspan_covers(c,100,400)                                             # bifurcation branches
def c_3F2(c,o): return near_pt(c,190,127,60)                                               # separatrix through saddle
def c_3G2(c,o): return any(local_minima(y)>=1 for _,y in xy_series(c) if len(y)>20)        # potential minima
def c_3H1(c,o): return len(xy_series(c)) >= 1
def c_3H2(c,o): return n_attractors(c,0.4)==1                                              # single stable state
def c_3H4(c,o): return n_attractors(c,0.4)>=2                                              # bistable
def c_3H5(c,o): return sustained_osc(c,3)                                                            # sustained oscillation
def c_3H6(c,o):
    e=[y for _,y in xy_series(c)]; f=[float(y[-1]) for y in e if len(y)]
    return any(v>0.3 for v in f) and any(abs(v)<0.05 for v in f)                          # some coexist, some ~0

# ---------- Part 4 ----------
def c_4A1(c,o): return osc(c,2)                              # r=-1.7 growing oscillation
def c_4A2(c,o): return osc(c,2)
def c_4B1(c,o): return sustained_osc(c,3)                              # sustained limit cycle at r=1.7
def c_4B2(c,o): return sustained_osc(c,3)
def c_4B3(c,o): return len(xy_series(c))>=1                  # delayed LV orbits
def c_4C1(c,o): return not sustained_osc(c,3)          # relaxes, no sustained oscillation
def c_4C2(c,o): return sustained_osc(c,3)                              # delay -> sustained oscillation
def c_4C3(c,o): return sustained_osc(c,3)                              # rings oscillate

# ---------- Part 5 ----------
def c_5A1(c,o): return osc(c,3)                              # sinusoid
def c_5A2(c,o):                                             # Euler energy GROWS (secular)
    for x,y in xy_series(c):
        if len(y)>20 and 1.5<y.mean()<50 and y[-len(y)//5:].mean() > 1.15*y[:len(y)//5].mean():
            return True
    return False
def c_5A3(c,o): return c_5A4(c,o)                            # leapfrog energy bounded
def c_5A5(c,o): return osc(c,3)                              # Verlet reproduces oscillation
def c_5B1(c,o): return near_pt(c,4,0,1) or len(pts(c))>50    # orbit in plane
def c_5B2(c,o): return len(xy_series(c))>=1
def c_5B3(c,o): return len(xy_series(c))>=1
def c_5C4(c,o): return produced(c,o)                                 # setup (lenient; audited by Opus)
def c_5C5(c,o): return any(s.get("kind")=="scatter" for s in c)   # particle positions
def c_5C7(c,o): return curve_max_at(c,1.0,0.35)             # g(r) first peak near r=1

# ---------- Part 6 ----------
def c_6A3(c,o):                                             # exponential samples
    for s in c:
        if s.get("kind")=="hist" and s.get("data"):
            d=np.asarray(s["data"],float)
            if d.size>=2000 and (d>=0).mean()>0.98 and 0.7<d.mean()<1.4: return True
    return False
def c_6B1(c,o): return len(xy_series(c))>=3                 # ensemble of walks
def c_6B3(c,o): return len(xy_series(c))>=1
def c_6B4(c,o): return len(pts(c))>50                       # 2D path
def c_6B5(c,o): return len(xy_series(c))>=1
def c_6C1(c,o): return len(xy_series(c))>=3                 # E-M brownian trajectories
def c_6C2(c,o): return len(xy_series(c))>=1                 # OU
def c_6C3(c,o): return produced(c,o)                                 # setup
def c_6C4(c,o):                                             # bistable: visits ~100 and ~300
    for x,y in xy_series(c):
        if len(y)>50 and (y<160).mean()>0.05 and (y>240).mean()>0.05: return True
    return False
def c_6D1(c,o):                                            # toggle hops both states
    for x,y in xy_series(c):
        if len(y)>50 and y.min()<150 and y.max()>250: return True
    return osc(c,2)
def c_6D2(c,o): return len(floats(o))>=1
def c_6D3(c,o): return len(floats(o))>=1

# ---------- Part 7 (non-morphology) ----------
def c_7A1(c,o): return produced(c,o)
def c_7A2(c,o): return len(xy_series(c))>=1                 # spreading distribution
def c_7B1(c,o): return len(xy_series(c))>=2                 # traveling fronts
def c_7B2(c,o): return len(xy_series(c))>=1                 # relaxes to Gaussian
def c_7C1(c,o): return produced(c,o)
def c_7E1(c,o): return produced(c,o)

# ---------- Part 8 ----------
def c_8A1(c,o): return hv(o,3.14159,rel=0.05)
def c_8A2(c,o): return hv(o,3.14159,rel=0.06)
def c_8A3(c,o): return hv(o,0.8862,rel=0.05)               # sqrt(pi)/2
def c_8A4(c,o): return hv(o,0.8862,rel=0.05)
def c_8B1(c,o): return hv(o,0.25,absol=0.06) or hv(o,0.75,absol=0.06)
def c_8B3(c,o): return hv(o,0.25,absol=0.06) or hv(o,0.75,absol=0.06)
def c_8B4(c,o): return hist_ok(c,0.0,0.707,mtol=0.15,stol=0.2) or hist_ok(c,0,1,0.15,0.3)
def c_8B5(c,o): return hv(o,-9,absol=0.5) or any(y.min()<=-8.5 for _,y in xy_series(c) if len(y))
def c_8C3(c,o):                                            # energy falls then plateaus
    for x,y in xy_series(c):
        if len(y)>20 and y[:len(y)//5].mean() > y[-len(y)//5:].mean()+abs(y[-len(y)//5:].mean())*0.1: return True
    return False
def c_8D1(c,o): return produced(c,o)
def c_8D3(c,o): return hv(o,160,rel=0.1)                    # mean stays 160
def c_8D4(c,o): return len(floats(o))>=2                   # SD comparison printed

# ---------- Part 9 ----------
def c_9A2(c,o): return any(near_pt(c,*m,0.6) for m in [(3,2),(-2.805,3.131),(-3.779,-3.283),(3.584,-1.848)])
def c_9A3(c,o): return hv(o,0,absol=1.0) or any(near_pt(c,*m,0.6) for m in [(3,2),(-2.805,3.131),(-3.779,-3.283),(3.584,-1.848)])
def c_9A4(c,o): return c_9A3(c,o)
def c_9A5(c,o): return c_9A3(c,o)
def c_9B1(c,o): return any(v>=6 for v in floats(o))        # alignment score printed (positive)
def c_9B2(c,o): return any(v>=6 for v in floats(o))
def c_9C1(c,o): return hv(o,0,absol=1.0)                   # Rastrigin best -> ~0
def c_9C2(c,o): return len(floats(o))>=1                   # tour length printed

# ---------- Part 10 ----------
def c_10A2(c,o): return produced(c,o)                               # PCA (audited)
def c_10A4(c,o): return len(pts(c))>10
def c_10A5(c,o): return len(pts(c))>10
def c_10B2(c,o): return produced(c,o)                               # HCA -> 3 clusters
def c_10B3(c,o): return any(s.get("kind")=="image" for s in c)   # gene heatmap
def c_10B4(c,o): return hv(o,3,absol=0.5)                  # GMM picks 3 components
def c_10C1(c,o): return hv(o,4.6,rel=0.15) or hv(o,5,absol=0.5) or hv(o,0.57,rel=0.2)
def c_10C2(c,o): return len(xy_series(c))>=1 or any(s.get("kind")=="hist" for s in c)
def c_10C3(c,o): return hv(o,0.37,rel=0.2) or hv(o,18,absol=1) or hv(o,16,absol=1)
def c_10C4(c,o): return produced(c,o)
def c_10C6(c,o):
    x=(o or "").replace(" ","")
    return "000" in x and "111" in x

def c_1A3(c,o): return "even" in (o or "").lower() and "odd" in (o or "").lower()
def c_1A5(c,o): return any(s.get("kind")=="image" for s in c) or len(pts(c))>4
def c_1C3(c,o): return hv(o, 1.41421, rel=0.002)
REGISTRY = {"1A.3":c_1A3,"1A.5":c_1A5,"1C.3":c_1C3,
            "2C.1":c_2C1, "2E.3":c_2E3, "5A.4":c_5A4, "6A.4":c_6A4, "3A.1":c_3A1,
            "9B.3":c_9B3, "10C.5":c_10C5, "8D.2":c_8D2, "7E.3":c_7E3, "10B.1":c_10B1,
            "2A.1":c_2A1, "2B.1":c_2B1, "2B.4":c_2B4, "2B.5":c_2B5, "2B.6":c_2B6,
            "2B.8":c_2B8, "2C.2":c_2C2, "2D.1":c_2D1, "2D.3":c_2D3, "2E.4":c_2E4,
            "2E.2":c_scurve, "2E.5":c_scurve, "2F.2":c_scurve, "2F.6":c_scurve,
            "3A.4":c_3A4,"3A.5":c_3A4,"3A.6":c_3A4,"3B.2":c_3B2,"3B.3":c_3B3,"3C.2":c_3C2,
            "3D.2":c_3D2,"3D.6":c_3D6,"3D.7":c_3D7,"3E.1":c_3E1,"3F.2":c_3F2,"3G.2":c_3G2,
            "3H.1":c_3H1,"3H.2":c_3H2,"3H.4":c_3H4,"3H.5":c_3H5,"3H.6":c_3H6,
            "4A.1":c_4A1,"4A.2":c_4A2,"4B.1":c_4B1,"4B.2":c_4B2,"4B.3":c_4B3,
            "4C.1":c_4C1,"4C.2":c_4C2,"4C.3":c_4C3,
            "5A.1":c_5A1,"5A.2":c_5A2,"5A.3":c_5A3,"5A.5":c_5A5,"5B.1":c_5B1,"5B.2":c_5B2,
            "5B.3":c_5B3,"5C.4":c_5C4,"5C.5":c_5C5,"5C.7":c_5C7,
            "6A.3":c_6A3,"6B.1":c_6B1,"6B.3":c_6B3,"6B.4":c_6B4,"6B.5":c_6B5,
            "6C.1":c_6C1,"6C.2":c_6C2,"6C.3":c_6C3,"6C.4":c_6C4,"6D.1":c_6D1,"6D.2":c_6D2,"6D.3":c_6D3,
            "7A.1":c_7A1,"7A.2":c_7A2,"7B.1":c_7B1,"7B.2":c_7B2,"7C.1":c_7C1,"7E.1":c_7E1,
            "8A.1":c_8A1,"8A.2":c_8A2,"8A.3":c_8A3,"8A.4":c_8A4,"8B.1":c_8B1,"8B.3":c_8B3,
            "8B.4":c_8B4,"8B.5":c_8B5,"8C.3":c_8C3,"8D.1":c_8D1,"8D.3":c_8D3,"8D.4":c_8D4,
            "9A.2":c_9A2,"9A.3":c_9A3,"9A.4":c_9A4,"9A.5":c_9A5,"9B.1":c_9B1,"9B.2":c_9B2,
            "9C.1":c_9C1,"9C.2":c_9C2,
            "10A.2":c_10A2,"10A.4":c_10A4,"10A.5":c_10A5,"10B.2":c_10B2,"10B.3":c_10B3,
            "10B.4":c_10B4,"10C.1":c_10C1,"10C.2":c_10C2,"10C.3":c_10C3,"10C.4":c_10C4,"10C.6":c_10C6}
REGISTRY["7C.2"]=c_7C2; REGISTRY["7C.3"]=c_7C3
MORPHOLOGY = {"7E.3", "7E.2", "7D.2"}   # true 2D morphology -> vision
# Vision verdicts (classified by eye from the captured 2D fields):
VISION = {"7E.2":[1,1,1,1,1],   # stripes/labyrinth (correct)
          "7E.3":[0,0,0,0,0],   # should be spots; all labyrinth (genuine failure)
          "7D.2":[1,1,1,0,0]}   # spiral wave: s1-3 yes, s4-5 just a gradient

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
        has_ref = any(n in O for n in s)
        op = f"{sum(1 for n in s if O.get(n,{}).get('pass'))}/5" if has_ref else "  -"
        if has_ref:
            m=sum(1 for n in s if verdicts[n]==bool(O.get(n,{}).get("pass")))
            agree+=m; tot+=len(s); mstr=f"{m}/{len(s)}"
        else: mstr="  -"
        print(f"{b:9} det {det}/{len(s):<3} opus {op:>4} match {mstr}")
    if tot: print(f"\nDeterministic vs Opus agreement (pilot boxes only): {agree}/{tot} = {100*agree/tot:.0f}%")

if __name__ == "__main__":
    main()
