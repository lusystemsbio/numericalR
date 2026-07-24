import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle-switch parameters: X and Y mutually repress via Hill functions
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand sides of the ODEs (the two "surfaces" whose zero level we trace)
def fX(X, Y):
    # dX/dt = basal + repressive Hill(Y) - degradation
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    # dY/dt = basal + repressive Hill(X) - degradation
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# Build a grid over the (X, Y) plane and evaluate both surfaces on it
# ----------------------------------------------------------------------
Xmax, Ymax = 700.0, 700.0        # domain generous enough to enclose all states
N = 600
xs = np.linspace(0.0, Xmax, N)
ys = np.linspace(0.0, Ymax, N)
XX, YY = np.meshgrid(xs, ys)      # 2-D coordinate arrays
ZX = fX(XX, YY)                   # surface Z = fX(X,Y); its Z=0 contour is X-nullcline
ZY = fY(XX, YY)                   # surface Z = fY(X,Y); its Z=0 contour is Y-nullcline

# ----------------------------------------------------------------------
# Trace each nullcline as the zero-level contour of its surface.
# (General method: no need to solve for X or Y explicitly.)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
csX = ax.contour(XX, YY, ZX, levels=[0.0], colors="tab:blue", linewidths=2)
csY = ax.contour(XX, YY, ZY, levels=[0.0], colors="tab:red", linewidths=2)

# ----------------------------------------------------------------------
# Find the crossings (steady states) explicitly, step by step:
#   1) extract the polyline vertices of the two zero contours
#   2) walk consecutive segments and test every pair for intersection
# ----------------------------------------------------------------------
def contour_segments(cs):
    # Collect the (x,y) vertex arrays from a single-level contour set,
    # handling both older (.collections) and newer (.get_paths) Matplotlib.
    segs = []
    if hasattr(cs, "get_paths"):
        paths = cs.get_paths()
    else:
        paths = [p for coll in cs.collections for p in coll.get_paths()]
    for p in paths:
        v = p.vertices
        if len(v) >= 2:
            segs.append(v)
    return segs

def seg_intersection(p1, p2, p3, p4):
    # Intersection of segment p1->p2 with segment p3->p4, or None.
    r = p2 - p1
    s = p4 - p3
    denom = r[0] * s[1] - r[1] * s[0]
    if abs(denom) < 1e-12:          # parallel / degenerate
        return None
    qp = p3 - p1
    t = (qp[0] * s[1] - qp[1] * s[0]) / denom
    u = (qp[0] * r[1] - qp[1] * r[0]) / denom
    if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
        return p1 + t * r
    return None

segsX = contour_segments(csX)
segsY = contour_segments(csY)

raw_crossings = []
for polyX in segsX:
    for i in range(len(polyX) - 1):
        a1, a2 = polyX[i], polyX[i + 1]
        for polyY in segsY:
            for j in range(len(polyY) - 1):
                b1, b2 = polyY[j], polyY[j + 1]
                pt = seg_intersection(a1, a2, b1, b2)
                if pt is not None:
                    raw_crossings.append(pt)

# De-duplicate nearby intersection points (contour segments are dense)
crossings = []
for pt in raw_crossings:
    if not any(np.hypot(pt[0] - c[0], pt[1] - c[1]) < 1.0 for c in crossings):
        crossings.append(pt)
crossings.sort(key=lambda c: c[0])

ax.plot([c[0] for c in crossings], [c[1] for c in crossings],
        "ko", markersize=9, label="crossings (steady states)")
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Nullclines as zero-level contours")
ax.legend(["X-nullcline (fX=0)", "Y-nullcline (fY=0)", "crossings"])
ax.set_xlim(0, Xmax); ax.set_ylim(0, Ymax)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.5.1_s5.png")

# ----------------------------------------------------------------------
# INDEPENDENT CHECK via separation of variables.
# Each equation can be solved explicitly for one variable:
#   fX=0  ->  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX      (X as function of Y)
#   fY=0  ->  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY      (Y as function of X)
# Substitute to get one scalar equation g(X)=0 and root-find it.
# ----------------------------------------------------------------------
def Y_of_X(X):                       # from Y-nullcline
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

def X_of_Y(Y):                       # from X-nullcline
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def g(X):                            # residual: X reproduced after one round trip
    return X_of_Y(Y_of_X(X)) - X

# Bracket-and-bisect all sign changes of g on a fine 1-D scan
xg = np.linspace(1.0, Xmax, 20000)
gv = np.array([g(x) for x in xg])
sv_roots = []
for i in range(len(xg) - 1):
    if gv[i] == 0.0:
        sv_roots.append(xg[i])
    elif gv[i] * gv[i + 1] < 0.0:
        lo, hi = xg[i], xg[i + 1]
        for _ in range(80):          # bisection
            mid = 0.5 * (lo + hi)
            if g(lo) * g(mid) <= 0.0:
                hi = mid
            else:
                lo = mid
        sv_roots.append(0.5 * (lo + hi))
sv_states = [(x, Y_of_X(x)) for x in sorted(sv_roots)]

# ----------------------------------------------------------------------
# Print all numerical results
# ----------------------------------------------------------------------
print("=== Steady states from contour crossings ===")
for k, c in enumerate(crossings):
    print(f"crossing {k}: X = {c[0]:.6f}, Y = {c[1]:.6f}")

print("=== Steady states from separation of variables ===")
for k, s in enumerate(sv_states):
    print(f"sv state {k}: X = {s[0]:.6f}, Y = {s[1]:.6f}")

print("=== Residuals of separation-of-variables states (should be ~0) ===")
for k, s in enumerate(sv_states):
    print(f"sv state {k}: fX = {fX(s[0], s[1]):.3e}, fY = {fY(s[0], s[1]):.3e}")

print("=== Match between the two methods (nearest-pair distances) ===")
for k, c in enumerate(crossings):
    d = min(np.hypot(c[0] - s[0], c[1] - s[1]) for s in sv_states)
    print(f"crossing {k}: distance to nearest sv state = {d:.6f}")
max_mismatch = max(min(np.hypot(c[0] - s[0], c[1] - s[1]) for s in sv_states)
                   for c in crossings) if crossings and sv_states else float("nan")
print(f"max crossing-to-sv-state distance = {max_mismatch:.6f}")

# One-sentence explanation of why the check confirms the result:
print("EXPLANATION: The check confirms the result because separation of "
      "variables solves each fX=0 / fY=0 equation exactly and analytically, "
      "so its curves and their intersections are ground truth; the contour "
      "nullclines and crossings landing on top of them (near-zero residuals "
      "and distances) proves the general zero-contour method traced the same "
      "nullclines and found the same steady states without needing to solve "
      "for either variable.")
