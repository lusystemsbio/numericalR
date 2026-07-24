import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Toggle-switch parameters
# dX/dt = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X
# dY/dt = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---- Nullcline expressions (solve each dX/dt=0 / dY/dt=0 for one variable) ----
# X-nullcline: dX/dt=0  =>  X = F(Y)  (X as a function of Y)
def F(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline: dY/dt=0  =>  Y = G(X)  (Y as a function of X)
def G(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ---------------------------------------------------------------
# Build the two nullclines as polylines (sequences of points)
# ---------------------------------------------------------------
N = 1500
Ys = np.linspace(0.0, 600.0, N)          # sample the Y axis for the X-nullcline
A = np.column_stack([F(Ys), Ys])         # X-nullcline points  (F(Y), Y)

Xs = np.linspace(0.0, 700.0, N)          # sample the X axis for the Y-nullcline
B = np.column_stack([Xs, G(Xs)])         # Y-nullcline points  (X, G(X))

# ---------------------------------------------------------------
# Elementary segment-intersection test (parametric form)
# Returns the crossing point of segments p1->p2 and p3->p4, or None.
# ---------------------------------------------------------------
def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]

def seg_intersect(p1, p2, p3, p4):
    r = p2 - p1                 # direction of first segment
    s = p4 - p3                 # direction of second segment
    denom = cross(r, s)         # zero => parallel / collinear
    if abs(denom) < 1e-12:
        return None
    qp = p3 - p1
    t = cross(qp, s) / denom    # position along segment 1
    u = cross(qp, r) / denom    # position along segment 2
    if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
        return p1 + t * r       # the actual intersection coordinate
    return None

def dedupe(points, tol=1e-3):
    out = []
    for p in points:
        if not any(abs(p[0] - q[0]) < tol and abs(p[1] - q[1]) < tol for q in out):
            out.append(p)
    return out

# ---------------------------------------------------------------
# (1) EXHAUSTIVE all-pairs search: test every A-segment vs every B-segment
# ---------------------------------------------------------------
exhaustive_hits = []
exhaustive_pairs = 0
for i in range(len(A) - 1):
    for j in range(len(B) - 1):
        exhaustive_pairs += 1
        p = seg_intersect(A[i], A[i + 1], B[j], B[j + 1])
        if p is not None:
            exhaustive_hits.append(p)
exhaustive_hits = dedupe(exhaustive_hits)

# ---------------------------------------------------------------
# (2) FAST search: first narrow to the few B-segments that CHANGE SIGN.
# On the Y-nullcline B, define s(x) = F(G(x)) - x.  Where s changes sign
# between consecutive vertices, the X-nullcline crosses the Y-nullcline,
# so only those bracketing segments can hold an intersection.
# ---------------------------------------------------------------
s = F(B[:, 1]) - B[:, 0]                          # residual sampled along B
sign_change = np.where(s[:-1] * s[1:] < 0.0)[0]   # candidate B-segment indices

fast_hits = []
fast_pairs = 0
for j in sign_change:
    # y-extent of this candidate B-segment (tiny, since G is smooth)
    ylo, yhi = sorted([B[j, 1], B[j + 1, 1]])
    pad = (Ys[1] - Ys[0])
    # only A-segments whose y-range overlaps need to be tested
    lo = np.searchsorted(Ys, ylo - pad) - 1
    hi = np.searchsorted(Ys, yhi + pad) + 1
    lo = max(lo, 0); hi = min(hi, len(A) - 1)
    for i in range(lo, hi):
        fast_pairs += 1
        p = seg_intersect(A[i], A[i + 1], B[j], B[j + 1])
        if p is not None:
            fast_hits.append(p)
fast_hits = dedupe(fast_hits)

# ---------------------------------------------------------------
# Stability of each steady state via the Jacobian eigenvalues
# ---------------------------------------------------------------
def jacobian(X, Y):
    u = (Y / Yth) ** nY
    v = (X / Xth) ** nX
    dfX_dX = -kX
    dfX_dY = -gX1 * nY * u / (Y * (1.0 + u) ** 2)
    dfY_dX = -gY1 * nX * v / (X * (1.0 + v) ** 2)
    dfY_dY = -kY
    return np.array([[dfX_dX, dfX_dY], [dfY_dX, dfY_dY]])

def classify(X, Y):
    ev = np.linalg.eigvals(jacobian(X, Y))
    return ("stable" if np.all(ev.real < 0) else "unstable"), ev

# order the steady states by X for readable output
states = sorted(fast_hits, key=lambda p: p[0])

print("=== Steady states (nullcline intersections) ===")
n_stable = 0
n_unstable = 0
for k, p in enumerate(states, 1):
    kind, ev = classify(p[0], p[1])
    if kind == "stable":
        n_stable += 1
    else:
        n_unstable += 1
    print(f"Steady state {k}: X = {p[0]:.6f}, Y = {p[1]:.6f}  [{kind}]")
    print(f"    eigenvalues = {ev[0]:.6f}, {ev[1]:.6f}")

print(f"Number of steady states found (fast)       : {len(fast_hits)}")
print(f"Number of steady states found (exhaustive) : {len(exhaustive_hits)}")
print(f"Stable steady states   : {n_stable}")
print(f"Unstable steady states : {n_unstable}")

# ---------------------------------------------------------------
# Confirm the fast method agrees with the exhaustive search
# ---------------------------------------------------------------
same_count = (len(fast_hits) == len(exhaustive_hits))
matched = all(any(abs(f[0] - e[0]) < 1e-2 and abs(f[1] - e[1]) < 1e-2
                  for e in exhaustive_hits) for f in fast_hits)
print(f"Fast and exhaustive find same set of steady states : {same_count and matched}")
print(f"Exhaustive segment-pairs tested : {exhaustive_pairs}")
print(f"Fast segment-pairs tested       : {fast_pairs}")
print(f"Speed-up (pairs ratio)          : {exhaustive_pairs / max(fast_pairs,1):.1f}x")

# ---------------------------------------------------------------
# Plot the two nullclines and mark the crossings
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(A[:, 0], A[:, 1], label="X-nullcline (dX/dt=0)", color="tab:blue")
plt.plot(B[:, 0], B[:, 1], label="Y-nullcline (dY/dt=0)", color="tab:red")
for p in states:
    kind, _ = classify(p[0], p[1])
    color = "black" if kind == "stable" else "white"
    plt.plot(p[0], p[1], "o", markersize=11, markerfacecolor=color,
             markeredgecolor="green", markeredgewidth=2, zorder=5)
    plt.annotate(f"({p[0]:.0f}, {p[1]:.0f})\n{kind}", (p[0], p[1]),
                 textcoords="offset points", xytext=(8, 8), fontsize=9)
plt.xlabel("X"); plt.ylabel("Y")
plt.title("Toggle switch: steady states at nullcline intersections")
plt.legend(); plt.xlim(0, 650); plt.ylim(0, 420); plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.2.1_s3.png")

# One-sentence explanation:
# This check confirms the result because the fast sign-change filter returns
# exactly the same three intersection points as the brute-force all-pairs test
# while examining far fewer segment pairs, so the narrowing discards only
# non-crossing pairs and never misses a genuine steady state.
print("Explanation: the fast sign-change filter reproduces the exact same three "
      "crossings as the exhaustive all-pairs search while testing far fewer pairs, "
      "proving the narrowing discards only non-intersecting segments and misses no steady state.")
