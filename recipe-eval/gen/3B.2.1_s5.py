import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -------------------- Model parameters (toggle switch) --------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# -------------------- Right-hand side and Jacobian --------------------
def rhs(X, Y):
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

def jac(X, Y):
    # u = (Y/Yth)^nY ; d/dY [1/(1+u)] = -nY*u / (Y*(1+u)^2)
    u = (Y / Yth) ** nY
    v = (X / Xth) ** nX
    a = -kX                                             # d(dX)/dX
    b = gX1 * (-nY * u / (Y * (1.0 + u) ** 2))          # d(dX)/dY
    c = gY1 * (-nX * v / (X * (1.0 + v) ** 2))          # d(dY)/dX
    d = -kY                                             # d(dY)/dY
    return np.array([[a, b], [c, d]])

# -------------------- Nullclines as explicit polylines --------------------
# X-nullcline: dX/dt = 0  ->  X = fx(Y), so parametrize by Y.
def fx(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX
# Y-nullcline: dY/dt = 0  ->  Y = fy(X), so parametrize by X.
def fy(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

N = 2000
Ygrid = np.linspace(0.0, 600.0, N)          # parameter for X-nullcline
Xgrid = np.linspace(0.0, 600.0, N)          # parameter for Y-nullcline
Xnull = np.column_stack([fx(Ygrid), Ygrid]) # points (X, Y) on X-nullcline
Ynull = np.column_stack([Xgrid, fy(Xgrid)]) # points (X, Y) on Y-nullcline

# -------------------- Segment-intersection primitive --------------------
def seg_intersect(p1, p2, p3, p4):
    # Solve p1 + t*(p2-p1) = p3 + s*(p4-p3) for t,s in [0,1].
    r = p2 - p1
    s_ = p4 - p3
    denom = r[0] * s_[1] - r[1] * s_[0]
    if abs(denom) < 1e-15:                  # parallel / degenerate
        return None
    q = p3 - p1
    t = (q[0] * s_[1] - q[1] * s_[0]) / denom
    s = (q[0] * r[1] - q[1] * r[0]) / denom
    if -1e-12 <= t <= 1 + 1e-12 and -1e-12 <= s <= 1 + 1e-12:
        return p1 + t * r
    return None

def newton_refine(pt):
    # Polish an approximate crossing to a true root of rhs using analytic Jacobian.
    X, Y = float(pt[0]), float(pt[1])
    for _ in range(50):
        F = rhs(X, Y)
        if np.linalg.norm(F) < 1e-12:
            break
        step = np.linalg.solve(jac(X, Y), F)
        X, Y = X - step[0], Y - step[1]
    return np.array([X, Y])

def dedupe(points, tol=1e-4):
    out = []
    for p in points:
        if not any(np.hypot(p[0] - q[0], p[1] - q[1]) < tol for q in out):
            out.append(p)
    return out

# -------------------- Exhaustive all-pairs intersection --------------------
def find_exhaustive(A, B):
    hits, pairs = [], 0
    for i in range(len(A) - 1):
        for j in range(len(B) - 1):
            pairs += 1                       # count every pair actually tested
            p = seg_intersect(A[i], A[i + 1], B[j], B[j + 1])
            if p is not None:
                hits.append(newton_refine(p))
    return dedupe(hits), pairs

# -------------------- Fast sign-change-narrowed intersection --------------------
def find_fast(A, B):
    # Residual along the Y-nullcline: horizontal gap to the X-nullcline at the same Y.
    # r(X) = X - fx(fy(X)); its sign changes exactly where the two curves cross.
    r = B[:, 0] - fx(B[:, 1])
    sign_change = np.where(r[:-1] * r[1:] <= 0.0)[0]   # candidate Y-null segments
    hits, pairs = [], 0
    for j in sign_change:
        # Only test X-null segments whose Y-range overlaps this candidate's Y-range.
        y_lo, y_hi = sorted([B[j, 1], B[j + 1, 1]])
        lo = max(np.searchsorted(Ygrid, y_lo) - 1, 0)
        hi = min(np.searchsorted(Ygrid, y_hi) + 1, len(A) - 1)
        for i in range(lo, hi):
            pairs += 1
            p = seg_intersect(A[i], A[i + 1], B[j], B[j + 1])
            if p is not None:
                hits.append(newton_refine(p))
    return dedupe(hits), pairs

# -------------------- Run both searches --------------------
ss_slow, pairs_slow = find_exhaustive(Xnull, Ynull)
ss_fast, pairs_fast = find_fast(Xnull, Ynull)
ss = sorted(ss_fast, key=lambda p: p[0])

# -------------------- Report steady states and stability --------------------
print(f"Number of steady states (exhaustive) : {len(ss_slow)}")
print(f"Number of steady states (fast)       : {len(ss_fast)}")

n_stable = n_unstable = 0
for idx, p in enumerate(ss, 1):
    ev = np.linalg.eigvals(jac(p[0], p[1]))
    stable = np.all(ev.real < 0)
    label = "stable" if stable else "unstable"
    n_stable += stable
    n_unstable += (not stable)
    print(f"Steady state {idx}: X = {p[0]:.6f}, Y = {p[1]:.6f}")
    print(f"   eigenvalues = {ev[0]:.6f}, {ev[1]:.6f}  ->  {label}")
    print(f"   residual |dX/dt,dY/dt| = {np.linalg.norm(rhs(p[0], p[1])):.3e}")

print(f"Stable steady states   : {n_stable}")
print(f"Unstable steady states : {n_unstable}")

# -------------------- Confirm agreement and speed-up --------------------
match = (len(ss_slow) == len(ss_fast) == 3)
for a in ss_slow:
    match = match and any(np.hypot(a[0]-b[0], a[1]-b[1]) < 1e-3 for b in ss_fast)
print(f"Fast and exhaustive agree on all steady states : {match}")
print(f"Segment pairs tested (exhaustive) : {pairs_slow}")
print(f"Segment pairs tested (fast)       : {pairs_fast}")
print(f"Speed-up factor (fewer pairs)     : {pairs_slow / max(pairs_fast, 1):.1f}x")

# -------------------- Plot nullclines and crossings --------------------
plt.figure(figsize=(7, 6))
plt.plot(Xnull[:, 0], Xnull[:, 1], label="X-nullcline (dX/dt=0)", color="tab:blue")
plt.plot(Ynull[:, 0], Ynull[:, 1], label="Y-nullcline (dY/dt=0)", color="tab:red")
for idx, p in enumerate(ss, 1):
    ev = np.linalg.eigvals(jac(p[0], p[1]))
    stable = np.all(ev.real < 0)
    plt.scatter(p[0], p[1], s=110, zorder=5,
                facecolor=("black" if stable else "white"),
                edgecolor="black",
                label=("stable steady state" if stable and idx == 1 else
                       "unstable steady state" if not stable else None))
    plt.annotate(f"({p[0]:.1f}, {p[1]:.1f})", (p[0], p[1]),
                 textcoords="offset points", xytext=(8, 8))
plt.xlabel("X"); plt.ylabel("Y")
plt.title("Toggle-switch nullclines and steady states")
plt.legend(); plt.xlim(0, 400); plt.ylim(0, 400); plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.2.1_s5.png")

# One-sentence explanation:
# The check confirms the result because the fast sign-change filter recovers the
# exact same three crossings as the brute-force all-pairs test while examining far
# fewer segment pairs, and the Jacobian eigenvalues show the expected two-stable/
# one-unstable bistable structure of a genuine toggle switch.
print("Explanation: agreement of the fast and exhaustive searches on the same three "
      "crossings, plus a Jacobian spectrum giving two stable and one unstable node, "
      "confirms these are the true bistable steady states of the toggle switch.")
