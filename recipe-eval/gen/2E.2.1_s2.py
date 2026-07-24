import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---- Model parameters ----
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

# f(X,k): basal + excitatory Hill - linear degradation
def f(X, k):
    h = (X / Xth) ** n
    return g0 + g1 * h / (1.0 + h) - k * X

# df/dX analytically: derivative of Hill term minus k
def dfdX(X, k):
    r = X / Xth
    # d/dX [ h/(1+h) ] with h=(X/Xth)^n  ->  n*r^(n-1)/Xth / (1+h)^2
    h = r ** n
    return g1 * (n * r ** (n - 1) / Xth) / (1.0 + h) ** 2 - k

# ---- Root finding at a given k: bracket-and-solve on a fine X grid ----
def find_roots(k, Xmax, npts=4000):
    Xg = np.linspace(1e-6, Xmax, npts)   # avoid X=0 singularity in r^(n-1)
    fg = f(Xg, k)
    roots = []
    for i in range(len(Xg) - 1):
        # sign change => a root lies in [Xg[i], Xg[i+1]]; solve with library root-finder
        if fg[i] == 0.0:
            roots.append(Xg[i])
        elif fg[i] * fg[i + 1] < 0.0:
            roots.append(brentq(f, Xg[i], Xg[i + 1], args=(k,)))
    # dedupe roots that are numerically identical
    uniq = []
    for r in roots:
        if not any(abs(r - u) < 1e-6 for u in uniq):
            uniq.append(r)
    return uniq

# ---- Sweep k on a grid; collect (k, X, stable?) points ----
Xmax = g0 / min(0.05, 1.0) + g1 / 0.05 + 10 * Xth  # generous upper bound on X
ks = np.linspace(0.05, 0.35, 400)

pts_k, pts_X, pts_stable = [], [], []
for k in ks:
    for X in find_roots(k, Xmax):
        slope = dfdX(X, k)
        stable = slope < 0.0          # steady state stable when df/dX < 0
        pts_k.append(k)
        pts_X.append(X)
        pts_stable.append(stable)

pts_k = np.array(pts_k)
pts_X = np.array(pts_X)
pts_stable = np.array(pts_stable)

# ---- Order the scattered points with a nearest-neighbor walk ----
# Work in a normalized (k, X) space so both axes contribute comparably.
kn = (pts_k - pts_k.min()) / (pts_k.max() - pts_k.min())
Xn = (pts_X - pts_X.min()) / (pts_X.max() - pts_X.min())
coords = np.column_stack([kn, Xn])

N = len(coords)
visited = np.zeros(N, dtype=bool)
order = [0]                 # start from the first point
visited[0] = True
for _ in range(N - 1):
    cur = coords[order[-1]]
    d = np.sum((coords - cur) ** 2, axis=1)
    d[visited] = np.inf     # ignore already-visited points
    nxt = int(np.argmin(d))
    order.append(nxt)
    visited[nxt] = True
order = np.array(order)

ok = pts_k[order]
oX = pts_X[order]
os = pts_stable[order]

# ---- Check: count steady states across k to confirm the S-shape ----
counts = {}
for k in ks:
    counts.setdefault(len(find_roots(k, Xmax)), 0)
    counts[len(find_roots(k, Xmax))] += 1
three_ks = [k for k in ks if len(find_roots(k, Xmax)) == 3]

print("Model parameters: g0=%.1f, g1=%.1f, Xth=%.1f, n=%.1f" % (g0, g1, Xth, n))
print("k grid: from %.4f to %.4f (%d values)" % (ks.min(), ks.max(), len(ks)))
print("Total steady-state points found: %d" % N)
print("Stable points: %d" % int(np.sum(pts_stable)))
print("Unstable points: %d" % int(np.sum(~pts_stable)))
print("Distribution of (#steady states -> #k values):")
for c in sorted(counts):
    print("  %d steady states: %d k-values" % (c, counts[c]))
if three_ks:
    print("Bistable (3 steady states) k range: %.4f to %.4f" % (min(three_ks), max(three_ks)))
    kmid = three_ks[len(three_ks) // 2]
    rts = sorted(find_roots(kmid, Xmax))
    print("Example at k=%.4f -> %d steady states: %s" % (kmid, len(rts), ", ".join("%.3f" % r for r in rts)))
    for r in rts:
        print("    X=%.3f  df/dX=%+.5f  -> %s" % (r, dfdX(r, kmid), "STABLE" if dfdX(r, kmid) < 0 else "UNSTABLE"))
n_stable_in_three = sum(1 for r in find_roots(three_ks[len(three_ks)//2], Xmax) if dfdX(r, three_ks[len(three_ks)//2]) < 0)
print("In the 3-state region: %d stable + %d unstable steady states" %
      (n_stable_in_three, 3 - n_stable_in_three))
print("Check explanation: a middle k-band with three steady states (two stable, one unstable in "
      "between) that collapses to one outside it is exactly the fold/saddle-node signature of an "
      "S-shaped bistable curve, so observing that count pattern confirms the result.")

# ---- Plot: steady-state X vs k, colored by stability ----
plt.figure(figsize=(8, 6))
# faint connecting line following the nearest-neighbor ordering (the traced curve)
plt.plot(ok, oX, color="0.75", lw=0.8, zorder=1)
plt.scatter(pts_k[pts_stable], pts_X[pts_stable], c="tab:blue", s=12,
            label="stable (df/dX < 0)", zorder=2)
plt.scatter(pts_k[~pts_stable], pts_X[~pts_stable], c="tab:red", s=12,
            label="unstable (df/dX > 0)", zorder=2)
plt.xlabel("k (control parameter)")
plt.ylabel("steady-state X")
plt.title("Bifurcation curve of self-activating gene (S-shaped bistability)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.2.1_s2.png")
