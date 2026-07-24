import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---- Model definition ----------------------------------------------------
# Self-activating gene: basal transcription + excitatory Hill term - linear degradation
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    # rate of change of X; steady states are the roots f(X,k)=0
    return g0 + g1 * (X / Xth) ** n / (1.0 + (X / Xth) ** n) - k * X

def dfdX(X, k):
    # analytic derivative of f w.r.t. X; its sign classifies stability
    u = (X / Xth) ** n
    dhill = g1 * (n / X) * u / (1.0 + u) ** 2   # d/dX of Hill term
    return dhill - k

# ---- Root finding on a k-grid -------------------------------------------
ks = np.linspace(0.05, 0.30, 400)      # sweep control parameter (bistable near 0.15)

pts_k, pts_X, pts_stable = [], [], []  # scattered steady-state points

for k in ks:
    # bracket the search: X>0, upper bound where degradation dominates supply
    Xmax = (g0 + g1) / k + 10.0
    Xgrid = np.linspace(1e-6, Xmax, 4000)
    fg = f(Xgrid, k)
    # scan for sign changes, refine each bracket with a library root-finder (brentq)
    for i in range(len(Xgrid) - 1):
        if fg[i] == 0.0:
            root = Xgrid[i]
        elif fg[i] * fg[i + 1] < 0.0:
            root = brentq(f, Xgrid[i], Xgrid[i + 1], args=(k,))
        else:
            continue
        stable = dfdX(root, k) < 0.0   # df/dX<0 => stable, >0 => unstable
        pts_k.append(k); pts_X.append(root); pts_stable.append(stable)

pts_k = np.array(pts_k); pts_X = np.array(pts_X); pts_stable = np.array(pts_stable)

# ---- Nearest-neighbor walk to order the scattered points ----------------
# Normalize coordinates so k and X are comparable distances
kn = (pts_k - pts_k.min()) / (pts_k.ptp())
Xn = (pts_X - pts_X.min()) / (pts_X.ptp())
coords = np.column_stack([kn, Xn])

N = len(coords)
visited = np.zeros(N, dtype=bool)
order = [0]                 # start at the first point
visited[0] = True
for _ in range(N - 1):
    cur = coords[order[-1]]
    d = np.sum((coords - cur) ** 2, axis=1)
    d[visited] = np.inf    # ignore already-visited points
    nxt = int(np.argmin(d))
    order.append(nxt); visited[nxt] = True
order = np.array(order)

ordered_k = pts_k[order]; ordered_X = pts_X[order]; ordered_stable = pts_stable[order]

# ---- S-shape / bistability check ----------------------------------------
# Count steady states at each k value; bistable region has 3 (2 stable, 1 unstable)
counts = {}
for k in ks:
    mask = np.isclose(pts_k, k)
    if mask.sum() > 0:
        counts[k] = (mask.sum(), int(pts_stable[mask].sum()), int((~pts_stable[mask]).sum()))

triple_ks = [k for k, c in counts.items() if c[0] == 3]
single_ks = [k for k, c in counts.items() if c[0] == 1]

print("Total steady-state points found:", N)
print("k range swept:", ks.min(), "to", ks.max())
if triple_ks:
    print("Bistable (3 steady states) k-range:", min(triple_ks), "to", max(triple_ks))
    # verify composition 2 stable + 1 unstable in the tristate region
    kmid = triple_ks[len(triple_ks) // 2]
    print("At k =", kmid, "-> (total, stable, unstable) =", counts[kmid])
else:
    print("No 3-steady-state region detected")
print("Number of k with 1 steady state:", len(single_ks))
print("Number of k with 3 steady states:", len(triple_ks))

# report the saddle-node fold points (edges of the bistable window)
if triple_ks:
    print("Lower fold (k) near:", min(triple_ks))
    print("Upper fold (k) near:", max(triple_ks))

# sample the three branches at the bistable midpoint
if triple_ks:
    m = np.isclose(pts_k, kmid)
    Xs = np.sort(pts_X[m])
    print("At k =", kmid, "steady-state X values (sorted):", Xs)

# ---- Plot the bifurcation curve -----------------------------------------
plt.figure(figsize=(8, 6))
st = pts_stable
plt.scatter(pts_k[st], pts_X[st], c="tab:blue", s=12, label="stable (df/dX < 0)")
plt.scatter(pts_k[~st], pts_X[~st], c="tab:red", s=12, label="unstable (df/dX > 0)")
# draw the nearest-neighbor ordered walk connecting the points
plt.plot(ordered_k, ordered_X, color="gray", lw=0.6, alpha=0.6, zorder=0,
         label="nearest-neighbor walk")
plt.xlabel("k (degradation / control parameter)")
plt.ylabel("steady-state X")
plt.title("Bifurcation curve of self-activating gene (S-shaped bistability)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Finding a contiguous middle k-range with exactly three "
      "steady states (two stable, one unstable) flanked by single-state regions "
      "is the defining signature of an S-shaped fold/hysteresis curve, so it "
      "confirms the bistable bifurcation.")
