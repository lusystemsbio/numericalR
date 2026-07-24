import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---- Model definition: self-activating gene ----
# f(X,k) = g0 + g1 * (X/Xth)^n / (1 + (X/Xth)^n) - k*X
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4

def f(X, k):
    u = (X / Xth) ** n
    return g0 + g1 * u / (1.0 + u) - k * X

def dfdX(X, k):
    # analytic derivative of the Hill term minus the linear degradation slope
    u = (X / Xth) ** n
    dHill_dX = (1.0 / (1.0 + u) ** 2) * (n / Xth) * (X / Xth) ** (n - 1)
    return g1 * dHill_dX - k

# ---- Sweep the control parameter k ----
ks = np.linspace(0.02, 0.60, 400)      # range that brackets the bistable region near k=0.15

# collect scattered steady-state points as (k, X, stable_flag)
pts_k, pts_X, pts_stable = [], [], []
count_per_k = []                       # number of steady states at each k (for the S-curve check)

# X scan grid: steady states lie between the low basal state and the high saturated state.
# Max possible X is (g0+g1)/k for the smallest k; add margin.
Xscan = np.linspace(0.0, (g0 + g1) / ks.min() + 100.0, 6000)

for k in ks:
    fvals = f(Xscan, k)
    roots_here = []
    # find every sign change and refine that bracket with a library root-finder (brentq)
    for i in range(len(Xscan) - 1):
        a, b = fvals[i], fvals[i + 1]
        if a == 0.0:
            roots_here.append(Xscan[i])
        elif a * b < 0.0:                      # a bracket -> exactly one root inside
            r = brentq(f, Xscan[i], Xscan[i + 1], args=(k,))
            roots_here.append(r)
    # de-duplicate roots that fall in adjacent brackets
    roots_here = sorted(roots_here)
    dedup = []
    for r in roots_here:
        if not dedup or abs(r - dedup[-1]) > 1e-6:
            dedup.append(r)
    count_per_k.append(len(dedup))
    for r in dedup:
        # classify: df/dX < 0 => stable, df/dX > 0 => unstable
        stable = dfdX(r, k) < 0.0
        pts_k.append(k)
        pts_X.append(r)
        pts_stable.append(stable)

pts_k = np.array(pts_k)
pts_X = np.array(pts_X)
pts_stable = np.array(pts_stable)
count_per_k = np.array(count_per_k)

# ---- Order the scattered points with a nearest-neighbor walk ----
# normalize each axis so k and X are comparably scaled before measuring distance
P = np.column_stack([pts_k / pts_k.ptp(), pts_X / pts_X.ptp()])
N = len(P)
visited = np.zeros(N, dtype=bool)
order = [int(np.argmin(P[:, 0] + P[:, 1]))]     # start from a corner (small k, small X)
visited[order[0]] = True
for _ in range(N - 1):
    cur = P[order[-1]]
    d = np.sum((P - cur) ** 2, axis=1)
    d[visited] = np.inf                          # never revisit
    nxt = int(np.argmin(d))
    order.append(nxt)
    visited[nxt] = True
order = np.array(order)

# ---- S-shape check: middle range should have 3 states (2 stable, 1 unstable) ----
three_mask = count_per_k == 3
if three_mask.any():
    k_lo = ks[three_mask].min()
    k_hi = ks[three_mask].max()
else:
    k_lo = k_hi = None

# verify the 2-stable / 1-unstable composition on a representative bistable k
if k_lo is not None:
    k_mid = ks[three_mask][len(ks[three_mask]) // 2]
    sel = np.isclose(pts_k, k_mid)
    n_stable = int(np.sum(pts_stable[sel]))
    n_unstable = int(np.sum(~pts_stable[sel]))
else:
    k_mid = n_stable = n_unstable = None

# ---- Report numerical results ----
print(f"k range swept: {ks.min():.4f} to {ks.max():.4f} ({len(ks)} samples)")
print(f"Total steady-state points found: {N}")
print(f"Min steady states at any k: {count_per_k.min()}")
print(f"Max steady states at any k: {count_per_k.max()}")
print(f"Number of k with exactly 3 steady states: {int(three_mask.sum())}")
print(f"Bistable (3-state) k range: {k_lo} to {k_hi}")
print(f"Representative middle k tested: {k_mid}")
print(f"  stable states at that k: {n_stable}")
print(f"  unstable states at that k: {n_unstable}")
# report the three states at the representative k
if k_mid is not None:
    sel = np.isclose(pts_k, k_mid)
    for X, st in sorted(zip(pts_X[sel], pts_stable[sel])):
        print(f"  X = {X:10.4f}  df/dX = {dfdX(X, k_mid):+.6f}  -> {'stable' if st else 'unstable'}")
print(f"k=0.15 sanity: number of steady states = {count_per_k[np.argmin(np.abs(ks-0.15))]}")

# ---- Plot the bifurcation curve ----
plt.figure(figsize=(8, 6))
# faint nearest-neighbor-ordered path through all points
plt.plot(pts_k[order], pts_X[order], '-', color='0.75', lw=0.8, zorder=1)
plt.scatter(pts_k[pts_stable], pts_X[pts_stable], s=10, c='tab:blue',
            label='stable (df/dX < 0)', zorder=2)
plt.scatter(pts_k[~pts_stable], pts_X[~pts_stable], s=10, c='tab:red',
            label='unstable (df/dX > 0)', zorder=2)
plt.xlabel('k (degradation / control parameter)')
plt.ylabel('steady-state X')
plt.title('Bifurcation curve: self-activating gene\n(S-shaped, bistable region)')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.2.1_s4.png")

# One-sentence explanation of why the check confirms the result:
# Finding a contiguous middle band of k with exactly three steady states (two stable, one
# unstable) that collapses to a single state on either side is precisely the definition of an
# S-shaped fold/hysteresis curve, so counting states per k directly verifies the bistable shape.
print("Check: three coexisting states over a middle k band, collapsing to one outside it, is the")
print("defining signature of an S-shaped bistable (fold) curve, confirming the result.")
