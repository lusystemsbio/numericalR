import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---- Model definition ------------------------------------------------------
# Self-activating gene: basal transcription + excitatory Hill activation - linear degradation
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    # rate of change of X: production (basal + Hill) minus degradation k*X
    return g0 + g1 * (X / Xth)**n / (1.0 + (X / Xth)**n) - k * X

def dfdX(X, k):
    # analytic derivative df/dX; sign classifies stability of a steady state
    u = (X / Xth)**n
    hill_deriv = g1 * (n / X) * u / (1.0 + u)**2  # d/dX of Hill term
    return hill_deriv - k

# ---- Sweep k on a grid, find every root of f(X,k)=0 at each k ----------------
k_grid = np.linspace(0.05, 0.35, 400)

# X search grid: bracket sign changes of f, then refine each bracket with brentq
X_scan = np.linspace(1e-3, 600.0, 4000)

k_pts, X_pts, stable_flags = [], [], []
for k in k_grid:
    vals = f(X_scan, k)
    # locate sub-intervals where f changes sign -> a root lives inside
    sign_change = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    for i in sign_change:
        try:
            Xr = brentq(f, X_scan[i], X_scan[i + 1], args=(k,))  # library root-finder
        except ValueError:
            continue
        k_pts.append(k)
        X_pts.append(Xr)
        # steady state stable if df/dX < 0 (perturbations decay)
        stable_flags.append(dfdX(Xr, k) < 0)

k_pts = np.array(k_pts)
X_pts = np.array(X_pts)
stable_flags = np.array(stable_flags)

# ---- Order the scattered points with a nearest-neighbor walk -----------------
# normalize coordinates so k and X contribute comparably to distance
kn = (k_pts - k_pts.min()) / (k_pts.ptp())
Xn = (X_pts - X_pts.min()) / (X_pts.ptp())
pts = np.column_stack([kn, Xn])

N = len(pts)
visited = np.zeros(N, dtype=bool)
order = [0]
visited[0] = True
for _ in range(N - 1):
    cur = order[-1]
    d = np.sum((pts - pts[cur])**2, axis=1)
    d[visited] = np.inf
    nxt = int(np.argmin(d))  # closest unvisited point
    order.append(nxt)
    visited[nxt] = True
order = np.array(order)

# ---- Check: count steady states per k to confirm the S-shape -----------------
counts = {}
for k in k_grid:
    m = np.isclose(k_pts, k)
    counts[k] = int(m.sum())
n_states = np.array([counts[k] for k in k_grid])

k_bistable = k_grid[n_states == 3]
print("Number of sampled k values:", len(k_grid))
print("Total steady-state points found:", N)
print("Max steady states at any single k:", int(n_states.max()))
if k_bistable.size > 0:
    print("Bistable k range (3 steady states) lower bound:", float(k_bistable.min()))
    print("Bistable k range (3 steady states) upper bound:", float(k_bistable.max()))
    print("Number of k values with 3 steady states:", int(k_bistable.size))
    print("Number of k values with 1 steady state:", int((n_states == 1).sum()))
else:
    print("No bistable region detected")

# report stability breakdown within the bistable region
mid_k = k_bistable[len(k_bistable) // 2] if k_bistable.size else k_grid[len(k_grid)//2]
m = np.isclose(k_pts, mid_k)
Xs_mid = np.sort(X_pts[m])
print("At mid k =", float(mid_k), "steady states X:", [float(v) for v in Xs_mid])
print("  stable count:", int(stable_flags[m].sum()), "unstable count:", int((~stable_flags[m]).sum()))

# ---- Plot the bifurcation curve ----------------------------------------------
plt.figure(figsize=(8, 6))
# nearest-neighbor ordered path (light guide line through the S-curve)
plt.plot(k_pts[order], X_pts[order], '-', color='0.8', lw=1, zorder=1, label='NN-ordered path')
st = stable_flags
plt.scatter(k_pts[st], X_pts[st], c='tab:blue', s=12, zorder=2, label='stable (df/dX<0)')
plt.scatter(k_pts[~st], X_pts[~st], c='tab:red', s=12, zorder=2, label='unstable (df/dX>0)')
plt.xlabel('k (degradation / control parameter)')
plt.ylabel('steady-state X')
plt.title('Bifurcation curve: self-activating gene')
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: Finding a contiguous middle k-range with exactly three "
      "steady states (two stable flanking one unstable) that collapses to a single "
      "steady state on both sides is the defining signature of an S-shaped "
      "(fold/saddle-node) bistable bifurcation, so observing it confirms the curve is S-shaped.")
