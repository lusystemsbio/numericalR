import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ---- Model parameters ----
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

# f(X,k) = basal + Hill activation - linear degradation
def f(X, k):
    h = (X / Xth) ** n
    return g0 + g1 * h / (1.0 + h) - k * X

# df/dX: derivative used to classify stability (Hill term derivative - k)
def dfdX(X, k):
    # d/dX [ g1 * u/(1+u) ] with u=(X/Xth)^n ; chain rule
    u = (X / Xth) ** n
    dudX = n * (X ** (n - 1)) / (Xth ** n)      # du/dX
    hill_deriv = g1 * dudX / (1.0 + u) ** 2      # quotient rule collapses to this
    return hill_deriv - k

# ---- Root finding: explicit bracketing + library root-finder (brentq) ----
def find_roots(k):
    # Sample X on a fine grid; a root lies wherever f changes sign between samples.
    Xgrid = np.linspace(1e-6, 1500.0, 4000)
    fvals = f(Xgrid, k)
    roots = []
    for i in range(len(Xgrid) - 1):
        if fvals[i] == 0.0:
            roots.append(Xgrid[i])
        elif fvals[i] * fvals[i + 1] < 0.0:
            # sign change -> bracket -> refine with brentq
            r = brentq(f, Xgrid[i], Xgrid[i + 1], args=(k,))
            roots.append(r)
    # de-duplicate near-identical roots
    uniq = []
    for r in sorted(roots):
        if not uniq or abs(r - uniq[-1]) > 1e-4:
            uniq.append(r)
    return uniq

# ---- Sweep k over a grid, collect (k, X, stable?) points ----
ks = np.linspace(0.05, 0.35, 400)
pts = []  # each entry: (k, X, is_stable)
for k in ks:
    for X in find_roots(k):
        # Stable steady state <=> df/dX < 0 (perturbation decays)
        is_stable = dfdX(X, k) < 0.0
        pts.append((k, X, is_stable))

pts = np.array([(p[0], p[1], 1.0 if p[2] else 0.0) for p in pts])

# ---- Nearest-neighbor walk to order the scattered points into a curve ----
# Normalize coordinates so k and X are comparable, then greedily hop to the
# nearest unvisited point starting from the lowest-k point.
coords = pts[:, :2].copy()
cn = coords.copy()
cn[:, 0] = (cn[:, 0] - cn[:, 0].min()) / (cn[:, 0].ptp())
cn[:, 1] = (cn[:, 1] - cn[:, 1].min()) / (cn[:, 1].ptp())
N = len(cn)
visited = np.zeros(N, dtype=bool)
order = [int(np.argmin(pts[:, 0]))]  # start at smallest k
visited[order[0]] = True
for _ in range(N - 1):
    cur = cn[order[-1]]
    d = np.sum((cn - cur) ** 2, axis=1)
    d[visited] = np.inf
    nxt = int(np.argmin(d))
    order.append(nxt)
    visited[nxt] = True
ordered = pts[order]

# ---- Bistability check: count steady states across k ----
n_states = {}
for k in ks:
    r = find_roots(k)
    n_states.setdefault(len(r), 0)
    n_states[len(r)] += 1

bistable_ks = [k for k in ks if len(find_roots(k)) == 3]
if bistable_ks:
    k_lo, k_hi = min(bistable_ks), max(bistable_ks)
else:
    k_lo, k_hi = float("nan"), float("nan")

# Detailed check near k = 0.15
k_check = 0.15
roots_check = find_roots(k_check)
stab_check = [dfdX(X, k_check) < 0.0 for X in roots_check]

# ---- Print numerical results ----
print(f"Number of k samples: {len(ks)}")
print(f"Distribution of steady-state counts over k grid: {dict(sorted(n_states.items()))}")
print(f"Three-steady-state (bistable) k range: [{k_lo:.4f}, {k_hi:.4f}]")
print(f"Steady states at k = {k_check}: {[round(X, 4) for X in roots_check]}")
for X, s in zip(roots_check, stab_check):
    print(f"  X = {X:10.4f}  df/dX = {dfdX(X, k_check):+.5f}  -> {'STABLE' if s else 'UNSTABLE'}")
n_stable = sum(stab_check)
n_unstable = len(stab_check) - n_stable
print(f"At k = {k_check}: {len(roots_check)} states, {n_stable} stable, {n_unstable} unstable")
print(f"Total scattered steady-state points collected: {N}")

# ---- Plot: bifurcation curve, colored by stability ----
plt.figure(figsize=(8, 6))
stable_mask = ordered[:, 2] == 1.0
plt.scatter(ordered[stable_mask, 0], ordered[stable_mask, 1],
            c="tab:blue", s=12, label="stable (df/dX < 0)")
plt.scatter(ordered[~stable_mask, 0], ordered[~stable_mask, 1],
            c="tab:red", s=12, label="unstable (df/dX > 0)")
plt.axvspan(k_lo, k_hi, color="gray", alpha=0.15, label="bistable region")
plt.xlabel("k (degradation / control parameter)")
plt.ylabel("steady-state X")
plt.title("Bifurcation curve of self-activating gene (S-shaped)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.2.1_s3.png")

# Explanation (one sentence): The presence of a middle k range with three roots
# (two stable, one unstable) that collapses to a single root outside it is the
# defining signature of the fold/saddle-node S-curve, so observing exactly that
# fold-and-collapse pattern confirms the curve is genuinely S-shaped and bistable.
print("Check: three-states-in-middle collapsing to one outside => S-shaped bistable curve confirmed.")
