import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Vector field f(state) = (dX/dt, dY/dt)
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# Wrappers for solve_ivp (signature t, y)
def forward_rhs(t, s):          # normal (forward) flow
    return f(s)

def reversed_rhs(t, s):         # time-reversed flow: dX/dt = -f(X)
    return -f(s)

# ------------------------------------------------------------------
# 1. Locate the saddle point near (190, 127) and get its Jacobian
# ------------------------------------------------------------------
saddle = fsolve(f, [190.0, 127.0], full_output=False)
print(f"Saddle point (X*, Y*): ({saddle[0]:.6f}, {saddle[1]:.6f})")

# Numerical Jacobian of f at the saddle (finite differences)
def jacobian(func, s, eps=1e-6):
    J = np.zeros((2, 2))
    for j in range(2):
        ds = np.zeros(2); ds[j] = eps
        J[:, j] = (func(s + ds) - func(s - ds)) / (2 * eps)
    return J

J = jacobian(f, saddle)
evals, evecs = np.linalg.eig(J)
print(f"Jacobian eigenvalues at saddle: {evals[0]:.6f}, {evals[1]:.6f}")

# A saddle has one positive and one negative eigenvalue.
# The separatrix is the STABLE manifold (eigenvalue < 0). Under the
# time-REVERSED flow that stable direction becomes unstable, so points
# seeded just beside the saddle are carried OUTWARD along the separatrix.
stable_idx = np.argmin(evals.real)          # most negative eigenvalue
stable_dir = evecs[:, stable_idx].real
stable_dir = stable_dir / np.linalg.norm(stable_dir)
print(f"Stable eigenvector (separatrix direction): "
      f"({stable_dir[0]:.6f}, {stable_dir[1]:.6f})")

# ------------------------------------------------------------------
# 2. Seed ten points just beside the saddle and integrate the
#    time-reversed system to trace the separatrix.
# ------------------------------------------------------------------
# Five offsets along +stable_dir and five along -stable_dir cover both
# branches of the separatrix that emanate from the saddle.
offsets = np.array([0.2, 0.5, 1.0, 2.0, 4.0])
seeds = []
for d in offsets:
    seeds.append(saddle + d * stable_dir)   # one branch
    seeds.append(saddle - d * stable_dir)   # other branch
seeds = np.array(seeds)
print(f"Number of seed points near saddle: {len(seeds)}")

t_span = (0.0, 4000.0)
t_eval = np.linspace(*t_span, 4000)
sep_branches = []
for k, s0 in enumerate(seeds):
    sol = solve_ivp(reversed_rhs, t_span, s0, t_eval=t_eval,
                    rtol=1e-9, atol=1e-9, max_step=1.0)
    sep_branches.append(sol.y)
    end = sol.y[:, -1]
    print(f"Reversed seed {k:2d} start=({s0[0]:8.3f},{s0[1]:8.3f}) "
          f"end=({end[0]:8.3f},{end[1]:8.3f})")

# ------------------------------------------------------------------
# 3. Find the two stable fixed points (attractors) for reference
# ------------------------------------------------------------------
attr_lo = fsolve(f, [80.0, 300.0])    # high-Y / low-X state
attr_hi = fsolve(f, [500.0, 40.0])    # high-X / low-Y state
print(f"Attractor A (X,Y): ({attr_lo[0]:.4f}, {attr_lo[1]:.4f})")
print(f"Attractor B (X,Y): ({attr_hi[0]:.4f}, {attr_hi[1]:.4f})")

# ------------------------------------------------------------------
# 4. Forward trajectories from random starts (basin check).
#    Each must settle into ONE attractor and never cross the separatrix.
# ------------------------------------------------------------------
rng = np.random.default_rng(0)
Xmax, Ymax = 600.0, 450.0
n_random = 40
fwd_trajs = []
labels = []
for i in range(n_random):
    s0 = np.array([rng.uniform(0, Xmax), rng.uniform(0, Ymax)])
    sol = solve_ivp(forward_rhs, (0.0, 3000.0), s0,
                    t_eval=np.linspace(0, 3000, 1500),
                    rtol=1e-8, atol=1e-8)
    end = sol.y[:, -1]
    # classify by which attractor it converged to
    dA = np.linalg.norm(end - attr_lo)
    dB = np.linalg.norm(end - attr_hi)
    labels.append(0 if dA < dB else 1)
    fwd_trajs.append(sol.y)
n_A = labels.count(0); n_B = labels.count(1)
print(f"Random forward starts -> Attractor A: {n_A}, Attractor B: {n_B}")

# ------------------------------------------------------------------
# 5. Explicit crossing check: does any forward trajectory cross the
#    separatrix curve? Build the separatrix as an ordered polyline and
#    test the sign of each forward point relative to it.
# ------------------------------------------------------------------
# Assemble the full separatrix polyline (both branches, ordered).
branch_pos = np.hstack([sep_branches[2*i]   for i in range(len(offsets))])
branch_neg = np.hstack([sep_branches[2*i+1] for i in range(len(offsets))])
# order each branch by distance from the saddle so it forms a line
def order_branch(b):
    d = np.linalg.norm(b - saddle[:, None], axis=0)
    return b[:, np.argsort(d)]
sep_line = np.hstack([order_branch(branch_neg)[:, ::-1],
                      order_branch(branch_pos)])

# Signed side of point p relative to the separatrix: use the nearest
# segment and the cross product to decide left/right.
def side_of(p, line):
    diffs = line[:, 1:] - line[:, :-1]
    mids = 0.5 * (line[:, 1:] + line[:, :-1])
    d2 = np.sum((mids - p[:, None]) ** 2, axis=0)
    j = np.argmin(d2)
    seg = diffs[:, j]
    rel = p - mids[:, j]
    return np.sign(seg[0] * rel[1] - seg[1] * rel[0])

# For each forward trajectory, the sign of its side should stay constant
# (never flips) -> it never crosses the separatrix.
crossings = 0
for y in fwd_trajs:
    signs = np.array([side_of(y[:, m], sep_line) for m in range(0, y.shape[1], 20)])
    signs = signs[signs != 0]
    if signs.size and np.any(signs[1:] != signs[:-1]):
        crossings += 1
print(f"Forward trajectories that cross the separatrix: {crossings} / {n_random}")

# Also confirm the two attractors sit on opposite sides of the curve.
print(f"Attractor A side: {side_of(attr_lo, sep_line):+.0f}")
print(f"Attractor B side: {side_of(attr_hi, sep_line):+.0f}")

# ------------------------------------------------------------------
# 6. Phase-plane plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 7))

# vector-field direction arrows
Xg, Yg = np.meshgrid(np.linspace(0, Xmax, 22), np.linspace(0, Ymax, 22))
U = np.zeros_like(Xg); V = np.zeros_like(Yg)
for i in range(Xg.shape[0]):
    for j in range(Xg.shape[1]):
        d = f([Xg[i, j], Yg[i, j]])
        n = np.hypot(*d) + 1e-12
        U[i, j], V[i, j] = d / n
ax.quiver(Xg, Yg, U, V, color="0.8", pivot="mid", width=0.002)

# forward trajectories, colored by basin
for y, lab in zip(fwd_trajs, labels):
    ax.plot(y[0], y[1], color=("tab:blue" if lab == 0 else "tab:orange"),
            lw=0.8, alpha=0.6)

# separatrix (from reversed integration)
ax.plot(sep_line[0], sep_line[1], "k-", lw=2.5, label="Separatrix (reversed flow)")

# markers
ax.plot(*saddle, "ks", ms=10, label="Saddle")
ax.plot(*attr_lo, "o", color="tab:blue", ms=11, mec="k", label="Attractor A")
ax.plot(*attr_hi, "o", color="tab:orange", ms=11, mec="k", label="Attractor B")
ax.plot(seeds[:, 0], seeds[:, 1], "r.", ms=6, label="Seeds near saddle")

ax.set_xlim(0, Xmax); ax.set_ylim(0, Ymax)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle switch: separatrix and basins of attraction")
ax.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3F.2.1_s3.png", dpi=130)
print("Saved phase-plane plot.")

# ------------------------------------------------------------------
# Why the check confirms the result:
# ------------------------------------------------------------------
print("Explanation: Because no forward trajectory ever crosses the "
      "reversed-flow curve and trajectories on opposite sides converge to "
      "different attractors, that curve is exactly the invariant boundary "
      "separating the two basins of attraction.")
