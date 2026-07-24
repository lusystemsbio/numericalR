import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Vector field f(state) = (dX/dt, dY/dt)
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# ----------------------------------------------------------------------
# 1. Locate the saddle: a fixed point f(state)=0 near (190, 127)
# ----------------------------------------------------------------------
saddle = fsolve(f, [190.0, 127.0], full_output=False)
print(f"Saddle fixed point X: {saddle[0]:.6f}")
print(f"Saddle fixed point Y: {saddle[1]:.6f}")
print(f"Residual dX/dt at saddle: {f(saddle)[0]:.3e}")
print(f"Residual dY/dt at saddle: {f(saddle)[1]:.3e}")

# ----------------------------------------------------------------------
# 2. Jacobian at the saddle (finite differences) and its eigenstructure.
#    The separatrix is the STABLE manifold of the saddle (eigenvalue<0).
# ----------------------------------------------------------------------
def jacobian(state, h=1e-4):
    J = np.zeros((2, 2))
    for j in range(2):
        d = np.zeros(2); d[j] = h
        J[:, j] = (f(state + d) - f(state - d)) / (2 * h)
    return J

J = jacobian(saddle)
evals, evecs = np.linalg.eig(J)
print(f"Jacobian eigenvalue 1: {evals[0]:.6f}")
print(f"Jacobian eigenvalue 2: {evals[1]:.6f}")

# stable eigenvector = eigenvector of the negative eigenvalue
stable_idx = int(np.argmin(evals.real))
v_stable = evecs[:, stable_idx].real
v_stable = v_stable / np.linalg.norm(v_stable)
print(f"Stable eigenvalue (separatrix direction): {evals[stable_idx].real:.6f}")
print(f"Stable eigenvector: [{v_stable[0]:.6f}, {v_stable[1]:.6f}]")

# ----------------------------------------------------------------------
# 3. Time-REVERSED flow: dX/dt = -f(X).
#    Under reversal the stable manifold becomes repelling, so points seeded
#    just beside the saddle along +/- v_stable are carried OUTWARD, tracing
#    the separatrix. Seed ten points (five each side, small offsets).
# ----------------------------------------------------------------------
def f_reversed(t, state):
    return -f(state)

offsets = np.array([0.1, 0.5, 1.0, 2.0, 3.0])   # distances from saddle
seeds = []
for d in offsets:
    seeds.append(saddle + d * v_stable)          # one branch
    seeds.append(saddle - d * v_stable)          # opposite branch
print(f"Number of seed points near saddle: {len(seeds)}")

t_span = (0.0, 4000.0)
t_eval = np.linspace(*t_span, 6000)
sep_branches = []
for s in seeds:
    sol = solve_ivp(f_reversed, t_span, s, t_eval=t_eval,
                    rtol=1e-9, atol=1e-9, dense_output=False)
    sep_branches.append(sol.y)
print(f"Integrated {len(sep_branches)} reversed trajectories to trace the separatrix.")

# ----------------------------------------------------------------------
# 4. Forward flow for sample trajectories and for the basin check.
# ----------------------------------------------------------------------
def f_forward(t, state):
    return f(state)

def integrate_forward(x0, T=3000.0, n=3000):
    sol = solve_ivp(f_forward, (0, T), x0, t_eval=np.linspace(0, T, n),
                    rtol=1e-8, atol=1e-8)
    return sol.y

# Two stable nodes (attractors), found from far-apart starts
attr_lowX = integrate_forward([10.0, 400.0])[:, -1]   # low X, high Y
attr_hiX  = integrate_forward([600.0, 10.0])[:, -1]   # high X, low Y
print(f"Attractor A (low X):  X={attr_lowX[0]:.3f}, Y={attr_lowX[1]:.3f}")
print(f"Attractor B (high X): X={attr_hiX[0]:.3f}, Y={attr_hiX[1]:.3f}")

# A handful of illustrative forward trajectories
sample_starts = [[50, 350], [30, 500], [500, 300], [550, 100],
                 [200, 300], [250, 80], [100, 200], [400, 250]]
sample_trajs = [integrate_forward(s0) for s0 in sample_starts]

# ----------------------------------------------------------------------
# 5. Basin check: many random starts. Classify each by which attractor it
#    reaches; none should cross the separatrix (the boundary between basins).
# ----------------------------------------------------------------------
rng = np.random.default_rng(0)
N = 300
starts = rng.uniform([0, 0], [700, 550], size=(N, 2))
ends = np.array([integrate_forward(s0, T=4000, n=2)[:, -1] for s0 in starts])
# distance to each attractor decides the basin label
d_low = np.linalg.norm(ends - attr_lowX, axis=1)
d_hi  = np.linalg.norm(ends - attr_hiX, axis=1)
basin = np.where(d_low < d_hi, 0, 1)   # 0 -> attractor A, 1 -> attractor B
print(f"Random starts landing in basin A (low X): {np.sum(basin == 0)}")
print(f"Random starts landing in basin B (high X): {np.sum(basin == 1)}")

# ----------------------------------------------------------------------
# 6. Plots
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# -- Left: separatrix + sample forward trajectories + nullclines --
Xg = np.linspace(0, 700, 400)
Yg = np.linspace(0, 550, 400)
XX, YY = np.meshgrid(Xg, Yg)
dXX = gX0 + gX1 / (1 + (YY / Yth) ** nY) - kX * XX
dYY = gY0 + gY1 / (1 + (XX / Xth) ** nX) - kY * YY
ax1.contour(XX, YY, dXX, levels=[0], colors='tab:blue', linewidths=1, alpha=0.6)
ax1.contour(XX, YY, dYY, levels=[0], colors='tab:orange', linewidths=1, alpha=0.6)

for b in sep_branches:
    ax1.plot(b[0], b[1], 'k-', lw=2.2, zorder=5)
ax1.plot([], [], 'k-', lw=2.2, label='separatrix (reversed flow)')
for tr in sample_trajs:
    ax1.plot(tr[0], tr[1], '-', color='gray', lw=1.0, alpha=0.8)
    ax1.plot(tr[0, 0], tr[1, 0], '.', color='green', ms=6)
ax1.plot([], [], '-', color='gray', lw=1.0, label='forward trajectories')
ax1.plot(*saddle, 'rs', ms=10, label='saddle')
ax1.plot(*attr_lowX, 'k*', ms=15, label='attractors')
ax1.plot(*attr_hiX, 'k*', ms=15)
ax1.set_xlim(0, 700); ax1.set_ylim(0, 550)
ax1.set_xlabel('X'); ax1.set_ylabel('Y')
ax1.set_title('Toggle switch: separatrix and sample trajectories')
ax1.legend(loc='upper right', fontsize=8)

# -- Right: basin check --
ax2.scatter(starts[basin == 0, 0], starts[basin == 0, 1],
            s=12, color='tab:cyan', label='ends at attractor A')
ax2.scatter(starts[basin == 1, 0], starts[basin == 1, 1],
            s=12, color='tab:pink', label='ends at attractor B')
for b in sep_branches:
    ax2.plot(b[0], b[1], 'k-', lw=2.5, zorder=5)
ax2.plot([], [], 'k-', lw=2.5, label='separatrix')
ax2.plot(*saddle, 'rs', ms=10)
ax2.set_xlim(0, 700); ax2.set_ylim(0, 550)
ax2.set_xlabel('X'); ax2.set_ylabel('Y')
ax2.set_title('Basin check: separatrix divides the two basins')
ax2.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3F.2.1_s1.png", dpi=130)

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# Because solution trajectories of an autonomous ODE are unique and cannot
# cross, the separatrix (the saddle's stable manifold) acts as an
# impassable wall: every random start on one side flows to one attractor
# and every start on the other side to the other, so the traced curve is
# exactly the boundary between the two basins of attraction.
print("Explanation: trajectory uniqueness forbids crossings, so the curve "
      "that cleanly separates starts reaching attractor A from those "
      "reaching attractor B is precisely the basin boundary.")
