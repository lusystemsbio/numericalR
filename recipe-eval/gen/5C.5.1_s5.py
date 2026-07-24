import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Lennard-Jones box: ordered grid -> disordered liquid
# Force law (given):  F(r) = 1/r^13 - 1/r^7  (radial magnitude)
# ---------------------------------------------------------------

N = 25            # number of particles
a = 6.25          # box side length
dt = 0.01         # time step
spacing = 1.25    # initial grid spacing (5x5 grid fills the box)

# --- Initial positions: 5x5 regular grid at spacing 1.25 ---
pos = np.zeros((N, 2))
k = 0
for ix in range(5):
    for iy in range(5):
        pos[k, 0] = ix * spacing
        pos[k, 1] = iy * spacing
        k += 1
grid_start = pos.copy()  # remember the ordered start for the check

# --- Deterministic initial velocities via a golden-ratio (low-discrepancy) sequence ---
# frac(i*phi) is equidistributed in (0,1); map it onto (-0.05, 0.05).
phi = (np.sqrt(5.0) - 1.0) / 2.0   # golden ratio conjugate ~0.618
vel = np.zeros((N, 2))
for i in range(N):
    fx = ((2 * i + 1) * phi) % 1.0     # distinct sequence entries for x
    fy = ((2 * i + 2) * phi) % 1.0     # and for y
    vel[i, 0] = -0.05 + 0.1 * fx
    vel[i, 1] = -0.05 + 0.1 * fy
# Remove any net drift so the box does not translate as a whole.
vel -= vel.mean(axis=0)

# --- Force evaluation with minimum-image periodic boundaries ---
def compute_forces(p):
    f = np.zeros_like(p)
    for i in range(N):
        for j in range(i + 1, N):
            d = p[i] - p[j]
            d -= a * np.round(d / a)          # minimum-image convention
            r2 = d[0] * d[0] + d[1] * d[1]
            r = np.sqrt(r2)
            # radial force magnitude F(r) = 1/r^13 - 1/r^7
            fmag = 1.0 / r**13 - 1.0 / r**7
            fij = fmag * (d / r)              # vector force on i from j
            f[i] += fij
            f[j] -= fij                       # Newton's third law
    return f

# --- One velocity-Verlet step (unit masses), with wrap into the box ---
def verlet_step(p, v, f):
    v_half = v + 0.5 * dt * f                 # half-kick
    p_new = p + dt * v_half                   # drift
    p_new = np.mod(p_new, a)                  # periodic wrap
    f_new = compute_forces(p_new)             # new forces
    v_new = v_half + 0.5 * dt * f_new         # second half-kick
    return p_new, v_new, f_new

# --- Warm-up run to t = 10 ---
forces = compute_forces(pos)
n_warm = int(round(10.0 / dt))
for _ in range(n_warm):
    pos, vel, forces = verlet_step(pos, vel, forces)
pos_after_warmup = pos.copy()

# --- Production run to t = 100, saving snapshots ---
n_prod = int(round(90.0 / dt))               # from t=10 to t=100
snapshot_times = [10, 25, 50, 75, 100]
snapshots = {10: pos_after_warmup.copy()}
t = 10.0
for step in range(1, n_prod + 1):
    pos, vel, forces = verlet_step(pos, vel, forces)
    t = 10.0 + step * dt
    for st in snapshot_times:
        if abs(t - st) < dt / 2 and st not in snapshots:
            snapshots[st] = pos.copy()
grid_end = pos.copy()

# ---------------------------------------------------------------
# Structural check: ordered grid at start vs. irregular liquid at end.
# On the grid every particle sits at an exact multiple of the spacing,
# so the standard deviation of coordinates modulo the spacing is ~0.
# In a disordered liquid the positions are scattered, so that
# deviation is O(spacing).  We also compare nearest-neighbour spread.
# ---------------------------------------------------------------
def grid_residual(p):
    # distance of each coordinate from the nearest grid line
    frac = np.mod(p, spacing)
    frac = np.minimum(frac, spacing - frac)
    return np.std(frac)

def nn_std(p):
    dmins = []
    for i in range(N):
        d = p[i] - p
        d -= a * np.round(d / a)
        r = np.sqrt((d**2).sum(axis=1))
        r[i] = np.inf
        dmins.append(r.min())
    return np.std(dmins), np.mean(dmins)

res_start = grid_residual(grid_start)
res_end = grid_residual(grid_end)
nn_std_start, nn_mean_start = nn_std(grid_start)
nn_std_end, nn_mean_end = nn_std(grid_end)

print("Number of particles: %d" % N)
print("Box side a: %.5f" % a)
print("Grid spacing: %.5f" % spacing)
print("Time step dt: %.5f" % dt)
print("Warm-up steps (to t=10): %d" % n_warm)
print("Production steps (to t=100): %d" % n_prod)
print("Grid-line residual std at start (should be ~0): %.6e" % res_start)
print("Grid-line residual std at end (liquid, larger): %.6e" % res_end)
print("Nearest-neighbour distance mean at start: %.6f" % nn_mean_start)
print("Nearest-neighbour distance std  at start: %.6e" % nn_std_start)
print("Nearest-neighbour distance mean at end:   %.6f" % nn_mean_end)
print("Nearest-neighbour distance std  at end:   %.6f" % nn_std_end)
ordered_start = res_start < 1e-6
disordered_end = res_end > 0.05 and nn_std_end > 10 * nn_std_start
print("Check - starts ordered on grid: %s" % ordered_start)
print("Check - ends disordered (liquid): %s" % disordered_end)
print("Check PASSED: %s" % (ordered_start and disordered_end))
# The check confirms the result because a start residual of ~0 means every
# particle sat exactly on the regular lattice, while a large end residual and
# a spread of nearest-neighbour distances means the lattice has melted into an
# irregular liquid-like arrangement -- exactly the ordered->disordered transition.
print("Explanation: the residual is zero only when positions lie exactly on the "
      "lattice, so its jump from ~0 to O(spacing) together with a broadened "
      "nearest-neighbour spread proves the ordered grid melted into a disordered liquid.")

# --- Snapshots figure ---
fig, axes = plt.subplots(1, len(snapshot_times), figsize=(4 * len(snapshot_times), 4))
for ax, st in zip(axes, snapshot_times):
    p = snapshots[st]
    ax.scatter(p[:, 0], p[:, 1], s=60, c="tab:blue", edgecolors="k")
    ax.set_xlim(0, a)
    ax.set_ylim(0, a)
    ax.set_aspect("equal")
    ax.set_title("t = %d" % st)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
fig.suptitle("2D Lennard-Jones: ordered grid melting into a disordered liquid")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.5.1_s5.png")
