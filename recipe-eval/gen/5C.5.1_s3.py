import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Lennard-Jones box: ordered grid -> disordered liquid
# Force magnitude model: F(r) = 1/r^13 - 1/r^7  (given in problem)
# Force vector on i from j = F(r)*(rvec/r) = (1/r^14 - 1/r^8)*rvec
# Minimum-image periodic boundaries, velocity-Verlet integrator.
# ---------------------------------------------------------------

# ----- parameters -----
N = 25              # number of particles
a = 6.25            # square box side
spacing = 1.25      # grid spacing (5x5 grid -> 5*1.25 = 6.25 = a)
dt = 0.01           # time step
t_warm = 10.0       # warm-up end time
t_prod = 100.0      # production end time
mass = 1.0          # unit mass

# ----- initial positions on a regular 5x5 grid -----
# place particles at cell centers so they tile the box periodically
pos0 = np.zeros((N, 2))
k = 0
for ix in range(5):
    for iy in range(5):
        pos0[k, 0] = (ix + 0.5) * spacing
        pos0[k, 1] = (iy + 0.5) * spacing
        k += 1
pos = pos0.copy()

# ----- deterministic initial velocities via a golden-ratio sequence -----
# low-discrepancy sequence frac(n*phi) in [0,1) -> mapped to (-0.05, 0.05)
phi = (1.0 + np.sqrt(5.0)) / 2.0
seq = np.mod(np.arange(1, 2 * N + 1) * phi, 1.0)   # 2N values, one per velocity component
vel = np.zeros((N, 2))
vel[:, 0] = -0.05 + 0.1 * seq[0:N]        # vx
vel[:, 1] = -0.05 + 0.1 * seq[N:2 * N]    # vy
# remove any net drift so the whole box does not translate
vel -= vel.mean(axis=0)

# ----- force routine with minimum-image convention -----
def compute_forces(p):
    f = np.zeros_like(p)
    for i in range(N):
        for j in range(i + 1, N):
            d = p[i] - p[j]
            # minimum image: wrap separation into (-a/2, a/2]
            d -= a * np.round(d / a)
            r2 = d[0] * d[0] + d[1] * d[1]
            r = np.sqrt(r2)
            # scalar factor so that force vector = fac * d  (= F(r)*d/r)
            fac = 1.0 / r**14 - 1.0 / r**8
            f[i] += fac * d
            f[j] -= fac * d
    return f

# ----- velocity-Verlet step with periodic wrap (from 5C.4) -----
def vv_step(p, v, f):
    # half-kick position update
    p = p + v * dt + 0.5 * (f / mass) * dt * dt
    p = np.mod(p, a)                      # periodic boundaries: keep in [0, a)
    f_new = compute_forces(p)
    v = v + 0.5 * (f + f_new) / mass * dt # velocity update using avg force
    return p, v, f_new

# ----- structural check helper: nearest-neighbour distances -----
def nn_distances(p):
    nn = np.empty(N)
    for i in range(N):
        best = np.inf
        for j in range(N):
            if i == j:
                continue
            d = p[i] - p[j]
            d -= a * np.round(d / a)
            rr = np.hypot(d[0], d[1])
            if rr < best:
                best = rr
        nn[i] = best
    return nn

# record the initial grid
pos_initial = pos.copy()
nn_init = nn_distances(pos_initial)

# ----- warm-up run to t = 10 -----
forces = compute_forces(pos)
n_warm = int(round(t_warm / dt))
for _ in range(n_warm):
    pos, vel, forces = vv_step(pos, vel, forces)
pos_warm = pos.copy()

# ----- production run to t = 100, capturing snapshots -----
snap_times = [25.0, 50.0, 75.0, 100.0]
snap_steps = {int(round(t / dt)): t for t in snap_times}
snapshots = {}
n_prod = int(round((t_prod - t_warm) / dt))
step_global = n_warm
for _ in range(n_prod):
    pos, vel, forces = vv_step(pos, vel, forces)
    step_global += 1
    if step_global in snap_steps:
        snapshots[snap_times[len(snapshots)]] = pos.copy()
pos_final = pos.copy()
nn_final = nn_distances(pos_final)

# ----- structural averages for the check -----
mean_nn_init = nn_init.mean()
std_nn_init = nn_init.std()
mean_nn_final = nn_final.mean()
std_nn_final = nn_final.std()

print(f"Number of particles: {N}")
print(f"Box side a: {a}")
print(f"Grid spacing: {spacing}")
print(f"Warm-up time: {t_warm}, Production time: {t_prod}, dt: {dt}")
print(f"Initial mean nearest-neighbour distance: {mean_nn_init:.6f}")
print(f"Initial std of nearest-neighbour distance: {std_nn_init:.6e}")
print(f"Final   mean nearest-neighbour distance: {mean_nn_final:.6f}")
print(f"Final   std of nearest-neighbour distance: {std_nn_final:.6f}")
print(f"Ordered-start check (std ~ 0): {std_nn_init < 1e-9}")
print(f"Disordered-end check (std grew): {std_nn_final > 10 * max(std_nn_init, 1e-12)}")

# ----- snapshots figure -----
panels = [("t = 0 (grid)", pos_initial),
          ("t = 10 (warm-up)", pos_warm),
          ("t = 25", snapshots[25.0]),
          ("t = 50", snapshots[50.0]),
          ("t = 75", snapshots[75.0]),
          ("t = 100 (final)", snapshots[100.0])]
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
for ax, (title, p) in zip(axes.ravel(), panels):
    ax.scatter(p[:, 0], p[:, 1], s=60, c="tab:blue", edgecolors="k")
    ax.set_xlim(0, a)
    ax.set_ylim(0, a)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
fig.suptitle("2D Lennard-Jones: ordered grid melting into a disordered liquid")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.5.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: the check confirms the result because a near-zero spread in "
      "nearest-neighbour distances certifies the identical, regular grid at the "
      "start while a much larger spread at the end certifies an irregular liquid-"
      "like arrangement, and these structural averages are reproducible across "
      "implementations even though chaotic dynamics make the exact final "
      "coordinates differ between R and Python.")
