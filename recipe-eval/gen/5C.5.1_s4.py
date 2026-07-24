import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------
N = 25            # number of particles
a = 6.25          # side length of the square, periodic box
n_side = 5        # 5 x 5 grid
spacing = 1.25    # grid spacing (n_side * spacing = a)
dt = 0.01         # integration time step
t_warm = 10.0     # warm-up time
t_prod = 100.0    # production time (total)
mass = 1.0        # unit mass for every particle

# -----------------------------------------------------------------------------
# Lennard-Jones force.  We are told the (scalar, radial) force law is
#   F(r) = 1/r^13 - 1/r^7
# This is the magnitude of the force along the line joining two particles.
# The vector force on particle i due to j is F(r) * (r_ij / r), where
# r_ij is the minimum-image separation vector and r = |r_ij|.
# -----------------------------------------------------------------------------
def lj_force_scalar(r):
    # radial force magnitude (positive = repulsive along r_ij direction)
    return 1.0 / r**13 - 1.0 / r**7

def compute_forces(pos):
    # returns the (N,2) array of total forces using minimum-image convention
    forces = np.zeros_like(pos)
    for i in range(N):
        for j in range(i + 1, N):
            # separation vector, then fold into [-a/2, a/2) via minimum image
            d = pos[i] - pos[j]
            d -= a * np.round(d / a)
            r = np.hypot(d[0], d[1])
            f = lj_force_scalar(r)          # scalar radial force
            fv = f * d / r                  # vector force on i from j
            forces[i] += fv                 # Newton's third law
            forces[j] -= fv
    return forces

# -----------------------------------------------------------------------------
# One velocity-Verlet step with periodic wrapping (the integrator from 5C.4)
# -----------------------------------------------------------------------------
def velocity_verlet_step(pos, vel, forces):
    # half-step: advance positions using current velocity and acceleration
    acc = forces / mass
    pos_new = pos + vel * dt + 0.5 * acc * dt**2
    pos_new = np.mod(pos_new, a)            # wrap back into the box [0, a)
    # recompute forces at the new positions
    forces_new = compute_forces(pos_new)
    acc_new = forces_new / mass
    # finish the velocity update with the average acceleration
    vel_new = vel + 0.5 * (acc + acc_new) * dt
    return pos_new, vel_new, forces_new

def run(pos, vel, forces, n_steps):
    for _ in range(n_steps):
        pos, vel, forces = velocity_verlet_step(pos, vel, forces)
    return pos, vel, forces

# -----------------------------------------------------------------------------
# Initial conditions: 5x5 ordered grid
# -----------------------------------------------------------------------------
pos = np.zeros((N, 2))
k = 0
for ix in range(n_side):
    for iy in range(n_side):
        pos[k] = [(ix + 0.5) * spacing, (iy + 0.5) * spacing]
        k += 1
pos_initial = pos.copy()

# Deterministic velocities spread over (-0.05, 0.05) via a golden-ratio sequence.
# The fractional parts of i*phi fill (0,1) very uniformly (low discrepancy);
# we map them linearly onto (-0.05, 0.05) for both velocity components.
phi = (1.0 + np.sqrt(5.0)) / 2.0
vel = np.zeros((N, 2))
for i in range(N):
    fx = (( (2 * i + 1) * phi) % 1.0)       # golden-ratio fractional part
    fy = (( (2 * i + 2) * phi) % 1.0)
    vel[i, 0] = -0.05 + 0.10 * fx
    vel[i, 1] = -0.05 + 0.10 * fy
# remove any net drift so the centre of mass stays put
vel -= vel.mean(axis=0)

# -----------------------------------------------------------------------------
# Warm-up run to t = 10, then production run to t = 100
# -----------------------------------------------------------------------------
forces = compute_forces(pos)
n_warm = int(round(t_warm / dt))
n_prod = int(round((t_prod - t_warm) / dt))

# capture a few intermediate snapshots along the way
snap_times = [0.0, t_warm, 40.0, 70.0, t_prod]
snapshots = {0.0: pos.copy()}

pos, vel, forces = run(pos, vel, forces, n_warm)
snapshots[t_warm] = pos.copy()

# production, stopping to grab snapshots at t = 40, 70, 100
for target in [40.0, 70.0, t_prod]:
    prev = max(tt for tt in snapshots.keys())
    steps = int(round((target - prev) / dt))
    pos, vel, forces = run(pos, vel, forces, steps)
    snapshots[target] = pos.copy()

pos_final = snapshots[t_prod]

# -----------------------------------------------------------------------------
# Structural check: regularity of nearest-neighbour distances.
# On a regular grid every particle has the same nearest-neighbour distance
# (= spacing), so the spread (std) of nearest-neighbour distances is ~0.
# In a disordered liquid the distances scatter, so the std is clearly nonzero.
# -----------------------------------------------------------------------------
def nearest_neighbour_distances(p):
    nn = np.zeros(N)
    for i in range(N):
        best = np.inf
        for j in range(N):
            if i == j:
                continue
            d = p[i] - p[j]
            d -= a * np.round(d / a)        # minimum image
            r = np.hypot(d[0], d[1])
            best = min(best, r)
        nn[i] = best
    return nn

nn_initial = nearest_neighbour_distances(pos_initial)
nn_final = nearest_neighbour_distances(pos_final)

std_initial = nn_initial.std()
std_final = nn_final.std()
mean_initial = nn_initial.mean()
mean_final = nn_final.mean()

# -----------------------------------------------------------------------------
# Snapshots figure
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, len(snap_times), figsize=(4 * len(snap_times), 4.2))
for ax, tt in zip(axes, snap_times):
    p = snapshots[tt]
    ax.scatter(p[:, 0], p[:, 1], s=60, c="tab:blue", edgecolors="k")
    ax.set_xlim(0, a)
    ax.set_ylim(0, a)
    ax.set_aspect("equal")
    ax.set_title(f"t = {tt:g}")
    ax.set_xlabel("x")
axes[0].set_ylabel("y")
fig.suptitle("2D Lennard-Jones: ordered grid -> disordered liquid")
fig.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.5.1_s4.png")

# -----------------------------------------------------------------------------
# Numerical results
# -----------------------------------------------------------------------------
print(f"Number of particles N: {N}")
print(f"Box side a: {a}")
print(f"Grid spacing: {spacing}")
print(f"Time step dt: {dt}")
print(f"Warm-up time: {t_warm}")
print(f"Production time: {t_prod}")
print(f"Initial mean nearest-neighbour distance: {mean_initial:.6f}")
print(f"Final   mean nearest-neighbour distance: {mean_final:.6f}")
print(f"Initial std of nearest-neighbour distances (regular grid): {std_initial:.6e}")
print(f"Final   std of nearest-neighbour distances (liquid):       {std_final:.6e}")
print(f"Ratio std_final / std_initial: {std_final / (std_initial + 1e-30):.6e}")
print(f"Starts on regular grid (std ~ 0): {std_initial < 1e-9}")
print(f"Ends disordered (std clearly > 0): {std_final > 1e-3}")
# Explanation:
print("Check explanation: a near-zero spread of nearest-neighbour distances "
      "means every particle is identically spaced (a regular grid), while a "
      "clearly nonzero spread means the spacings scatter (a disordered liquid), "
      "so this order-independent structural statistic confirms the ordered->liquid "
      "transition even though the exact chaotic final configuration is not reproducible.")
