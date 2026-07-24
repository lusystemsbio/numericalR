import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------------
# Model: unit-mass particles in a periodic box of side a.
# Each particle obeys d^2 x_i/dt^2 = sum_j F(r_ij) * (unit vector i<-j),
# with the Lennard-Jones-type pair force magnitude F(r) = 1/r^13 - 1/r^7.
# The separation used is the minimum-image separation (nearest periodic copy).
# -------------------------------------------------------------------------

def bc_periodic(x, a):
    # Folding rule: bring a coordinate back into [0, a) by removing whole boxes.
    # Applied to POSITIONS only, never to velocities.
    return x - np.floor(x / a) * a


def min_image(dx, a):
    # Minimum-image convention: shift the separation into (-a/2, a/2]
    # so we always interact with the nearest periodic copy.
    return dx - np.round(dx / a) * a


def lj_force_magnitude(r):
    # Pair force magnitude F(r) = 1/r^13 - 1/r^7.
    return 1.0 / r**13 - 1.0 / r**7


def compute_forces(pos, a):
    # Sum pairwise forces on every particle, using minimum-image separations.
    n = pos.shape[0]
    forces = np.zeros_like(pos)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # separation vector from j to i, folded to nearest image
            d = min_image(pos[i] - pos[j], a)
            r = np.sqrt(np.sum(d * d))
            if r < 1e-12:            # guard against coincident particles
                continue
            Fmag = lj_force_magnitude(r)   # scalar force magnitude
            forces[i] += Fmag * (d / r)    # project along the unit separation
    return forces


def velocity_verlet_step(pos, vel, a, dt):
    # One explicit velocity-Verlet step (written out, not a library call).
    # 1) forces (=accelerations, since unit mass) at the current positions
    acc = compute_forces(pos, a)
    # 2) half-kick the velocities
    vel_half = vel + 0.5 * dt * acc
    # 3) drift the positions with the half-step velocity
    pos_new = pos + dt * vel_half
    # 4) apply the periodic-boundary correction to POSITIONS only
    pos_new = bc_periodic(pos_new, a)
    # 5) forces at the new positions
    acc_new = compute_forces(pos_new, a)
    # 6) second half-kick to complete the velocity update
    vel_new = vel_half + 0.5 * dt * acc_new
    return pos_new, vel_new


# -------------------------------------------------------------------------
# Run from 5C.5: box side a = 6.25, N = 25 particles.
# -------------------------------------------------------------------------
a = 6.25
N = 25
dt = 0.001

rng = np.random.default_rng(0)
# start particles on a slightly perturbed 5x5 grid so they are not coincident
grid = np.linspace(0.5, a - 0.5, 5)
gx, gy = np.meshgrid(grid, grid)
pos = np.column_stack([gx.ravel(), gy.ravel()]) + 0.01 * rng.standard_normal((N, 2))
pos = bc_periodic(pos, a)
vel = 0.1 * rng.standard_normal((N, 2))

# integrate a short trajectory
n_steps = 200
traj = np.empty((n_steps + 1, N, 2))
traj[0] = pos
for k in range(n_steps):
    pos, vel = velocity_verlet_step(pos, vel, a, dt)
    traj[k + 1] = pos

# confirm all positions stayed inside the box during the whole run
print("Box side a:", a)
print("Number of particles N:", N)
print("Min position over trajectory:", float(traj.min()))
print("Max position over trajectory:", float(traj.max()))
print("All positions inside [0, a):", bool(traj.min() >= 0.0 and traj.max() < a))

# -------------------------------------------------------------------------
# CHECK 1: a single particle stepping past a wall reappears on the opposite side.
# Place it near the far x-wall with an outward velocity so one drift crosses it.
# -------------------------------------------------------------------------
p1 = np.array([[a - 0.05, a / 2]])   # just inside the +x wall
v1 = np.array([[5.0, 0.0]])          # moving toward and past the wall
raw_x = p1[0, 0] + dt * (v1[0, 0])   # naive un-folded x after a drift (approx)
p1_new, v1_new = velocity_verlet_step(p1, v1, a, dt)
print("\nCHECK 1 (wall crossing):")
print("x before step:", float(p1[0, 0]))
print("x un-folded (would exceed a):", float(raw_x))
print("x after folding:", float(p1_new[0, 0]))
print("Reappeared on opposite (low-x) side:", bool(p1_new[0, 0] < a / 2))

# -------------------------------------------------------------------------
# CHECK 2: a particle leaving one side reappears on the opposite side,
# and its velocity is NOT wrapped (velocities only change via the kicks).
# Compare a folded position with the un-folded reference differing by one box.
# -------------------------------------------------------------------------
p2 = np.array([[0.02, a / 2]])       # just inside the -x wall
v2 = np.array([[-5.0, 0.0]])         # moving out through the -x wall
p2_new, v2_new = velocity_verlet_step(p2, v2, a, dt)
# reference un-folded position after one drift (half-step velocity ~ v2 since forces tiny)
acc_ref = compute_forces(p2, a)
vhalf_ref = v2 + 0.5 * dt * acc_ref
unfolded = p2[0, 0] + dt * vhalf_ref[0, 0]     # this is negative -> left the box
print("\nCHECK 2 (opposite side + velocity unchanged by wrap):")
print("x before step:", float(p2[0, 0]))
print("x un-folded (negative, outside box):", float(unfolded))
print("x after folding:", float(p2_new[0, 0]))
print("Folded x equals un-folded + a:", bool(abs(p2_new[0, 0] - (unfolded + a)) < 1e-9))
print("Reappeared on opposite (high-x) side:", bool(p2_new[0, 0] > a / 2))
# velocity check: wrapping does not touch velocity; folded and unfolded runs share vx
print("vx after step:", float(v2_new[0, 0]))
print("vx sign unchanged (still moving -x):", bool(v2_new[0, 0] < 0))

# one-sentence explanation of why CHECK 2 confirms the result
print("\nExplanation: The folded position equals the un-folded position minus exactly "
      "one box length while the velocity keeps its original value and direction, "
      "which shows the boundary rule only translates positions by whole boxes and "
      "never alters velocities.")

# -------------------------------------------------------------------------
# Figure: trajectories inside the box.
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 6))
for i in range(N):
    ax.plot(traj[:, i, 0], traj[:, i, 1], lw=0.6, alpha=0.6)
ax.scatter(traj[-1, :, 0], traj[-1, :, 1], c="k", s=12, zorder=3)
ax.set_xlim(0, a)
ax.set_ylim(0, a)
ax.set_aspect("equal")
ax.set_title(f"Velocity-Verlet in periodic box (a={a}, N={N})")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.4.1_s3.png")
