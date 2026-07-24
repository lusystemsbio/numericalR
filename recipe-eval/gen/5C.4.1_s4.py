import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Periodic-boundary folding rule (applied to positions only)
# bc_periodic(x, a) = x - floor(x/a)*a  maps any coordinate into [0, a).
# ---------------------------------------------------------------
def bc_periodic(x, a):
    return x - np.floor(x / a) * a

# ---------------------------------------------------------------
# Minimum-image separation between two coordinate arrays for box side a.
# The displacement is wrapped into (-a/2, a/2] so each particle interacts
# with the nearest periodic image of every other particle.
# ---------------------------------------------------------------
def min_image(dx, a):
    return dx - a * np.round(dx / a)

# ---------------------------------------------------------------
# Total force on every particle from the Lennard-Jones pair force
# F(r) = 1/r^13 - 1/r^7, projected along the minimum-image separation.
# forces[i] = sum_j F(r_ij) * (r_i - r_j)/r_ij   (unit-mass particles)
# ---------------------------------------------------------------
def compute_forces(pos, a):
    n = pos.shape[0]
    forces = np.zeros_like(pos)
    for i in range(n):
        for j in range(i + 1, n):
            dr = min_image(pos[i] - pos[j], a)      # minimum-image separation vector
            r = np.sqrt(np.sum(dr * dr))
            if r < 1e-12:                            # avoid singular overlap
                continue
            fmag = 1.0 / r**13 - 1.0 / r**7          # scalar LJ pair force F(r)
            fij = fmag * dr / r                      # projected along separation
            forces[i] += fij                         # Newton's third law
            forces[j] -= fij
    return forces

# ---------------------------------------------------------------
# One explicit velocity-Verlet step (unit mass), written out in stages:
#   1) half-kick velocities using current force
#   2) drift positions
#   3) FOLD positions back into the box (velocities are NOT wrapped)
#   4) recompute force at new positions
#   5) half-kick velocities using new force
# ---------------------------------------------------------------
def verlet_step(pos, vel, forces, dt, a):
    vel = vel + 0.5 * dt * forces                    # 1) half velocity update
    pos = pos + dt * vel                             # 2) position update
    pos = bc_periodic(pos, a)                        # 3) periodic correction on positions only
    new_forces = compute_forces(pos, a)             # 4) forces at new positions
    vel = vel + 0.5 * dt * new_forces                # 5) second half velocity update
    return pos, vel, new_forces

# ===============================================================
# 5C.5 run: box side a = 6.25 with 25 particles
# ===============================================================
a = 6.25
n = 25
dt = 1e-3
np.random.seed(0)

# Place 25 particles on a 5x5 grid inside the box, small random velocities.
grid = np.linspace(0.5, a - 0.5, 5)
pos = np.array([[x, y] for x in grid for y in grid], dtype=float)
vel = 0.1 * (np.random.rand(n, 2) - 0.5)

print(f"Box side a = {a}")
print(f"Number of particles = {n}")
print(f"Time step dt = {dt}")

forces = compute_forces(pos, a)
nsteps = 200
for _ in range(nsteps):
    pos, vel, forces = verlet_step(pos, vel, forces, dt, a)

print(f"Steps integrated = {nsteps}")
print(f"All positions inside [0, a)? {bool(np.all((pos >= 0) & (pos < a)))}")
print(f"Min position component after run = {pos.min():.6f}")
print(f"Max position component after run = {pos.max():.6f}")

# ===============================================================
# CHECK 1: a particle stepping past a wall reappears on the opposite side.
# Start one free particle just inside the right wall moving to the right;
# with no forces its position should exceed a, then fold to near 0.
# ===============================================================
xc = np.array([[a - 0.01]])          # position just inside the right wall
vc = np.array([[1.0]])               # velocity pointing outward (+x)
dtc = 0.05
xc_before_fold = xc + dtc * vc       # drift only (single free particle, no force)
xc_after_fold = bc_periodic(xc_before_fold, a)
print("\n--- CHECK 1: stepping past a wall ---")
print(f"Position before fold (past wall) = {xc_before_fold[0,0]:.6f}")
print(f"Position after  fold (wrapped)   = {xc_after_fold[0,0]:.6f}")
print(f"Wrapped to opposite side? {bool(xc_before_fold[0,0] >= a and xc_after_fold[0,0] < a)}")

# ===============================================================
# CHECK 2: particle leaves one side and reappears on the other while
# its velocity is never wrapped. Integrate a single free particle and
# confirm the velocity is unchanged by the fold.
# ===============================================================
xp = np.array([[a - 0.02]])
vp = np.array([[0.5]])
vp_initial = vp.copy()
fp = np.zeros_like(xp)               # free particle: force stays zero
crossed = False
for _ in range(500):
    xp, vp, fp = verlet_step(xp, vp, fp, dtc, a)
    if xp[0, 0] < a - 0.02:          # detect it wrapped back near the left edge
        crossed = True
print("\n--- CHECK 2: reappears opposite side, velocity not wrapped ---")
print(f"Initial velocity = {vp_initial[0,0]:.6f}")
print(f"Final velocity   = {vp[0,0]:.6f}")
print(f"Velocity unchanged by folding? {bool(np.isclose(vp[0,0], vp_initial[0,0]))}")
print(f"Final position inside [0, a)? {bool(0 <= xp[0,0] < a)}")
print(f"Particle crossed to opposite side? {crossed}")

# Explanation (one sentence):
print("\nWhy this confirms the result: because the fold subtracts an integer "
      "number of box lengths from the position only, a particle exiting one "
      "face re-enters the opposite face with its velocity exactly preserved, "
      "which is precisely what a correct periodic boundary must do.")

# ---------------------------------------------------------------
# Visual record of the 25-particle configuration and the wrap check.
# ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
ax1.scatter(pos[:, 0], pos[:, 1], c="tab:blue")
ax1.set_xlim(0, a); ax1.set_ylim(0, a)
ax1.set_title("25 particles after run (folded in box)")
ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.set_aspect("equal")

ax2.axvline(0, color="k", ls="--"); ax2.axvline(a, color="k", ls="--")
ax2.plot([xc_before_fold[0,0]], [0.5], "rx", ms=12, label="before fold (past wall)")
ax2.plot([xc_after_fold[0,0]], [0.5], "go", ms=10, label="after fold (opposite side)")
ax2.set_xlim(-0.5, a + 0.5); ax2.set_ylim(0, 1)
ax2.set_title("Wall-crossing wrap check"); ax2.set_xlabel("x"); ax2.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.4.1_s4.png")
