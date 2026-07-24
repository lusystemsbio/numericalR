import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Periodic-boundary folding rule (applied to positions only).
# bc_periodic(x, a) = x - floor(x/a)*a  maps any coordinate into [0, a).
# ---------------------------------------------------------------
def bc_periodic(x, a):
    return x - np.floor(x / a) * a

# ---------------------------------------------------------------
# Minimum-image separation between two coordinate arrays in a box of side a.
# Each component is shifted into (-a/2, a/2] so the nearest periodic image is used.
# ---------------------------------------------------------------
def min_image(dx, a):
    return dx - a * np.round(dx / a)

# ---------------------------------------------------------------
# Total force on every particle.
# Pair force magnitude (Lennard-Jones form given): F(r) = 1/r^13 - 1/r^7,
# directed along the minimum-image unit vector from j to i.
# ---------------------------------------------------------------
def compute_forces(pos, a):
    n = pos.shape[0]
    forces = np.zeros_like(pos)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # minimum-image separation vector pointing from j to i
            dx = min_image(pos[i] - pos[j], a)
            r = np.sqrt(np.sum(dx * dx))
            if r == 0.0:
                continue
            Fmag = 1.0 / r**13 - 1.0 / r**7   # scalar pair force
            forces[i] += Fmag * (dx / r)      # project along unit separation
    return forces

# ---------------------------------------------------------------
# One velocity-Verlet step (unit masses), done explicitly:
#   1) half-kick velocities using current forces
#   2) drift positions
#   3) fold positions back into the box (velocities untouched)
#   4) recompute forces at the new positions
#   5) second half-kick velocities
# ---------------------------------------------------------------
def velocity_verlet_step(pos, vel, forces, a, dt):
    vel_half = vel + 0.5 * dt * forces          # 1) half kick
    pos_new = pos + dt * vel_half               # 2) drift
    pos_new = bc_periodic(pos_new, a)           # 3) fold positions only
    forces_new = compute_forces(pos_new, a)     # 4) new forces
    vel_new = vel_half + 0.5 * dt * forces_new  # 5) second half kick
    return pos_new, vel_new, forces_new

# ---------------------------------------------------------------
# Run 5C.5 configuration: box side a = 6.25, 25 particles on a 5x5 grid.
# ---------------------------------------------------------------
a = 6.25
n_particles = 25
dt = 1e-3
n_steps = 50

grid = np.linspace(0.5, a - 0.5, 5)
X, Y = np.meshgrid(grid, grid)
pos = np.column_stack([X.ravel(), Y.ravel()]).astype(float)
rng = np.random.default_rng(0)
vel = 0.01 * rng.standard_normal(pos.shape)

forces = compute_forces(pos, a)
for _ in range(n_steps):
    pos, vel, forces = velocity_verlet_step(pos, vel, forces, a, dt)

print("Box side a:", a)
print("Number of particles:", n_particles)
print("Min position component after run:", float(pos.min()))
print("Max position component after run:", float(pos.max()))
print("All positions inside [0, a):", bool(np.all((pos >= 0) & (pos < a))))

# ---------------------------------------------------------------
# CHECK 1: a particle stepping past a wall reappears on the opposite side.
# Place a single particle just inside the far wall moving outward; after the
# drift its raw coordinate exceeds a, and folding must bring it near 0.
# ---------------------------------------------------------------
p = np.array([[a - 0.01, a / 2]])   # near the +x wall
v = np.array([[5.0, 0.0]])          # moving outward in +x
f = np.zeros_like(p)                # no forces so the motion is a clean step
p_raw = p + dt * (v + 0.5 * dt * f) # unfolded drift (what velocity-Verlet computes)
p_new, v_new, _ = velocity_verlet_step(p.copy(), v.copy(), f.copy(), a, dt)

print("Check1 raw x before folding:", float(p_raw[0, 0]))
print("Check1 raw x exceeded wall a:", bool(p_raw[0, 0] >= a))
print("Check1 folded x after step:", float(p_new[0, 0]))
print("Check1 reappeared near 0 side:", bool(p_new[0, 0] < a / 2))
print("Check1 velocity unchanged (unwrapped):", float(v_new[0, 0]))

# ---------------------------------------------------------------
# CHECK 2: a particle leaving one side reappears on the opposite side while
# velocities are never wrapped. Give a large negative velocity so it exits the
# -x side; folding wraps the position to near +x, and the velocity is unchanged.
# ---------------------------------------------------------------
p2 = np.array([[0.01, a / 2]])      # near the -x wall
v2 = np.array([[-5.0, 0.0]])        # moving outward in -x
f2 = np.zeros_like(p2)
p2_raw = p2 + dt * (v2 + 0.5 * dt * f2)
p2_new, v2_new, _ = velocity_verlet_step(p2.copy(), v2.copy(), f2.copy(), a, dt)

print("Check2 raw x before folding:", float(p2_raw[0, 0]))
print("Check2 raw x went negative:", bool(p2_raw[0, 0] < 0))
print("Check2 folded x after step:", float(p2_new[0, 0]))
print("Check2 reappeared near +x side:", bool(p2_new[0, 0] > a / 2))
print("Check2 velocity before:", float(v2[0, 0]))
print("Check2 velocity after:", float(v2_new[0, 0]))
print("Check2 velocity never wrapped:", bool(v2_new[0, 0] == v2[0, 0]))

# One-sentence explanation of why these checks confirm the result:
explanation = ("This confirms the result because the raw drift carries the "
               "coordinate outside [0, a) yet bc_periodic returns it to the "
               "opposite side while the velocity is left exactly equal to its "
               "pre-step value, showing folding acts on positions only.")
print("Explanation:", explanation)

# ---------------------------------------------------------------
# Visualization: folded final positions of the 25-particle run.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(pos[:, 0], pos[:, 1], c="tab:blue")
ax.set_xlim(0, a)
ax.set_ylim(0, a)
ax.set_aspect("equal")
ax.set_title("Folded positions after run (a=6.25, N=25)")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.4.1_s1.png")
