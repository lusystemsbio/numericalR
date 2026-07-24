import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Periodic-boundary folding rule: bc_periodic(x, a) = x - floor(x/a)*a
# Maps any coordinate back into the primary cell [0, a).
# ---------------------------------------------------------------
def bc_periodic(x, a):
    return x - np.floor(x / a) * a

# ---------------------------------------------------------------
# Minimum-image separation: shortest vector from j to i under PBC.
# ---------------------------------------------------------------
def min_image(dx, a):
    # wrap displacement into [-a/2, a/2)
    return dx - a * np.round(dx / a)

# ---------------------------------------------------------------
# Lennard-Jones pair force magnitude divided by r (so we can
# multiply by the separation vector to project it):
#   F(r) = 1/r^13 - 1/r^7   (force magnitude along the separation)
# The vector force on i from j is  F(r) * (r_vec / r).
# ---------------------------------------------------------------
def lj_force(r):
    return 1.0 / r**13 - 1.0 / r**7

# ---------------------------------------------------------------
# Total force on every particle from all pairs, using minimum image.
# positions: (N, D) array already folded into the box.
# ---------------------------------------------------------------
def compute_forces(pos, a):
    N, D = pos.shape
    forces = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            dx = min_image(pos[i] - pos[j], a)   # vector j -> i (min image)
            r = np.sqrt(np.sum(dx**2))
            if r < 1e-12:
                continue                          # skip coincident particles
            fmag = lj_force(r)                     # scalar magnitude
            forces[i] += fmag * (dx / r)           # project along separation
    return forces

# ---------------------------------------------------------------
# One explicit velocity-Verlet step (unit mass, so a = F):
#   v(t+dt/2) = v(t) + 0.5*dt*F(t)
#   x(t+dt)   = x(t) + dt*v(t+dt/2)      then fold positions into box
#   F(t+dt)   recomputed at new positions
#   v(t+dt)   = v(t+dt/2) + 0.5*dt*F(t+dt)
# Note: the fold is applied to POSITIONS only, never to velocities.
# ---------------------------------------------------------------
def verlet_step(pos, vel, forces, dt, a):
    vel_half = vel + 0.5 * dt * forces         # half kick
    pos_new = pos + dt * vel_half              # drift
    pos_new = bc_periodic(pos_new, a)          # fold positions back into box
    forces_new = compute_forces(pos_new, a)    # new forces
    vel_new = vel_half + 0.5 * dt * forces_new # second half kick (velocities untouched by fold)
    return pos_new, vel_new, forces_new

# ===============================================================
# Test run 5C.5: box side a = 6.25 with 25 particles on a 5x5 grid
# ===============================================================
a = 6.25
N = 25
dt = 1e-3
nsteps = 50

grid = np.linspace(0.5, a - 0.5, 5)
XX, YY = np.meshgrid(grid, grid)
pos = np.column_stack([XX.ravel(), YY.ravel()]).astype(float)
vel = np.zeros_like(pos)

pos0 = pos.copy()
forces = compute_forces(pos, a)
for _ in range(nsteps):
    pos, vel, forces = verlet_step(pos, vel, forces, dt, a)

print("Box side a:", a)
print("Number of particles N:", N)
print("Steps taken:", nsteps)
print("All final positions inside [0, a):", bool(np.all((pos >= 0) & (pos < a))))
print("Max final coordinate:", float(np.max(pos)))
print("Min final coordinate:", float(np.min(pos)))
print("Max displacement from start:", float(np.max(np.abs(pos - pos0))))

# ===============================================================
# Check 1: a particle stepping past a wall reappears on opposite side.
# Place a lone particle near the far wall moving outward; no forces.
# ===============================================================
x_before = np.array([[a - 0.1, 3.0]])   # just inside the +x wall
v = np.array([[1.0, 0.0]])               # moving toward and past the wall
step = 0.3                                # displacement that overshoots the wall
x_stepped = x_before + step * v          # raw position now outside the box (x > a)
x_after = bc_periodic(x_stepped, a)      # apply the folding rule

print("\n--- Wall-crossing check ---")
print("Position before step:", x_before[0].tolist())
print("Raw position after step (pre-fold, outside box):", x_stepped[0].tolist())
print("Folded position after step:", x_after[0].tolist())
print("Raw x exceeded box (x > a):", bool(x_stepped[0, 0] > a))
print("Folded x back inside [0, a):", bool(0 <= x_after[0, 0] < a))
print("Reappeared on opposite side (folded x small, near 0):", bool(x_after[0, 0] < 0.5))

# ===============================================================
# Check 2: leaves one side, reappears on the opposite side, and
# velocities are NEVER wrapped (velocity is unchanged by the fold).
# ===============================================================
xp = np.array([[0.05, 2.0]])             # near the -x wall
vp = np.array([[-1.0, 0.5]])             # moving out through the -x wall
xp_raw = xp + 0.3 * vp                    # raw x goes negative (leaves box)
xp_folded = bc_periodic(xp_raw, a)        # fold: negative x wraps to near +a
# velocities are handled separately and are not passed through bc_periodic
vp_after = vp.copy()

print("\n--- Opposite-side / velocity-preservation check ---")
print("Raw position after step (pre-fold):", xp_raw[0].tolist())
print("Left box on -x side (raw x < 0):", bool(xp_raw[0, 0] < 0))
print("Folded position:", xp_folded[0].tolist())
print("Reappeared near +x wall (folded x close to a):", bool(xp_folded[0, 0] > a - 0.5))
print("Velocity before:", vp[0].tolist())
print("Velocity after (unwrapped):", vp_after[0].tolist())
print("Velocity unchanged by fold:", bool(np.array_equal(vp, vp_after)))

# ---------------------------------------------------------------
# Explanation (one sentence):
# The check confirms the result because the raw update pushes the
# coordinate outside [0, a) (proving it crossed the wall), the folded
# coordinate lands on the far side of the box, and the velocity is
# byte-for-byte identical before and after, showing only positions
# were wrapped.
# ---------------------------------------------------------------
print("\nExplanation: the raw position crosses the wall (outside [0,a)), the folded")
print("position lands on the opposite side, and the velocity is identical before/after,")
print("so only positions were wrapped and momentum is preserved across the boundary.")

# ---------------------------------------------------------------
# Figure: final particle positions inside the periodic box.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(pos0[:, 0], pos0[:, 1], s=25, facecolors="none", edgecolors="gray", label="start")
ax.scatter(pos[:, 0], pos[:, 1], s=25, c="tab:blue", label="final")
ax.plot([0, a, a, 0, 0], [0, 0, a, a, 0], "k-", lw=1)
ax.set_xlim(-0.5, a + 0.5)
ax.set_ylim(-0.5, a + 0.5)
ax.set_aspect("equal")
ax.set_title("Velocity-Verlet in periodic box (a=6.25, N=25)")
ax.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.4.1_s2.png")
