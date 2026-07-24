import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Periodic-boundary folding rule -----------------------------------------
def bc_periodic(x, a):
    # Fold coordinate(s) back into [0, a): x - floor(x/a)*a
    return x - np.floor(x / a) * a

# --- Minimum-image separation -----------------------------------------------
def min_image(dx, a):
    # Shift each component into (-a/2, a/2] so we use the nearest periodic image
    return dx - np.round(dx / a) * a

# --- Pair force: F(r) = 1/r^13 - 1/r^7, projected along separation ----------
def forces(pos, a):
    n = len(pos)
    F = np.zeros_like(pos)
    for i in range(n):
        for j in range(i + 1, n):
            dvec = min_image(pos[i] - pos[j], a)   # minimum-image vector i<-j
            r = np.hypot(*dvec)
            if r < 1e-12:
                continue
            fmag = 1.0 / r**13 - 1.0 / r**7        # scalar LJ force magnitude
            fvec = fmag * (dvec / r)               # project along the separation
            F[i] += fvec
            F[j] -= fvec                           # Newton's third law
    return F

# --- One velocity-Verlet step (unit mass), written out explicitly -----------
def verlet_step(pos, vel, acc, dt, a):
    # 1) half-kick + drift for positions using current acceleration
    pos_new = pos + vel * dt + 0.5 * acc * dt**2
    # 2) apply periodic-boundary correction to POSITIONS only (not velocities)
    pos_new = bc_periodic(pos_new, a)
    # 3) recompute acceleration (= force, since mass = 1) at the new positions
    acc_new = forces(pos_new, a)
    # 4) full velocity update using average of old and new accelerations
    vel_new = vel + 0.5 * (acc + acc_new) * dt
    return pos_new, vel_new, acc_new

# ============================================================================
# Run 5C.5: box side a = 6.25 with 25 particles on a 5x5 grid
# ============================================================================
a = 6.25
n = 25
grid = np.linspace(0.5, a - 0.5, 5)
pos = np.array([[x, y] for x in grid for y in grid], dtype=float)
vel = np.zeros_like(pos)
acc = forces(pos, a)

dt = 1e-3
nsteps = 200
for _ in range(nsteps):
    pos, vel, acc = verlet_step(pos, vel, acc, dt, a)

print("Box side a:", a)
print("Number of particles:", n)
print("Steps taken:", nsteps)
print("All final positions inside [0,a)?:", bool(np.all((pos >= 0) & (pos < a))))
print("Min final coordinate:", float(pos.min()))
print("Max final coordinate:", float(pos.max()))

# ============================================================================
# CHECK 1: a single particle stepping past a wall reappears on the opposite side
# ============================================================================
x_before = np.array([a - 0.05, 3.0])      # just inside the right wall
v = np.array([5.0, 0.0])                   # moving right, fast enough to exit
x_stepped = x_before + v * 0.05            # raw drift steps past x = a
x_after = bc_periodic(x_stepped, a)        # folding rule applied to position
print("Check1 raw x (before fold):", float(x_stepped[0]))
print("Check1 folded x (after fold):", float(x_after[0]))
print("Check1 reappeared near opposite (left) wall?:", bool(x_after[0] < 1.0))

# ============================================================================
# CHECK 2: particle leaves one side, reappears on opposite side; velocity kept
# ============================================================================
p = np.array([[0.02, a / 2]])              # near the LEFT wall
vp = np.array([[-3.0, 0.0]])               # moving LEFT, will exit through x=0
ap = np.zeros_like(p)                      # isolate PBC behaviour: no force
v_in = vp.copy()
p, vp, ap = verlet_step(p, vp, ap, 0.05, a)
print("Check2 position after crossing left wall:", p[0].tolist())
print("Check2 reappeared near opposite (right) wall?:", bool(p[0, 0] > a - 1.0))
print("Check2 velocity unchanged (never wrapped)?:", bool(np.allclose(vp, v_in)))

# --- Visualization ----------------------------------------------------------
plt.figure(figsize=(5, 5))
plt.scatter(pos[:, 0], pos[:, 1], c="tab:blue", label="final particle positions")
plt.gca().add_patch(plt.Rectangle((0, 0), a, a, fill=False, edgecolor="k"))
plt.xlim(-0.5, a + 0.5)
plt.ylim(-0.5, a + 0.5)
plt.gca().set_aspect("equal")
plt.title("Velocity-Verlet, periodic box (a=6.25, N=25)")
plt.legend(loc="upper right", fontsize=8)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.4.1_s5.png")

# Explanation: the check confirms the result because a position pushed just past
# a wall (x >= a or x < 0) is mapped by bc_periodic to the congruent point one
# box-length away, so it emerges on the opposite face while its velocity, which
# is never passed through the fold, is left exactly unchanged.
print("Explanation: folding maps a coordinate crossing one wall to its image on the opposite wall while leaving velocity untouched, which is exactly what periodic boundaries require.")
