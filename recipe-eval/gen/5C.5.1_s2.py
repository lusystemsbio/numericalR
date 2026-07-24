import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 2D Lennard-Jones box: ordered grid -> disordered liquid
# ----------------------------------------------------------------------

# ---- System parameters ----
N = 25                 # number of particles
side = 5               # 5x5 grid
a = 6.25               # box side length (periodic)
spacing = 1.25         # grid spacing (a / side)
dt = 0.01              # time step
mass = 1.0             # unit mass

# ---- Lennard-Jones force magnitude along the separation direction ----
# F(r) = 1/r^13 - 1/r^7  (given in the problem)
def lj_force_mag(r):
    return 1.0 / r**13 - 1.0 / r**7

# ---- Minimum-image separation under periodic boundaries ----
def min_image(delta, box):
    # wrap each component into (-box/2, box/2]
    return delta - box * np.round(delta / box)

# ---- Total force on every particle (pairwise, min-image) ----
def compute_forces(pos, box):
    f = np.zeros_like(pos)
    for i in range(N):
        for j in range(i + 1, N):
            dr = min_image(pos[i] - pos[j], box)   # vector i<-j
            r = np.hypot(dr[0], dr[1])
            fmag = lj_force_mag(r)                  # magnitude
            fvec = fmag * dr / r                    # directed force on i
            f[i] += fvec
            f[j] -= fvec                            # Newton's third law
    return f

# ---- Velocity-Verlet step with periodic wrap (from 5C.4) ----
def velocity_verlet_step(pos, vel, force, box):
    # half-kick + drift
    vel = vel + 0.5 * dt * force / mass
    pos = pos + dt * vel
    pos = np.mod(pos, box)                          # wrap into [0, box)
    # recompute force at new positions, then finish the kick
    new_force = compute_forces(pos, box)
    vel = vel + 0.5 * dt * new_force / mass
    return pos, vel, new_force

# ---- Initial ordered configuration: 5x5 grid at spacing 1.25 ----
pos0 = np.zeros((N, 2))
k = 0
for ix in range(side):
    for iy in range(side):
        pos0[k] = [(ix + 0.5) * spacing, (iy + 0.5) * spacing]
        k += 1
pos0 = np.mod(pos0, a)

# ---- Deterministic initial velocities via golden-ratio sequence ----
# fractional parts of i*phi are equidistributed in (0,1); map to (-0.05, 0.05)
phi = (1.0 + np.sqrt(5.0)) / 2.0
vel0 = np.zeros((N, 2))
for i in range(N):
    fx = ((i + 1) * phi) % 1.0
    fy = ((i + 1) * phi * phi) % 1.0
    vel0[i, 0] = (fx - 0.5) * 0.10   # spread over (-0.05, 0.05)
    vel0[i, 1] = (fy - 0.5) * 0.10
vel0 -= vel0.mean(axis=0)            # remove net drift

# ---- Nearest-neighbor distance spread (structural order metric) ----
def nn_stats(pos, box):
    nn = np.empty(N)
    for i in range(N):
        dmin = np.inf
        for j in range(N):
            if i == j:
                continue
            dr = min_image(pos[i] - pos[j], box)
            d = np.hypot(dr[0], dr[1])
            if d < dmin:
                dmin = d
        nn[i] = dmin
    return nn.mean(), nn.std()

# ----------------------------------------------------------------------
# Warm-up run to t = 10
# ----------------------------------------------------------------------
pos = pos0.copy()
vel = vel0.copy()
force = compute_forces(pos, a)
n_warm = int(round(10.0 / dt))
for _ in range(n_warm):
    pos, vel, force = velocity_verlet_step(pos, vel, force, a)
pos_warm = pos.copy()

# ----------------------------------------------------------------------
# Production run to t = 100, saving snapshots along the way
# ----------------------------------------------------------------------
n_prod = int(round(90.0 / dt))          # t = 10 -> t = 100
snap_times = [10, 30, 55, 100]          # production snapshot times
snap_steps = {int(round((t - 10) / dt)): t for t in snap_times}
snapshots = {}
snapshots[10] = pos.copy()

for step in range(1, n_prod + 1):
    pos, vel, force = velocity_verlet_step(pos, vel, force, a)
    if step in snap_steps:
        snapshots[snap_steps[step]] = pos.copy()
pos_final = pos.copy()

# ----------------------------------------------------------------------
# Structural check: ordered start vs. disordered end
# ----------------------------------------------------------------------
mean0, std0 = nn_stats(pos0, a)
meanf, stdf = nn_stats(pos_final, a)

print("Nearest-neighbor mean distance (initial grid): %.6f" % mean0)
print("Nearest-neighbor std  distance (initial grid): %.6e" % std0)
print("Nearest-neighbor mean distance (final liquid): %.6f" % meanf)
print("Nearest-neighbor std  distance (final liquid): %.6f" % stdf)
print("NN-distance std ratio (final/initial): %.3e" % (stdf / (std0 + 1e-30)))
print("Initial configuration ordered (NN std < 1e-6): %s" % (std0 < 1e-6))
print("Final configuration disordered (NN std > 0.01): %s" % (stdf > 0.01))

# ----------------------------------------------------------------------
# Snapshot figure: initial grid -> final disordered state
# ----------------------------------------------------------------------
panels = [("t = 0 (initial grid)", pos0)]
for t in snap_times:
    panels.append(("t = %d" % t, snapshots[t]))

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.ravel()
for ax, (title, p) in zip(axes, panels):
    ax.scatter(p[:, 0], p[:, 1], s=60, c="tab:blue", edgecolors="k")
    ax.set_xlim(0, a)
    ax.set_ylim(0, a)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
# leftover panel off
for ax in axes[len(panels):]:
    ax.axis("off")
fig.suptitle("2D Lennard-Jones: ordered grid -> disordered liquid", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.5.1_s2.png")

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# The nearest-neighbor distance has zero spread on the perfect grid but a
# finite spread once the system melts, so a jump in that spread from ~0 to a
# sizeable value confirms the ordered->disordered (solid->liquid) transition
# independently of the exact, chaos-sensitive final coordinates.
print("Explanation: NN-distance spread is ~0 for the ordered grid and becomes "
      "finite in the liquid, so its growth confirms melting regardless of the "
      "chaotic, run-specific final positions.")
