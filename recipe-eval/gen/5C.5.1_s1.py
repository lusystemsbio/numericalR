import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# 2D Lennard-Jones box: ordered grid -> disordered liquid
# 25 unit-mass particles, square box side a = 6.25,
# force F(r) = 1/r^13 - 1/r^7, minimum-image periodic boundaries,
# velocity-Verlet integrator.
# ---------------------------------------------------------------

N = 25            # number of particles
a = 6.25          # box side length
dt = 0.01         # time step
mass = 1.0        # unit mass

# --- initial positions: 5x5 grid at spacing 1.25 -----------------
spacing = 1.25
grid_1d = (np.arange(5) + 0.5) * spacing          # centered cells: 0.625,1.875,...
X, Y = np.meshgrid(grid_1d, grid_1d)
pos0 = np.column_stack([X.ravel(), Y.ravel()])    # (25,2) ordered grid
pos = pos0.copy()

# --- deterministic velocities via golden-ratio sequence ----------
# Spread 2*N = 50 values evenly in [0,1) using fractional parts of
# k*phi, then map linearly onto (-0.05, 0.05).
phi = (1.0 + np.sqrt(5.0)) / 2.0
k = np.arange(1, 2 * N + 1)
frac = (k * phi) % 1.0                             # low-discrepancy in [0,1)
vel = ((frac - 0.5) * 0.1).reshape(N, 2)           # -> (-0.05, 0.05)

def accelerations(p):
    """Acceleration on each particle from F(r)=1/r^13 - 1/r^7,
    summed over all pairs using the minimum-image convention."""
    acc = np.zeros_like(p)
    for i in range(N):
        for j in range(i + 1, N):
            d = p[i] - p[j]                        # separation vector
            d -= a * np.round(d / a)               # minimum image
            r = np.sqrt(d @ d)
            # scalar force magnitude (positive = repulsive)
            Fmag = 1.0 / r**13 - 1.0 / r**7
            fvec = Fmag * (d / r)                   # force on i from j
            acc[i] += fvec / mass                   # a = F/m
            acc[j] -= fvec / mass                   # Newton's third law
    return acc

def verlet_step(p, v, acc):
    """One velocity-Verlet step with periodic wrap-around."""
    p_new = p + v * dt + 0.5 * acc * dt**2          # position update
    p_new = np.mod(p_new, a)                        # periodic boundaries
    acc_new = accelerations(p_new)                  # force at new position
    v_new = v + 0.5 * (acc + acc_new) * dt          # velocity update
    return p_new, v_new, acc_new

def run(p, v, t_total):
    """Integrate for t_total, returning final state."""
    nsteps = int(round(t_total / dt))
    acc = accelerations(p)
    for _ in range(nsteps):
        p, v, acc = verlet_step(p, v, acc)
    return p, v

# --- warm-up run to t = 10 ---------------------------------------
pos, vel = run(pos, vel, 10.0)
pos_warm = pos.copy()

# --- production run to t = 100, capturing snapshots --------------
snap_times = [10.0, 30.0, 55.0, 100.0]             # times to record (from t=10)
snaps = {10.0: pos_warm.copy()}
prev_t = 10.0
for t in snap_times[1:]:
    pos, vel = run(pos, vel, t - prev_t)
    snaps[t] = pos.copy()
    prev_t = t
pos_final = pos.copy()

# ---------------------------------------------------------------
# Snapshots figure: initial grid -> final disordered state
# ---------------------------------------------------------------
panels = [("t = 0 (grid)", pos0)] + [(f"t = {t:g}", snaps[t]) for t in snap_times]
fig, axes = plt.subplots(2, 3, figsize=(11, 7.2))
axes = axes.ravel()
for ax, (title, p) in zip(axes, panels):
    ax.scatter(p[:, 0], p[:, 1], s=60, c="tab:blue", edgecolors="k")
    ax.set_xlim(0, a); ax.set_ylim(0, a)
    ax.set_aspect("equal"); ax.set_title(title)
axes[-1].axis("off")   # unused sixth panel
fig.suptitle("2D Lennard-Jones: ordered grid melting into a liquid")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.5.1_s1.png")

# ---------------------------------------------------------------
# Separate structural check: grid at start, liquid at end.
# We measure how far particles sit from the perfect grid, and the
# spread of nearest-neighbour distances (0 spread = ordered lattice,
# nonzero spread = disordered liquid). These structural averages are
# robust to chaotic divergence between R and Python.
# ---------------------------------------------------------------
def min_image_dist(p, q):
    d = p - q
    d -= a * np.round(d / a)
    return np.sqrt(d @ d)

def grid_deviation(p):
    """Max min-image distance of each particle from its grid site."""
    return max(min_image_dist(p[i], pos0[i]) for i in range(N))

def nn_stats(p):
    """Mean and std of nearest-neighbour distances (structural average)."""
    nn = []
    for i in range(N):
        dmin = min(min_image_dist(p[i], p[j]) for j in range(N) if j != i)
        nn.append(dmin)
    nn = np.array(nn)
    return nn.mean(), nn.std()

dev_start = grid_deviation(pos0)
dev_end = grid_deviation(pos_final)
nn_mean_start, nn_std_start = nn_stats(pos0)
nn_mean_end, nn_std_end = nn_stats(pos_final)

ordered_start = dev_start < 1e-9 and nn_std_start < 1e-9
disordered_end = dev_end > 0.3 and nn_std_end > 1e-3

print(f"Number of particles: {N}")
print(f"Box side a: {a}")
print(f"Initial grid spacing: {spacing}")
print(f"Max deviation from grid at start: {dev_start:.6e}")
print(f"Max deviation from grid at end:   {dev_end:.6f}")
print(f"Nearest-neighbour distance mean at start: {nn_mean_start:.6f}")
print(f"Nearest-neighbour distance std  at start: {nn_std_start:.6e}")
print(f"Nearest-neighbour distance mean at end:   {nn_mean_end:.6f}")
print(f"Nearest-neighbour distance std  at end:   {nn_std_end:.6f}")
print(f"Starts on regular grid (ordered): {ordered_start}")
print(f"Ends in irregular liquid arrangement (disordered): {disordered_end}")
print(f"CHECK PASSED: {ordered_start and disordered_end}")
# One-sentence explanation:
print("Explanation: the check confirms the result because a zero grid-deviation "
      "and zero nearest-neighbour spread at the start proves a perfect lattice, "
      "while a large grid-deviation and nonzero nearest-neighbour spread at the end "
      "proves a disordered liquid, and these structural averages are meaningful even "
      "though the chaotic trajectory itself is not reproducible between R and Python.")
