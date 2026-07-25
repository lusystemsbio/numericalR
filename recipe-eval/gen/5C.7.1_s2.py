import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Reproduce the kind of trajectory used in 5C.5: a 2D Lennard-Jones fluid
# of N = 25 particles in a periodic square box of side a = 6.25, integrated
# with velocity-Verlet up to t = 100. We then measure g(r) on the sampled
# configurations.  (sigma = epsilon = mass = 1, so the LJ minimum sits at
# r = 2^(1/6) ~ 1.12 and the "contact" separation is ~1.)
# ----------------------------------------------------------------------
np.random.seed(0)

N = 25
a = 6.25
rho0 = N / a**2          # mean number density N/a^2
dt = 0.005
t_final = 100.0
n_steps = int(t_final / dt)
T_set = 1.0              # target temperature (reduced units)

# --- initial positions on a 5x5 grid, small random velocities ---------
grid = np.linspace(0.5 * a / 5, a - 0.5 * a / 5, 5)
gx, gy = np.meshgrid(grid, grid)
pos = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
vel = np.random.normal(size=(N, 2))
vel -= vel.mean(axis=0)                       # zero net momentum
vel *= np.sqrt(T_set / (0.5 * (vel**2).sum() / N))  # scale to target T

def forces(p):
    # minimum-image pair forces for the LJ potential in a periodic box
    d = p[:, None, :] - p[None, :, :]         # all displacement vectors
    d -= a * np.round(d / a)                  # minimum image convention
    r2 = (d**2).sum(axis=2)
    np.fill_diagonal(r2, np.inf)              # skip self term
    inv2 = 1.0 / r2
    inv6 = inv2**3
    # LJ force magnitude/r: 48/r^14 - 24/r^8  = (24/r^2)(2 r^-12 - r^-6)
    fmag = 24.0 * inv2 * (2.0 * inv6**2 - inv6)
    f = (fmag[:, :, None] * d).sum(axis=1)
    return f

f = forces(pos)

# --- histogram setup --------------------------------------------------
dr = 0.05
r_max = a / 2.0          # 3.125 (largest meaningful minimum-image radius)
edges = np.arange(0.0, r_max + dr, dr)
centers = 0.5 * (edges[:-1] + edges[1:])
hist = np.zeros(len(centers))
n_frames = 0
equil_steps = n_steps // 5    # discard early transient before sampling

# --- velocity-Verlet integration, accumulating pair-distance histogram -
for step in range(n_steps):
    vel += 0.5 * dt * f
    pos += dt * vel
    pos %= a                                  # wrap into the box
    f = forces(pos)
    vel += 0.5 * dt * f

    # sample the trajectory (after equilibration) every 20 steps
    if step >= equil_steps and step % 20 == 0:
        d = pos[:, None, :] - pos[None, :, :]
        d -= a * np.round(d / a)              # minimum-image displacements
        rij = np.sqrt((d**2).sum(axis=2))
        iu = np.triu_indices(N, k=1)          # each unordered pair once
        h, _ = np.histogram(rij[iu], bins=edges)
        hist += h
        n_frames += 1

# --- normalize: each bin count -> average particles per shell, then by
#     shell area 2*pi*r*dr and mean density rho0 ------------------------
# each unordered pair contributes to 2 reference particles; averaged over
# N reference particles and n_frames sampled frames:
dN = 2.0 * hist / (N * n_frames)
shell_area = 2.0 * np.pi * centers * dr
g = dN / (shell_area * rho0)

# --- locate the first two peaks --------------------------------------
is_peak = (g[1:-1] > g[:-2]) & (g[1:-1] > g[2:])
peak_idx = np.where(is_peak)[0] + 1
peak_idx = peak_idx[g[peak_idx] > 1.05]       # ignore tiny ripples
first_peak_r = centers[peak_idx[0]]
first_peak_g = g[peak_idx[0]]
second_peak_r = centers[peak_idx[1]]
second_peak_g = g[peak_idx[1]]

# --- classic-liquid-shape check --------------------------------------
low_r_mask = centers < 0.8
near_zero_below_08 = np.max(g[low_r_mask])     # should be ~0
large_r_mask = centers > 2.6
tail_mean = np.mean(g[large_r_mask])           # should settle toward 1
check = (near_zero_below_08 < 0.2 and
         0.8 < first_peak_r < 1.3 and first_peak_g > 1.5 and
         1.6 < second_peak_r < 2.4 and
         0.7 < tail_mean < 1.3)

# --- report -----------------------------------------------------------
print(f"N                         : {N}")
print(f"box side a                : {a}")
print(f"mean density rho0 = N/a^2  : {rho0}")
print(f"bin width dr              : {dr}")
print(f"max radius a/2            : {r_max}")
print(f"sampled frames            : {n_frames}")
print(f"first peak  r             : {first_peak_r:.4f}")
print(f"first peak  g(r)          : {first_peak_g:.4f}")
print(f"second peak r             : {second_peak_r:.4f}")
print(f"second peak g(r)          : {second_peak_g:.4f}")
print(f"max g(r) for r < 0.8      : {near_zero_below_08:.4f}")
print(f"mean g(r) for r > 2.6     : {tail_mean:.4f}")
print(f"classic liquid shape?     : {check}")
print("Why this check confirms the result: reproducing the near-zero "
      "excluded region below r~0.8, a sharp first peak near the LJ "
      "equilibrium separation r~1, a weaker second-neighbor peak near "
      "r~2, and g(r)->1 at large r is the textbook signature of a liquid, "
      "so matching it validates that the histogramming and shell-area/"
      "density normalization were done correctly.")

# --- plot -------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(centers, g, '-', color='tab:blue', label='g(r)')
plt.axhline(1.0, color='gray', ls='--', lw=1, label='g = 1 (ideal)')
plt.plot(first_peak_r, first_peak_g, 'ro')
plt.annotate(f"1st peak\nr={first_peak_r:.2f}",
             (first_peak_r, first_peak_g),
             textcoords="offset points", xytext=(10, -5))
plt.plot(second_peak_r, second_peak_g, 'go')
plt.annotate(f"2nd peak\nr={second_peak_r:.2f}",
             (second_peak_r, second_peak_g),
             textcoords="offset points", xytext=(10, 5))
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial distribution function (2D LJ box, N=25, a=6.25, t=100)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/5C.7.1_s2.png")
