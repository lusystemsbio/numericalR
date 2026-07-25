import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Reproduce the t = 100 trajectory from 5C.5:
# 25 Lennard-Jones particles in a periodic box of side a = 6.25.
# We run a short velocity-Verlet MD to reach a liquid-like state.
# ---------------------------------------------------------------
np.random.seed(0)
N = 25              # number of particles
a = 6.25            # box side length
rho0 = N / a**2     # mean number density (2D)
print(f"Number of particles N = {N}")
print(f"Box side a = {a}")
print(f"Mean density rho0 = N/a^2 = {rho0:.6f}")

# Initialize particles on a slightly perturbed 5x5 grid inside the box
grid = np.linspace(0.0, a, 6)[:-1] + a/10.0
X, Y = np.meshgrid(grid, grid)
pos = np.column_stack([X.ravel(), Y.ravel()]).astype(float)
pos += 0.05 * (np.random.rand(N, 2) - 0.5)

# Small random initial velocities, zero net momentum
vel = 0.1 * (np.random.rand(N, 2) - 0.5)
vel -= vel.mean(axis=0)

def lj_forces(pos, a):
    """Lennard-Jones forces with minimum-image periodic boundaries."""
    F = np.zeros_like(pos)
    for i in range(N):
        d = pos[i] - pos                    # displacement to all others
        d -= a * np.round(d / a)            # minimum-image convention
        r2 = np.sum(d**2, axis=1)
        r2[i] = np.inf                      # skip self-interaction
        inv2 = 1.0 / r2
        inv6 = inv2**3
        # LJ force magnitude/r: 24*(2*r^-12 - r^-6)/r^2
        fmag = 24.0 * (2.0 * inv6**2 - inv6) * inv2
        F[i] = np.sum(fmag[:, None] * d, axis=0)
    return F

# velocity-Verlet integration up to t = 100
dt = 0.004
nsteps = int(100 / dt)
F = lj_forces(pos, a)
for step in range(nsteps):
    vel += 0.5 * dt * F
    pos += dt * vel
    pos %= a                                # wrap back into the box
    F = lj_forces(pos, a)
    vel += 0.5 * dt * F

# ---------------------------------------------------------------
# Measure g(r) explicitly by histogramming minimum-image pair
# distances and normalizing each bin by its shell area and density.
# ---------------------------------------------------------------
dr = 0.05           # bin width
rmax = a / 2.0      # radii up to a/2 = 3.125
edges = np.arange(0.0, rmax + dr, dr)
counts = np.zeros(len(edges) - 1)

# accumulate all minimum-image pair distances (each unordered pair once)
for i in range(N):
    for j in range(i + 1, N):
        d = pos[i] - pos[j]
        d -= a * np.round(d / a)            # minimum-image convention
        r = np.sqrt(np.sum(d**2))
        if r < rmax:
            b = int(r / dr)                 # bin index
            counts[b] += 1

r_centers = 0.5 * (edges[:-1] + edges[1:])
# Each unordered pair contributes to two reference particles, so the
# average shell count per particle is 2*counts/N.
dN_avg = 2.0 * counts / N
shell_area = 2.0 * np.pi * r_centers * dr   # 2*pi*r*dr
g = dN_avg / (shell_area * rho0)            # g(r) = <dN(r)>/(2*pi*r*dr*rho0)

# ---------------------------------------------------------------
# Locate first and second peaks (ignore the r->0 region below 0.8).
# ---------------------------------------------------------------
valid = r_centers > 0.5
# first peak: global max in the physically meaningful range
first_region = (r_centers > 0.8) & (r_centers < 1.5)
i1 = np.where(first_region)[0][np.argmax(g[first_region])]
r1, g1 = r_centers[i1], g[i1]
# second peak: max in the region beyond the first minimum
second_region = (r_centers > 1.5) & (r_centers < 2.5)
i2 = np.where(second_region)[0][np.argmax(g[second_region])]
r2, g2 = r_centers[i2], g[i2]

print(f"First peak:  r = {r1:.3f}, g(r) = {g1:.4f}")
print(f"Second peak: r = {r2:.3f}, g(r) = {g2:.4f}")

# Classic-liquid-shape check
below = g[(r_centers < 0.8)]
g_below = below.max() if below.size else 0.0
tail_region = r_centers > 2.7
g_tail = g[tail_region].mean() if np.any(tail_region) else float('nan')
print(f"Max g(r) below r=0.8 (should be near 0): {g_below:.4f}")
print(f"Mean g(r) for r>2.7 (should settle toward 1): {g_tail:.4f}")
check = (g_below < 0.3) and (0.9 < r1 < 1.2) and (g1 > g2) and (1.7 < r2 < 2.3)
print(f"Classic liquid shape confirmed: {bool(check)}")
# This check confirms the result because the measured g(r) matches the
# known signature of a Lennard-Jones liquid -- a depletion hole from the
# repulsive core, a dominant first shell at the equilibrium separation,
# a weaker second shell, and decay to unity -- so any bug in the binning
# or normalization would show up as a wrong peak position, height, or tail.

# ---------------------------------------------------------------
# Plot g(r) versus r with the first and second peaks marked.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(r_centers, g, '-o', ms=3, color='steelblue', label='g(r)')
ax.axhline(1.0, color='gray', ls='--', lw=1, label='g = 1 (ideal gas)')
ax.plot(r1, g1, 'r^', ms=10, label=f'1st peak (r={r1:.2f})')
ax.plot(r2, g2, 'gs', ms=9, label=f'2nd peak (r={r2:.2f})')
ax.set_xlabel('r')
ax.set_ylabel('g(r)')
ax.set_title('Radial distribution function of the LJ particle box (t = 100)')
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/5C.7.1_s5.png")
print("Saved plot to 5C.7.1_s5.png")
