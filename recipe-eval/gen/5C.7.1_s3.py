import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Reproduce the "t = 100 trajectory from 5C.5":
# 25 Lennard-Jones particles in a 2D periodic box of side a = 6.25.
# We run a short velocity-Verlet MD to equilibrium and collect the
# late-time configurations that make up the trajectory near t = 100.
# ---------------------------------------------------------------
np.random.seed(0)

N = 25                       # number of particles
a = 6.25                     # box side
rho0 = N / a**2              # mean (number) density, rho0 = N/a^2
dt = 0.005                   # MD time step
t_final = 100.0             # integrate up to t = 100
nsteps = int(t_final / dt)
rc = 2.5                     # LJ cutoff
T_target = 1.0               # target temperature (liquid regime)

# --- initialise particles on a 5x5 grid, random velocities ---
grid = np.linspace(0, a, int(np.sqrt(N)), endpoint=False)
X, Y = np.meshgrid(grid, grid)
pos = np.column_stack([X.ravel(), Y.ravel()]).astype(float)
pos += 0.1 * (np.random.rand(N, 2) - 0.5)          # tiny displacement off lattice
vel = np.random.randn(N, 2)
vel -= vel.mean(axis=0)                              # zero net momentum

def forces(pos):
    # minimum-image pairwise LJ forces with cutoff rc
    d = pos[:, None, :] - pos[None, :, :]           # all displacement vectors
    d -= a * np.round(d / a)                          # minimum-image convention
    r2 = np.sum(d**2, axis=2)
    np.fill_diagonal(r2, np.inf)                     # ignore self-interaction
    mask = r2 < rc**2
    inv_r2 = np.where(mask, 1.0 / r2, 0.0)
    inv_r6 = inv_r2**3
    # LJ force magnitude/r^2 factor: 24*(2*r^-12 - r^-6)/r^2
    fmag = 24.0 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0)
    fmag = np.where(mask, fmag, 0.0)
    fvec = (fmag[:, :, None] * d).sum(axis=1)
    return fvec

# --- velocity Verlet with occasional velocity rescaling (thermostat) ---
f = forces(pos)
gr_accum = None
dr = 0.05
rmax = a / 2.0                                       # radii up to a/2 = 3.125
edges = np.arange(0.0, rmax + dr, dr)
centers = 0.5 * (edges[:-1] + edges[1:])
frame_count = 0

for step in range(nsteps):
    vel += 0.5 * dt * f
    pos += dt * vel
    pos %= a                                          # periodic wrap
    f = forces(pos)
    vel += 0.5 * dt * f

    # thermostat: rescale velocities toward T_target during equilibration
    if step % 100 == 0:
        ke = 0.5 * np.sum(vel**2)
        T_now = ke / N                                # 2D: KE = N*T (kB=1, 2 dof)
        if T_now > 0:
            vel *= np.sqrt(T_target / T_now)

    # collect trajectory frames near t = 100 (last 20 time units) for g(r)
    if step * dt >= t_final - 20.0:
        # -------- explicit g(r) measurement for this frame --------
        d = pos[:, None, :] - pos[None, :, :]
        d -= a * np.round(d / a)                      # minimum-image pair vectors
        rij = np.sqrt(np.sum(d**2, axis=2))
        iu = np.triu_indices(N, k=1)                  # unique pairs i<j
        dists = rij[iu]
        dists = dists[dists < rmax]                   # keep distances up to a/2
        counts, _ = np.histogram(dists, bins=edges)   # dN in each shell
        if gr_accum is None:
            gr_accum = np.zeros_like(centers)
        # each pair counted once -> multiply by 2 so every particle is a reference
        gr_accum += 2.0 * counts
        frame_count += 1

# --- normalise: <dN(r)> / (2*pi*r*dr*rho0), averaged over N references & frames ---
avg_dN = gr_accum / (frame_count * N)                 # <dN(r)> per reference particle
shell_area = 2.0 * np.pi * centers * dr               # shell area 2*pi*r*dr (2D)
g = avg_dN / (shell_area * rho0)

# --- locate the first and second peaks of g(r) ---
first_region = (centers > 0.8) & (centers < 1.5)
i1 = np.where(first_region)[0][np.argmax(g[first_region])]
second_region = (centers > 1.5) & (centers < 2.5)
i2 = np.where(second_region)[0][np.argmax(g[second_region])]

# --- classic-liquid-shape check ---
below = g[centers < 0.8].max()                        # should be ~0
tail = g[centers > 2.7].mean()                        # should settle toward 1
check = (below < 0.5) and (0.9 < centers[i1] < 1.3) and \
        (g[i1] > 1.5) and (1.7 < centers[i2] < 2.3) and \
        (g[i1] > g[i2]) and (0.7 < tail < 1.4)

print(f"mean density rho0 = N/a^2 = {rho0:.6f}")
print(f"number of trajectory frames averaged = {frame_count}")
print(f"first peak position  r1 = {centers[i1]:.4f}")
print(f"first peak height  g(r1) = {g[i1]:.4f}")
print(f"second peak position r2 = {centers[i2]:.4f}")
print(f"second peak height g(r2) = {g[i2]:.4f}")
print(f"max g(r) below r=0.8 = {below:.4f}")
print(f"mean g(r) at large r (r>2.7) = {tail:.4f}")
print(f"classic liquid shape check passed = {bool(check)}")

# --- plot g(r) vs r ---
plt.figure(figsize=(8, 5))
plt.plot(centers, g, "-", color="steelblue", lw=1.5, label="g(r)")
plt.axhline(1.0, color="gray", ls="--", lw=0.8, label="g=1 (ideal gas)")
plt.plot(centers[i1], g[i1], "o", color="red",
         label=f"1st peak r={centers[i1]:.2f}, g={g[i1]:.2f}")
plt.plot(centers[i2], g[i2], "s", color="green",
         label=f"2nd peak r={centers[i2]:.2f}, g={g[i2]:.2f}")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial distribution function of the simulated particle box")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/5C.7.1_s3.png")

# The check confirms the result because reproducing the textbook liquid signature
# (g(r) ~ 0 at short range, a sharp first peak near the LJ minimum r~1, a weaker
# second peak near r~2, and g(r) -> 1 at large r) is exactly what a correctly
# measured radial distribution function of a Lennard-Jones liquid must look like.
