import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Reproduce a "5C.5"-style trajectory: 25 Lennard-Jones particles in a 2D box
# of side a = 6.25, integrated with velocity-Verlet out to t = 100.  We then
# measure the radial distribution function g(r) from the trajectory frames.
# -----------------------------------------------------------------------------
np.random.seed(5)

N = 25                       # number of particles
a = 6.25                     # box side length
rho0 = N / a**2              # mean (2D) number density: rho0 = N / a^2
dt = 0.005                   # integration time step
t_final = 100.0              # final trajectory time (from 5C.5)
nsteps = int(round(t_final / dt))
rc = a / 2.0                 # potential cutoff = a/2 (also max g(r) radius)

# --- initial condition: place particles on a 5x5 grid, give random velocities
grid = np.linspace(0.0, a, int(np.sqrt(N)), endpoint=False)
xs, ys = np.meshgrid(grid, grid)
pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
pos += 0.01 * (np.random.rand(N, 2) - 0.5)          # tiny jitter off the lattice
vel = np.random.randn(N, 2)
vel -= vel.mean(axis=0)                              # remove net momentum

def forces(pos):
    """Lennard-Jones forces with minimum-image convention; returns forces + PE."""
    F = np.zeros_like(pos)
    U = 0.0
    for i in range(N - 1):
        d = pos[i + 1:] - pos[i]                     # displacements to j>i
        d -= a * np.round(d / a)                     # minimum image
        r2 = np.sum(d * d, axis=1)
        mask = r2 < rc * rc                          # apply cutoff
        r2 = r2[mask]
        dm = d[mask]
        inv2 = 1.0 / r2
        inv6 = inv2**3
        inv12 = inv6**2
        # LJ force magnitude/r: 24*(2/r^12 - 1/r^6)/r^2 ; times displacement vector
        fmag = (24.0 * (2.0 * inv12 - inv6) * inv2)[:, None]
        fij = fmag * dm
        F[i] += fij.sum(axis=0)
        idx = np.where(mask)[0] + (i + 1)
        np.add.at(F, idx, -fij)
        U += np.sum(4.0 * (inv12 - inv6))
    return F, U

# --- histogram setup for g(r)
dr = 0.05                                            # bin width
nbins = int(rc / dr)
edges = np.linspace(0.0, nbins * dr, nbins + 1)
r_centers = 0.5 * (edges[:-1] + edges[1:])
hist = np.zeros(nbins)                               # accumulated pair counts
nframes = 0

def pair_distances(pos):
    """All unique minimum-image pair distances in the box."""
    dists = []
    for i in range(N - 1):
        d = pos[i + 1:] - pos[i]
        d -= a * np.round(d / a)                     # minimum image
        dists.append(np.sqrt(np.sum(d * d, axis=1)))
    return np.concatenate(dists)

# --- velocity-Verlet integration; thermostat (rescale) during equilibration
F, U = forces(pos)
T_target = 1.0
equil = nsteps // 2                                  # first half = equilibration
for step in range(nsteps):
    vel += 0.5 * dt * F                              # half kick
    pos += dt * vel                                  # drift
    pos %= a                                         # wrap into box
    F, U = forces(pos)
    vel += 0.5 * dt * F                              # half kick
    if step < equil and step % 50 == 0:             # crude velocity rescaling
        ke = 0.5 * np.sum(vel * vel)
        T_now = ke / N                               # 2D: <KE> = N*k_B*T (k_B=1)
        if T_now > 0:
            vel *= np.sqrt(T_target / T_now)
    if step >= equil and step % 10 == 0:            # sample after equilibration
        d = pair_distances(pos)
        h, _ = np.histogram(d, bins=edges)
        hist += h
        nframes += 1

# --- normalize: <dN(r)> = 2*pairs / (N*nframes); each pair seen by both members.
avg_dN = 2.0 * hist / (N * nframes)                  # avg # neighbors per particle per shell
shell_area = 2.0 * np.pi * r_centers * dr            # 2D annulus area
g = avg_dN / (shell_area * rho0)                     # g(r) = <dN> / (2*pi*r*dr*rho0)

# --- locate first and second peaks
valid = r_centers > 0.3                              # ignore the empty core
gi = np.where(valid, g, 0.0)
i1 = np.argmax(gi)
r1, g1 = r_centers[i1], g[i1]
# second peak: search beyond the first minimum after the first peak
after = np.arange(i1 + 1, nbins)
min_after = after[np.argmin(g[after])] if len(after) else i1
region = np.arange(min_after, nbins)
i2 = region[np.argmax(g[region])] if len(region) else i1
r2, g2 = r_centers[i2], g[i2]

# --- classic-liquid-shape checks
g_below = g[r_centers < 0.8].max() if np.any(r_centers < 0.8) else 0.0
g_tail = g[r_centers > 2.5].mean() if np.any(r_centers > 2.5) else np.nan

print(f"Number of particles N            = {N}")
print(f"Box side a                       = {a}")
print(f"Mean density rho0 = N/a^2        = {rho0:.6f}")
print(f"Bin width dr                     = {dr}")
print(f"Max radius a/2                   = {rc}")
print(f"Frames averaged                  = {nframes}")
print(f"First peak:  r1 = {r1:.3f}, g(r1) = {g1:.3f}")
print(f"Second peak: r2 = {r2:.3f}, g(r2) = {g2:.3f}")
print(f"Max g(r) for r < 0.8            = {g_below:.3f}")
print(f"Mean g(r) for r > 2.5 (tail)    = {g_tail:.3f}")
print(f"Check near-zero core (<0.05)?   = {g_below < 0.05}")
print(f"Check first peak near r=1?      = {0.9 <= r1 <= 1.2}")
print(f"Check second peak near r=2?     = {1.7 <= r2 <= 2.3}")
print(f"Check tail settles toward 1?    = {abs(g_tail - 1.0) < 0.3}")

# --- plot g(r) vs r with the two peaks marked
plt.figure(figsize=(8, 5))
plt.plot(r_centers, g, "-o", ms=3, lw=1.2, label="g(r)")
plt.axhline(1.0, color="gray", ls="--", lw=1, label="g = 1 (ideal)")
plt.plot(r1, g1, "r^", ms=10, label=f"1st peak (r={r1:.2f})")
plt.plot(r2, g2, "gs", ms=10, label=f"2nd peak (r={r2:.2f})")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial distribution function g(r) (25 LJ particles, a=6.25, t=100)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5C.7.1_s1.png")

# Explanation of the check:
# The classic-liquid shape confirms the result because g(r) is built purely from
# geometric pair-distance statistics, so recovering a near-zero excluded core,
# a sharp peak at the LJ equilibrium separation r~1, a weaker shell near r~2, and
# decay to 1 shows the code reproduces the known structural signature of a
# Lennard-Jones liquid rather than random noise or an ordered/gaseous artifact.
print("Check explanation: the near-zero core, sharp first peak at the LJ "
      "equilibrium separation (~1), weaker second shell (~2), and decay to 1 "
      "are the known structural signature of an LJ liquid, so matching them "
      "confirms g(r) is computed correctly.")
