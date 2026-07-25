import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Regenerate the t = 100 trajectory from 5C.5: a 2D Lennard-Jones fluid of
# N = 25 particles in a periodic box of side a = 6.25.  We run velocity-Verlet
# molecular dynamics to t = 100 and accumulate g(r) over the trajectory frames.
# ---------------------------------------------------------------------------
rng = np.random.default_rng(12345)

N = 25                      # number of particles
a = 6.25                    # box side length
rho0 = N / a**2             # mean number density  rho0 = N/a^2
T_target = 1.0              # reduced temperature (liquid state)
dt = 0.004                  # MD time step
rc = 2.5                    # LJ interaction cutoff
rc2 = rc * rc

# start particles on a regular grid (slightly jittered) with random velocities
side = int(np.ceil(np.sqrt(N)))
xs = np.linspace(0, a, side, endpoint=False)
grid = np.array([[x, y] for y in xs for x in xs])[:N]
pos = (grid + 0.01 * rng.standard_normal((N, 2))) % a
vel = rng.standard_normal((N, 2))
vel -= vel.mean(axis=0)             # remove net momentum


def forces(pos):
    """Minimum-image Lennard-Jones forces with periodic boundaries."""
    f = np.zeros_like(pos)
    for i in range(N - 1):
        d = pos[i] - pos[i + 1:]
        d -= a * np.round(d / a)                    # minimum image convention
        r2 = np.sum(d * d, axis=1)
        m = r2 < rc2                                # apply cutoff
        inv2 = 1.0 / r2[m]
        inv6 = inv2**3
        # LJ force / r :  24*(2 r^-12 - r^-6) / r^2
        fmag = (48.0 * inv6 * inv6 - 24.0 * inv6) * inv2
        fpair = fmag[:, None] * d[m]
        f[i] += fpair.sum(axis=0)                    # Newton's third law
        idx = np.arange(i + 1, N)[m]
        np.add.at(f, idx, -fpair)
    return f


# ---------------------------------------------------------------------------
# g(r) histogram set-up:  bins of width dr out to r = a/2 = 3.125.
# ---------------------------------------------------------------------------
dr = 0.05
rmax = a / 2.0                          # 3.125
nbins = int(rmax / dr)
edges = np.arange(nbins + 1) * dr
rmid = 0.5 * (edges[:-1] + edges[1:])   # bin-center radii
hist = np.zeros(nbins)                  # accumulated shell counts
nframes = 0

# ---------------------------------------------------------------------------
# Velocity-Verlet integration up to t = 100.  After equilibration we sample
# minimum-image pair distances into the histogram every 20 steps so that
# <dN(r)> is averaged over many independent liquid configurations.
# ---------------------------------------------------------------------------
f = forces(pos)
nsteps = int(100 / dt)
equil = nsteps // 2                     # discard first half as equilibration
for step in range(nsteps):
    vel += 0.5 * dt * f                 # velocity Verlet, half kick
    pos = (pos + dt * vel) % a          # drift + periodic wrap
    f = forces(pos)
    vel += 0.5 * dt * f                 # second half kick
    # simple velocity-rescaling thermostat (2D: T = KE/N by equipartition)
    ke = 0.5 * np.sum(vel * vel)
    T_now = ke / N
    vel *= np.sqrt(T_target / T_now)
    # sample g(r) after equilibration
    if step >= equil and step % 20 == 0:
        for i in range(N - 1):
            d = pos[i] - pos[i + 1:]
            d -= a * np.round(d / a)                 # minimum-image distances
            rij = np.sqrt(np.sum(d * d, axis=1))
            # each pair contributes to BOTH reference particles -> count twice
            h, _ = np.histogram(np.concatenate([rij, rij]), bins=edges)
            hist += h
        nframes += 1

# ---------------------------------------------------------------------------
# Normalize: <dN(r)> per reference particle, divided by the shell area
# 2*pi*r*dr and the mean density rho0, i.e. g(r)=<dN(r)>/(2*pi*r*dr*rho0).
# ---------------------------------------------------------------------------
avg_dN = hist / (N * nframes)                # average count per reference particle
shell_area = 2.0 * np.pi * rmid * dr         # area of each annular shell
g = avg_dN / (shell_area * rho0)

# locate first peak (near r=1) and second peak (near r=2)
i1 = np.argmax(g * (rmid < 1.6))
sr = (rmid > 1.5) & (rmid < 2.6)
i2 = np.where(sr)[0][np.argmax(g[sr])]

# ---------------------------------------------------------------------------
# Report numerical results.
# ---------------------------------------------------------------------------
print(f"Mean density rho0 = N/a^2 = {rho0:.5f}")
print(f"Number of trajectory frames averaged = {nframes}")
print(f"First  peak: r = {rmid[i1]:.3f}, g(r) = {g[i1]:.3f}")
print(f"Second peak: r = {rmid[i2]:.3f}, g(r) = {g[i2]:.3f}")
print(f"g(r) in first bin (r = {rmid[0]:.3f}) = {g[0]:.3f}")
print(f"Mean g(r) for r < 0.8 (excluded-volume core) = {g[rmid < 0.8].mean():.3f}")
print(f"Mean g(r) for r > 2.6 (large-r tail, should approach 1) = {g[rmid > 2.6].mean():.3f}")

# ---------------------------------------------------------------------------
# Classic-liquid-shape check.
# ---------------------------------------------------------------------------
core_ok = g[rmid < 0.8].mean() < 0.3
first_ok = 0.8 < rmid[i1] < 1.4 and g[i1] > 1.5
second_ok = 1.7 < rmid[i2] < 2.4 and g[i2] > 1.0 and g[i2] < g[i1]
tail_ok = abs(g[rmid > 2.6].mean() - 1.0) < 0.3
print(f"Check - near-zero core (r<0.8):           {core_ok}")
print(f"Check - sharp first peak near r=1:        {first_ok}")
print(f"Check - weaker second peak near r=2:      {second_ok}")
print(f"Check - settles toward 1 at large r:      {tail_ok}")
print(f"Classic liquid g(r) confirmed: {core_ok and first_ok and second_ok and tail_ok}")
# One-sentence explanation: the check confirms the result because the
# excluded-volume gap, a dominant first peak at the LJ equilibrium separation,
# a weaker second shell, and decay to 1 are exactly the short-range order
# (and absence of long-range order) that defines a liquid, so reproducing that
# signature means our histogram-and-normalize g(r) captured the true structure.
print("Explanation: reproducing the near-zero core, a sharp first peak at the "
      "LJ equilibrium separation, a weaker second peak, and decay toward 1 "
      "matches the unique short-range-order signature of a liquid, so the "
      "computed g(r) faithfully represents the box's real structure.")

# ---------------------------------------------------------------------------
# Plot g(r) versus r with the first and second peaks marked.
# ---------------------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rmid, g, '-o', ms=3, color='C0', label='g(r)')
plt.axhline(1.0, color='gray', ls='--', lw=0.8, label='ideal gas (g=1)')
plt.plot(rmid[i1], g[i1], 'rs', ms=9, label=f'1st peak (r={rmid[i1]:.2f})')
plt.plot(rmid[i2], g[i2], 'g^', ms=9, label=f'2nd peak (r={rmid[i2]:.2f})')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Radial distribution function (2D LJ liquid, N=25, a=6.25)')
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/5C.7.1_s4.png")
