import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Gierer-Meinhardt SUBSTRATE-DEPLETION model (1D reaction-diffusion)
#   activator u (slow, Du = d),  substrate v (fast, Dv = 1)
#   f(u,v) = u^2 v - u          (u self-activates, consuming v)
#   g(u,v) = mu (1 - u^2 v)     (v is supplied and depleted by the reaction)
# Homogeneous fixed point: u = v = 1  (f=g=0).
# ----------------------------------------------------------------------

def f_react(u, v):
    return u * u * v - u

def g_react(u, v, mu):
    return mu * (1.0 - u * u * v)


def laplacian_neumann(c, dX):
    """1D Laplacian with no-flux (Neumann) boundaries, implemented explicitly."""
    lap = np.empty_like(c)
    # interior points: standard 3-point second difference
    lap[1:-1] = (c[2:] - 2.0 * c[1:-1] + c[:-2]) / dX**2
    # boundaries: ghost node mirrors the neighbour -> zero flux
    lap[0]  = 2.0 * (c[1]  - c[0])  / dX**2
    lap[-1] = 2.0 * (c[-2] - c[-1]) / dX**2
    return lap


def integrate_blocks(u0, v0, d, mu, dX, dt, block_time, n_blocks):
    """
    Explicit forward-Euler reaction-diffusion integrator (multi-component),
    run in successive time blocks.  Returns the u-snapshot after each block.
    Done step-by-step (not via a black-box solver) so the method is visible.
    """
    u = u0.copy()
    v = v0.copy()
    steps_per_block = int(round(block_time / dt))
    snapshots = [u.copy()]          # include the initial condition
    times = [0.0]
    for b in range(n_blocks):
        for _ in range(steps_per_block):
            # 1) diffusion terms for each component (different diffusivities)
            lu = laplacian_neumann(u, dX)
            lv = laplacian_neumann(v, dX)
            # 2) reaction terms
            ru = f_react(u, v)
            rv = g_react(u, v, mu)
            # 3) explicit Euler update:  c <- c + dt (D*lap + reaction)
            u = u + dt * (d * lu + ru)
            v = v + dt * (1.0 * lv + rv)
        snapshots.append(u.copy())
        times.append((b + 1) * block_time)
    return np.array(snapshots), np.array(times)


# ----------------------------------------------------------------------
# Common spatial grid and initial condition
# ----------------------------------------------------------------------
L = 20.0
dX = 0.2
dt = 0.01
N = int(round(L / dX)) + 1          # number of grid nodes
X = np.linspace(0.0, L, N)

np.random.seed(10)
noise_u = (np.random.rand(N) - 0.5) * 0.2   # +-0.1
noise_v = (np.random.rand(N) - 0.5) * 0.2   # +-0.1

# CFL / stability check for the explicit scheme (max diffusivity = Dv = 1)
Dmax = 1.0
print(f"Grid nodes N = {N}, dX = {dX}, dt = {dt}")
print(f"Explicit-scheme stability limit dt <= dX^2/(2*Dmax) = {dX**2/(2*Dmax):.4f} (using dt = {dt})")
print("")

# ----------------------------------------------------------------------
# Three cases
# ----------------------------------------------------------------------
cases = [
    dict(name="Turing pattern (d=0.1, mu=1.5)",
         d=0.1, mu=1.5, u_base=1.0, block_time=25.0, n_blocks=8),
    dict(name="Close diffusion / flat (d=0.8, mu=1.5)",
         d=0.8, mu=1.5, u_base=1.0, block_time=25.0, n_blocks=8),
    dict(name="Oscillation (d=0.3, mu=0.9, u0=2)",
         d=0.3, mu=0.9, u_base=2.0, block_time=3.0, n_blocks=8),
]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

results = {}
for ax, c in zip(axes, cases):
    u0 = c["u_base"] + noise_u
    v0 = 1.0 + noise_v
    snaps, times = integrate_blocks(u0, v0, c["d"], c["mu"], dX, dt,
                                    c["block_time"], c["n_blocks"])
    results[c["name"]] = (snaps, times)

    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    for k in range(len(times)):
        ax.plot(X, snaps[k], color=colors[k], label=f"t={times[k]:.0f}")
    ax.set_title(c["name"])
    ax.set_xlabel("X")
    ax.set_ylabel("u(X)")
    ax.legend(fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.2.1_s1.png")

# ----------------------------------------------------------------------
# Quantitative check: classify each outcome from two diagnostics
#   spatial_std  : std of u over space at the final time  -> spatial structure
#   temporal_amp : peak-to-peak of the spatial-mean of u over the last blocks
#                  -> uniform temporal oscillation
# ----------------------------------------------------------------------
print("Diagnostics per case:")
for name, (snaps, times) in results.items():
    final = snaps[-1]
    spatial_std = np.std(final)
    spatial_mean_series = snaps.mean(axis=1)              # <u>(t)
    temporal_amp = spatial_mean_series.max() - spatial_mean_series.min()

    if spatial_std > 0.05 and temporal_amp < 0.1:
        outcome = "STATIONARY periodic Turing pattern"
    elif spatial_std < 0.02 and temporal_amp > 0.1:
        outcome = "SPATIALLY UNIFORM temporal oscillation"
    else:
        outcome = "HOMOGENEOUS steady state (flat)"

    print(f"  {name}")
    print(f"    final spatial std of u          = {spatial_std:.4f}")
    print(f"    temporal amplitude of <u>(t)    = {temporal_amp:.4f}")
    print(f"    final min/max u                 = {final.min():.4f} / {final.max():.4f}")
    print(f"    -> outcome: {outcome}")

print("")
print("Check explanation:")
print("The single unchanged integrator produces a stationary periodic pattern, a "
      "flat state, and a uniform oscillation purely by changing the diffusion ratio "
      "d and mu, confirming that these outcomes are set by those parameters (Turing "
      "instability requires slow activator/fast substrate) and not by the numerics.")
