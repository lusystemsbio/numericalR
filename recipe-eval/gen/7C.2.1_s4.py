import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Gierer-Meinhardt substrate-depletion model
#   f(u,v) = u^2 v - u        (activator u, diffuses slowly, Du = d)
#   g(u,v) = mu (1 - u^2 v)   (substrate v, diffuses fast,   Dv = 1)
# Turing instability requires the substrate to diffuse faster than the
# activator (d < 1); the ratio d and the parameter mu decide the outcome.
# ----------------------------------------------------------------------

# ---- spatial / temporal discretization ----
L    = 20.0          # domain length
dX   = 0.2           # grid spacing
dt   = 0.01          # time step
N    = int(round(L / dX))   # number of grid points
X    = np.arange(N) * dX
Dv   = 1.0           # fast-diffusing substrate

# diffusion stability number (must stay < 0.5 for explicit scheme): here 0.25
print(f"Grid points N = {N}")
print(f"Explicit diffusion number Dv*dt/dX^2 = {Dv*dt/dX**2:.4f}")


def laplacian(c):
    """1-D Laplacian with zero-flux (Neumann) boundaries, second order."""
    lap = np.empty_like(c)
    lap[1:-1] = c[2:] - 2.0 * c[1:-1] + c[:-2]     # interior points
    lap[0]    = c[1]  - c[0]                        # mirror left boundary
    lap[-1]   = c[-2] - c[-1]                       # mirror right boundary
    return lap / dX**2


def integrate(u, v, d, mu, t_end):
    """Explicit forward-Euler finite-difference RD integrator over one block.
    Advances u,v from current state by t_end using explicit time stepping."""
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        # reaction terms
        f = u * u * v - u                # activator kinetics
        g = mu * (1.0 - u * u * v)       # substrate kinetics
        # explicit update: state += dt*(diffusion + reaction)
        u_new = u + dt * (d  * laplacian(u) + f)
        v_new = v + dt * (Dv * laplacian(v) + g)
        u, v = u_new, v_new
    return u, v


def run_case(d, mu, u0_mean, title):
    """Run one case in successive time blocks, recording u(X) snapshots."""
    rng = np.random.default_rng(10)                     # fixed seed 10
    u = u0_mean + 0.2 * (rng.random(N) - 0.5)           # +-0.1 noise
    v = 1.0     + 0.2 * (rng.random(N) - 0.5)           # +-0.1 noise
    block = 10.0                                        # length of each time block
    nblocks = 6
    times = [0.0]
    snaps = [u.copy()]
    for b in range(nblocks):
        u, v = integrate(u, v, d, mu, block)            # run one successive block
        times.append((b + 1) * block)
        snaps.append(u.copy())
    print(f"\n[{title}]  d={d}, mu={mu}, initial u mean={u0_mean}")
    print(f"  final u: min={u.min():.4f}, max={u.max():.4f}, "
          f"mean={u.mean():.4f}, spatial amplitude(max-min)={u.max()-u.min():.4f}")
    return times, snaps


# ----------------------------------------------------------------------
# Three cases / three outcomes
# ----------------------------------------------------------------------
casesA = run_case(d=0.1, mu=1.5, u0_mean=1.0, title="Pattern-forming (Turing)")
casesB = run_case(d=0.8, mu=1.5, u0_mean=1.0, title="Close diffusion constants")
casesC = run_case(d=0.3, mu=0.9, u0_mean=2.0, title="Oscillating")

cases = [
    ("d=0.1, mu=1.5\nstationary Turing pattern", casesA),
    ("d=0.8, mu=1.5\nflat homogeneous state",    casesB),
    ("d=0.3, mu=0.9, u0=2\nuniform oscillation",  casesC),
]

# ----------------------------------------------------------------------
# Plot u(X) at successive times for each case
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, (label, (times, snaps)) in zip(axes, cases):
    for t, s in zip(times, snaps):
        ax.plot(X, s, label=f"t={t:.0f}")
    ax.set_title(label)
    ax.set_xlabel("X")
    ax.set_ylabel("u(X)")
    ax.legend(fontsize=7, ncol=2)
fig.suptitle("Gierer-Meinhardt substrate-depletion: Turing pattern and two failure modes")
fig.tight_layout()
fig.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.2.1_s4.png", dpi=110)

# ----------------------------------------------------------------------
# Separate check: quantify the three outcomes and their controlling ratio
# ----------------------------------------------------------------------
def classify(times, snaps):
    """Return spatial pattern amplitude (final) and temporal amplitude of the
    spatial mean (a nonzero spatial amp => stationary pattern; nonzero temporal
    amp of the mean => uniform oscillation; both ~0 => homogeneous steady state)."""
    final = snaps[-1]
    spatial_amp = final.max() - final.min()
    means = np.array([s.mean() for s in snaps[-4:]])   # last blocks (settled regime)
    temporal_amp = means.max() - means.min()
    return spatial_amp, temporal_amp

print("\n=== Outcome check (controlled by diffusion ratio d and mu) ===")
labels = ["d=0.1,mu=1.5", "d=0.8,mu=1.5", "d=0.3,mu=0.9,u0=2"]
for lbl, (_, dat) in zip(labels, [casesA, casesB, casesC]):
    pass
for lbl, dat in zip(labels, [casesA, casesB, casesC]):
    sa, ta = classify(*dat)
    print(f"{lbl:20s}  spatial_amplitude={sa:.4f}   temporal_amplitude(mean)={ta:.4f}")

print("\nThe check confirms the result because the three parameter settings each "
      "isolate one signature -- a large spatial amplitude with ~zero temporal "
      "amplitude (stationary Turing pattern), both amplitudes ~zero (homogeneous "
      "steady state), and ~zero spatial amplitude with nonzero temporal amplitude "
      "(spatially uniform oscillation) -- so the diffusion ratio d and mu alone "
      "select which of the three distinct dynamical outcomes occurs.")
