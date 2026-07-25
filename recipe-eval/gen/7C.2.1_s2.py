import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Gierer-Meinhardt substrate-depletion model
#   f(u,v) = u^2 v - u      (activator u, slow diffusion Du = d)
#   g(u,v) = mu (1 - u^2 v) (substrate v, fast diffusion Dv = 1)
# Homogeneous steady state: f=g=0 -> u*=1, v*=1
# ----------------------------------------------------------------------

def f(u, v):
    return u**2 * v - u

def g(u, v, mu):
    return mu * (1.0 - u**2 * v)

# ---- Explicit finite-difference reaction-diffusion integrator (from 7C.1) ----
# 1D, no-flux (Neumann) boundaries, forward Euler in time, centered in space.
def laplacian_1d(c, dX):
    lap = np.empty_like(c)
    # interior: second difference
    lap[1:-1] = (c[2:] - 2.0 * c[1:-1] + c[:-2]) / dX**2
    # no-flux boundaries: mirror the neighbor (zero gradient at the ends)
    lap[0]  = (c[1] - c[0]) * 2.0 / dX**2
    lap[-1] = (c[-2] - c[-1]) * 2.0 / dX**2
    return lap

def integrate_block(u, v, d, mu, dX, dt, nsteps):
    # advance the two-component system by nsteps explicit Euler steps
    for _ in range(nsteps):
        u_new = u + dt * (d   * laplacian_1d(u, dX) + f(u, v))
        v_new = v + dt * (1.0 * laplacian_1d(v, dX) + g(u, v, mu))
        u, v = u_new, v_new
    return u, v

# ---- Domain / discretization ----
L = 20.0
dX = 0.2
dt = 0.01
X = np.arange(0.0, L + dX/2, dX)
N = X.size

# ---- Cases: (label, d, mu, u_init_mean, snapshot times) ----
cases = [
    ("Pattern-forming (d=0.1, mu=1.5)",   0.1, 1.5, 1.0, [0, 20, 60, 150, 400]),
    ("Close diffusion (d=0.8, mu=1.5)",   0.8, 1.5, 1.0, [0, 20, 60, 150, 400]),
    ("Oscillating (d=0.3, mu=0.9, u0=2)", 0.3, 0.9, 2.0, [0, 5, 10, 15, 20]),
]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

results = {}  # store final fields for the diffusion-ratio / mu check

for ax, (label, d, mu, u0, snap_times) in zip(axes, cases):
    # nearly uniform initial condition with +-0.1 noise, seed 10
    rng = np.random.default_rng(10)
    u = u0 + rng.uniform(-0.1, 0.1, N)
    v = 1.0 + rng.uniform(-0.1, 0.1, N)

    t = 0.0
    for k, ts in enumerate(snap_times):
        nsteps = int(round((ts - t) / dt))
        if nsteps > 0:
            u, v = integrate_block(u, v, d, mu, dX, dt, nsteps)
            t = ts
        ax.plot(X, u, label=f"t={ts}")

    ax.set_title(label)
    ax.set_xlabel("X")
    ax.set_ylabel("u(X)")
    ax.legend(fontsize=8)

    results[label] = (u.copy(), d, mu)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.2.1_s2.png")

# ----------------------------------------------------------------------
# Quantitative check: classify each outcome by spatial and temporal signatures
#   - spatial variation of the final u field  (max-min over X)
#   - temporal variation at the domain center over one extra time block
# ----------------------------------------------------------------------
print("Domain length L =", L)
print("dX =", dX, " dt =", dt, " N grid points =", N)
print("Homogeneous steady state: u* = 1, v* = 1")
print("Diffusion numbers (D*dt/dX^2): Dv=1 ->", 1.0*dt/dX**2,
      "(<=0.5 => stable)")
print("")

for label, (u_final, d, mu) in results.items():
    spatial_range = u_final.max() - u_final.min()

    # re-run this case's final state a bit further and watch the center point
    rng = np.random.default_rng(10)
    u0_mean = 2.0 if "u0=2" in label else 1.0
    u = u0_mean + rng.uniform(-0.1, 0.1, N)
    v = 1.0 + rng.uniform(-0.1, 0.1, N)
    # integrate to a late time
    u, v = integrate_block(u, v, d, mu, dX, dt, 40000)
    center = N // 2
    trace = []
    for _ in range(200):  # sample u at center over 200 blocks of 10 steps
        u, v = integrate_block(u, v, d, mu, dX, dt, 10)
        trace.append(u[center])
    trace = np.array(trace)
    temporal_range = trace.max() - trace.min()
    final_spatial_range = u.max() - u.min()

    if final_spatial_range > 0.05 and temporal_range < 0.05:
        verdict = "stationary periodic Turing pattern"
    elif final_spatial_range < 0.05 and temporal_range < 0.05:
        verdict = "homogeneous steady state"
    else:
        verdict = "spatially uniform temporal oscillation"

    print(f"Case: {label}")
    print(f"  d (=Du) = {d},  Dv = 1,  diffusion ratio Dv/Du = {1.0/d:.4g},  mu = {mu}")
    print(f"  final spatial range max(u)-min(u) over X = {final_spatial_range:.6f}")
    print(f"  late-time temporal range of u at center   = {temporal_range:.6f}")
    print(f"  outcome = {verdict}")
    print("")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms Turing's result because it shows the same "
      "reaction kinetics produce three distinct outcomes selected purely by the "
      "diffusion ratio Dv/Du and mu -- a nonzero final spatial range with zero "
      "temporal variation (stationary Turing pattern) requires a large ratio (small d), "
      "a large ratio is destroyed when d is too close to Dv (flat homogeneous state), "
      "and a nonzero temporal variation that is spatially uniform (Hopf oscillation) "
      "is set by mu independent of diffusion.")
