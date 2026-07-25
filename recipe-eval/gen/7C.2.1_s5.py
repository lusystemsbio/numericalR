import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Gierer-Meinhardt substrate-depletion model
#   f(u,v) = u^2 v - u        (activator u, slow diffusion Du = d)
#   g(u,v) = mu (1 - u^2 v)    (substrate v, fast diffusion Dv = 1)
# Turing patterns need the substrate to diffuse faster than u.
# ---------------------------------------------------------------

# --- Reaction terms ---
def f(u, v):
    return u**2 * v - u

def g(u, v, mu):
    return mu * (1.0 - u**2 * v)

# --- Explicit 1D Laplacian with no-flux (Neumann) boundaries ---
def laplacian(c, dX):
    lap = np.empty_like(c)
    lap[1:-1] = (c[2:] - 2.0*c[1:-1] + c[:-2]) / dX**2   # interior
    lap[0]    = (c[1]  - c[0]) / dX**2 * 2.0             # left no-flux (ghost = c[1])
    lap[-1]   = (c[-2] - c[-1]) / dX**2 * 2.0            # right no-flux (ghost = c[-2])
    return lap

# --- Multi-component finite-difference RD integrator (from 7C.1) ---
# Explicit forward-Euler step for the two coupled fields.
def rd_step(u, v, d, mu, dX, dt):
    u_new = u + dt * (d   * laplacian(u, dX) + f(u, v))      # slow activator
    v_new = v + dt * (1.0 * laplacian(v, dX) + g(u, v, mu))  # fast substrate
    return u_new, v_new

# --- Run the integrator in successive time blocks, recording snapshots ---
def run_case(u0, v0, d, mu, dX, dt, snap_times):
    u, v = u0.copy(), v0.copy()
    snap_times = sorted(snap_times)
    snaps = {}
    steps_per_time = {t: int(round(t/dt)) for t in snap_times}
    t = 0.0
    step = 0
    max_step = max(steps_per_time.values())
    # record initial state if requested
    if 0 in [int(round(x/dt)) for x in snap_times]:
        pass
    target_steps = sorted(set(steps_per_time.values()))
    ti = 0
    while step <= max_step:
        # record snapshots whose step index matches the current step
        for tval, s in steps_per_time.items():
            if s == step:
                snaps[tval] = u.copy()
        if step == max_step:
            break
        u, v = rd_step(u, v, d, mu, dX, dt)   # one explicit block-step
        step += 1
    return snaps, u, v

# --- Domain / discretization ---
L   = 20.0
dX  = 0.2
dt  = 0.01
N   = int(round(L/dX)) + 1
X   = np.linspace(0.0, L, N)

# Stability check for explicit scheme (Dv is the largest diffusivity)
diff_number = 1.0 * dt / dX**2
print(f"Explicit diffusion number (Dv*dt/dX^2): {diff_number:.4f}  (needs < 0.5)")

# ---------------------------------------------------------------
# Case 1: Turing pattern     d = 0.1, mu = 1.5, u=v=1 +- 0.1 noise
# Case 2: Flat / uniform     d = 0.8, mu = 1.5, u=v=1 +- 0.1 noise
# Case 3: Uniform oscillation d = 0.3, mu = 0.9, initial u = 2
# ---------------------------------------------------------------

# Initial conditions for cases 1 & 2 (near uniform steady state u=v=1)
np.random.seed(10)
u_init_12 = 1.0 + 0.1*(2*np.random.rand(N) - 1)
v_init_12 = 1.0 + 0.1*(2*np.random.rand(N) - 1)

# Initial condition for case 3 (u = 2)
np.random.seed(10)
u_init_3 = 2.0 + 0.1*(2*np.random.rand(N) - 1)
v_init_3 = 1.0 + 0.1*(2*np.random.rand(N) - 1)

cases = [
    ("Turing pattern (d=0.1, mu=1.5)", u_init_12, v_init_12, 0.1, 1.5,
     [0, 20, 50, 100, 200]),
    ("Close diffusion (d=0.8, mu=1.5)", u_init_12, v_init_12, 0.8, 1.5,
     [0, 20, 50, 100, 200]),
    ("Oscillation (d=0.3, mu=0.9, u0=2)", u_init_3, v_init_3, 0.3, 0.9,
     [0, 5, 10, 15, 20, 25]),
]

fig, axes = plt.subplots(3, 1, figsize=(9, 12))

results = {}
for ax, (title, u0, v0, d, mu, snaps_t) in zip(axes, cases):
    snaps, u_fin, v_fin = run_case(u0, v0, d, mu, dX, dt, snaps_t)
    results[title] = (snaps, u_fin)
    for tval in sorted(snaps):
        ax.plot(X, snaps[tval], label=f"t={tval}")
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("u")
    ax.legend(fontsize=8, ncol=3)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.2.1_s5.png")

# ---------------------------------------------------------------
# Quantitative check: classify each outcome by
#   (a) spatial variation of the final u(X)  -> pattern vs flat
#   (b) temporal variation of the mean u     -> oscillation vs steady
# ---------------------------------------------------------------

# Measure temporal behaviour of the spatial mean of u over a final window.
def temporal_amplitude(u0, v0, d, mu, dX, dt, n_burn, n_watch):
    u, v = u0.copy(), v0.copy()
    for _ in range(n_burn):
        u, v = rd_step(u, v, d, mu, dX, dt)
    means = []
    for _ in range(n_watch):
        u, v = rd_step(u, v, d, mu, dX, dt)
        means.append(u.mean())
    means = np.array(means)
    return means.max() - means.min()

print("\n--- Outcome classification ---")
for title, u0, v0, d, mu, snaps_t in cases:
    snaps, u_fin = results[title]
    spatial_var = u_fin.max() - u_fin.min()          # spatial amplitude at final time
    temp_var = temporal_amplitude(u0, v0, d, mu, dX, dt,
                                  n_burn=3000, n_watch=500)  # temporal amplitude of mean(u)
    print(f"\n{title}")
    print(f"  final spatial amplitude (max-min of u(X)): {spatial_var:.4e}")
    print(f"  temporal amplitude of mean(u):             {temp_var:.4e}")
    if temp_var > 1e-2:
        outcome = "spatially UNIFORM TEMPORAL OSCILLATION"
    elif spatial_var > 1e-2:
        outcome = "stationary periodic TURING PATTERN"
    else:
        outcome = "HOMOGENEOUS STEADY STATE"
    print(f"  => {outcome}")

# One-sentence explanation of why this check confirms the result:
print("\nExplanation:")
print("This check confirms the result because measuring the spatial amplitude "
      "of the final u(X) and the temporal amplitude of its mean cleanly "
      "separates the three regimes: only the pattern case retains large spatial "
      "variation with no temporal change, the close-diffusion case decays to a "
      "flat homogeneous steady state, and the oscillating case shows a large "
      "time-varying but spatially uniform mean.")
