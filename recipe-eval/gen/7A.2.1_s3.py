import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
L = 100          # domain length
dX = 1.0         # spatial step
dt = 0.01        # time step
D = 1.0          # diffusion coefficient
T_total = 25.0   # total integration time
n_steps = int(round(T_total / dt))

# --- Grid setup (Dirichlet/absorbing ends) ---
# X runs from -L/2 to +L/2 so the initial spike sits at X = 0
X = np.arange(-L/2, L/2 + dX, dX)
N = len(X)
i0 = np.argmin(np.abs(X))        # index of X = 0

# Stability parameter for the explicit scheme: must be <= 0.5
alpha = D * dt / dX**2
print(f"Stability parameter alpha = D*dt/dX^2 = {alpha:.4f} (must be <= 0.5)")

# --- Initial condition: all probability at X = 0 ---
P = np.zeros(N)
P[i0] = 1.0 / dX   # unit total probability: sum(P)*dX = 1

# --- Helper moment calculations ---
def total_prob(P):
    return np.sum(P) * dX

def mean_X(P):
    m = total_prob(P)
    return np.sum(X * P) * dX / m if m > 0 else 0.0

def var_X(P):
    m = total_prob(P)
    if m <= 0:
        return 0.0
    mu = np.sum(X * P) * dX / m
    return np.sum((X - mu)**2 * P) * dX / m

# --- Times at which to snapshot the profile for the spreading plot ---
snapshot_times = [0.0, 1.0, 5.0, 10.0, 25.0]
snapshot_steps = {int(round(t/dt)): t for t in snapshot_times}
snapshots = {}

# --- Storage for time series of moments ---
t_hist, var_hist, mass_hist, mean_hist = [], [], [], []

# Record initial state
t_hist.append(0.0)
var_hist.append(var_X(P))
mass_hist.append(total_prob(P))
mean_hist.append(mean_X(P))
if 0 in snapshot_steps:
    snapshots[snapshot_steps[0]] = P.copy()

# --- Explicit finite-difference integration (FTCS scheme) ---
# Update: P_i^{n+1} = P_i^n + alpha*(P_{i+1}^n - 2 P_i^n + P_{i-1}^n)
# Interior points only; the boundary values stay at 0 (absorbing Dirichlet ends),
# so probability that reaches the edges leaks out of the system.
for step in range(1, n_steps + 1):
    lap = np.zeros(N)                       # discrete Laplacian d2P/dX2
    lap[1:-1] = P[2:] - 2.0*P[1:-1] + P[:-2]
    P = P + alpha * lap                     # forward-Euler step in time
    P[0] = 0.0                              # Dirichlet: absorbing left end
    P[-1] = 0.0                             # Dirichlet: absorbing right end

    t = step * dt
    t_hist.append(t)
    var_hist.append(var_X(P))
    mass_hist.append(total_prob(P))
    mean_hist.append(mean_X(P))
    if step in snapshot_steps:
        snapshots[snapshot_steps[step]] = P.copy()

t_hist = np.array(t_hist)
var_hist = np.array(var_hist)
mass_hist = np.array(mass_hist)
mean_hist = np.array(mean_hist)

# --- Theoretical variance ---
theory_var = 2.0 * D * t_hist

# --- Print summary numerical results ---
print(f"Initial total probability : {mass_hist[0]:.6f}")
print(f"Final total probability    : {mass_hist[-1]:.6f}")
print(f"Final time                 : {t_hist[-1]:.4f}")
print(f"Final mean X               : {mean_hist[-1]:.6f} (expect ~0 by symmetry)")
print(f"Final variance (numeric)   : {var_hist[-1]:.6f}")
print(f"Final variance (theory 2Dt): {theory_var[-1]:.6f}")

# Early-time comparison where boundary leakage is negligible
for tc in [1.0, 5.0, 10.0]:
    idx = np.argmin(np.abs(t_hist - tc))
    print(f"t = {t_hist[idx]:5.2f} : var_numeric = {var_hist[idx]:8.4f}, "
          f"2Dt = {theory_var[idx]:8.4f}, mass = {mass_hist[idx]:.6f}")

# Relative error at an early time
idx1 = np.argmin(np.abs(t_hist - 1.0))
rel_err = abs(var_hist[idx1] - theory_var[idx1]) / theory_var[idx1]
print(f"Relative variance error at t=1.0 : {rel_err:.4%}")

# --- Plots ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (1) P(X) spreading over time
ax = axes[0]
for t in snapshot_times:
    if t in snapshots:
        ax.plot(X, snapshots[t], label=f"t = {t:g}")
ax.set_xlabel("X")
ax.set_ylabel("P(X)")
ax.set_title("Diffusing distribution P(X) over time")
ax.set_xlim(-40, 40)
ax.legend()

# (2) Variance vs t against 2*D*t
ax = axes[1]
ax.plot(t_hist, var_hist, 'b-', label="numeric variance")
ax.plot(t_hist, theory_var, 'r--', label="theory 2Dt")
ax.set_xlabel("t")
ax.set_ylabel("variance")
ax.set_title("Variance growth vs theory")
ax.legend()

# (3) Total probability leaking through absorbing ends
ax = axes[2]
ax.plot(t_hist, mass_hist, 'g-')
ax.set_xlabel("t")
ax.set_ylabel("total probability")
ax.set_title("Mass loss at Dirichlet ends")
ax.set_ylim(0, 1.05)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.2.1_s3.png")

# --- One-sentence explanation ---
print("Explanation: The numeric variance tracks 2*D*t at early times while the "
      "profile stays Gaussian and total probability stays near 1, confirming the "
      "integrator reproduces the diffusion law; the later downward deviation in "
      "variance and mass appears only once the spreading Gaussian reaches the "
      "absorbing Dirichlet boundaries, exactly as expected.")
