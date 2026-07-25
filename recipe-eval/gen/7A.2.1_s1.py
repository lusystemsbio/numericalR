import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 1D explicit finite-difference diffusion integrator (from 7A.1)
#   dP/dt = D * d2P/dX2
# Discretized explicitly:
#   P_i^{n+1} = P_i^n + (D*dt/dX^2) * (P_{i+1} - 2*P_i + P_{i-1})
# Dirichlet (absorbing) boundaries: ends held at 0, so probability
# that reaches the edges is removed from the system.
# ---------------------------------------------------------------

# --- Parameters ---
L = 100          # domain length
dX = 1.0         # spatial step
dt = 0.01        # time step
D = 1.0          # diffusion coefficient

# Grid of positions, centered so that X = 0 sits in the middle.
N = int(L / dX) + 1                 # number of grid points
X = np.linspace(-L / 2, L / 2, N)   # positions from -50 to 50
i0 = np.argmin(np.abs(X))           # index of X = 0

# Stability factor for the explicit scheme (must be <= 0.5 for stability).
alpha = D * dt / dX**2
print(f"Stability factor alpha = D*dt/dX^2 = {alpha}")

# --- Initial condition: all probability concentrated at X = 0 ---
# Represent a point mass as unit probability in one cell (area = P*dX = 1).
P = np.zeros(N)
P[i0] = 1.0 / dX   # so that sum(P)*dX = 1

# --- Time integration ---
t_final = 100.0
n_steps = int(t_final / dt)

# Helper to compute moments of the (normalized) distribution.
def moments(P):
    total = np.sum(P) * dX                 # total probability
    if total <= 0:
        return total, 0.0, 0.0
    mean = np.sum(X * P) * dX / total       # mean position
    var = np.sum((X - mean)**2 * P) * dX / total  # variance
    return total, mean, var

# Times at which to snapshot the profile for the spreading plot.
snapshot_times = [0.0, 5.0, 20.0, 50.0, 100.0]
snapshot_steps = {int(round(t / dt)): t for t in snapshot_times}
snapshots = {}

# Record time series of the tracked quantities.
ts, totals, means, variances = [], [], [], []

# Capture t = 0 before stepping.
if 0 in snapshot_steps:
    snapshots[0.0] = P.copy()
tot0, mean0, var0 = moments(P)
ts.append(0.0); totals.append(tot0); means.append(mean0); variances.append(var0)

for step in range(1, n_steps + 1):
    # Explicit update of interior points using the discrete Laplacian.
    lap = np.zeros(N)
    lap[1:-1] = P[2:] - 2.0 * P[1:-1] + P[:-2]
    P[1:-1] = P[1:-1] + alpha * lap
    # Dirichlet (absorbing) boundary conditions: ends pinned to 0.
    P[0] = 0.0
    P[-1] = 0.0

    t = step * dt
    tot, mean, var = moments(P)
    ts.append(t); totals.append(tot); means.append(mean); variances.append(var)

    if step in snapshot_steps:
        snapshots[snapshot_steps[step]] = P.copy()

ts = np.array(ts)
totals = np.array(totals)
means = np.array(means)
variances = np.array(variances)

# --- Numerical summary at selected times ---
for t in snapshot_times:
    idx = np.argmin(np.abs(ts - t))
    print(f"t = {ts[idx]:6.2f} | total P = {totals[idx]:.6f} | "
          f"mean = {means[idx]:+.4f} | variance = {variances[idx]:.4f} | "
          f"theory 2*D*t = {2*D*ts[idx]:.4f}")

# Early-time agreement of variance with theory (before boundaries bite).
mask_early = ts <= 20.0
rel_err_early = np.max(np.abs(variances[mask_early] - 2*D*ts[mask_early])
                       / (2*D*ts[mask_early] + 1e-12))
print(f"Max relative variance error (t <= 20) = {rel_err_early:.4e}")
print(f"Final total probability at t = {ts[-1]:.2f} = {totals[-1]:.6f} "
      f"(leaked = {1.0 - totals[-1]:.6f})")

# --- Plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Left: P(X) spreading over time.
for t in snapshot_times:
    if t in snapshots:
        ax1.plot(X, snapshots[t], label=f"t = {t:g}")
ax1.set_xlabel("X")
ax1.set_ylabel("P(X)")
ax1.set_title("Diffusion of a point distribution")
ax1.set_xlim(-40, 40)
ax1.legend()

# Right: variance vs t against theory 2*D*t.
ax2.plot(ts, variances, "b-", label="numerical variance")
ax2.plot(ts, 2*D*ts, "r--", label="theory 2*D*t")
ax2.set_xlabel("t")
ax2.set_ylabel("variance")
ax2.set_title("Variance growth vs theory")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.2.1_s1.png")

# --- One-sentence explanation of the check ---
print("Check explanation: The point mass spreads into a Gaussian whose "
      "variance tracks 2*D*t at early times (confirming the integrator "
      "reproduces the diffusion law) while total probability stays ~1 until "
      "the spreading distribution reaches the absorbing Dirichlet ends and "
      "slowly leaks out, exactly the two behaviors the theory predicts.")
