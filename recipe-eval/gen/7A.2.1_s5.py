import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 1D explicit finite-difference diffusion integrator (from 7A.1)
#   dP/dt = D * d2P/dX2
# Discretized explicitly (forward Euler in time, centered in space):
#   P_i^{n+1} = P_i^n + (D*dt/dX^2) * (P_{i+1}^n - 2*P_i^n + P_{i-1}^n)
# with Dirichlet ends (P=0 at both boundaries), which absorb probability.
# ---------------------------------------------------------------

# Parameters
L  = 100      # domain length
dX = 1.0      # spatial step
dt = 0.01     # time step
D  = 1.0      # diffusion coefficient
T  = 400.0    # total integration time

# Spatial grid centered so that X = 0 is a grid point
N = int(L / dX) + 1                      # number of grid points
X = np.linspace(-L/2.0, L/2.0, N)        # positions from -50 to +50
i0 = np.argmin(np.abs(X))                # index of X = 0

# Diffusion (CFL) number; must be <= 0.5 for stability of explicit scheme
alpha = D * dt / dX**2
print(f"Diffusion number alpha = D*dt/dX^2 = {alpha:.6f} (stable if <= 0.5)")

# Initial condition: all probability concentrated at X = 0
P = np.zeros(N)
P[i0] = 1.0 / dX          # unit total probability: sum(P*dX) = 1

# Helper: moments of the distribution (probability = P*dX per cell)
def moments(P):
    prob = P * dX
    tot  = prob.sum()                     # total probability
    if tot > 0:
        mean = (X * prob).sum() / tot     # normalized mean
        var  = ((X - mean)**2 * prob).sum() / tot  # normalized variance
    else:
        mean = var = 0.0
    return tot, mean, var

# Time integration, recording moments and a few snapshots
nsteps = int(round(T / dt))
times, totals, means, variances = [], [], [], []
snapshot_times = [0.0, 10.0, 50.0, 100.0, 200.0]
snapshot_steps = {int(round(t/dt)): t for t in snapshot_times}
snapshots = {}

for n in range(nsteps + 1):
    t = n * dt
    # record moments (subsample to keep arrays manageable)
    if n % 100 == 0:
        tot, mean, var = moments(P)
        times.append(t); totals.append(tot); means.append(mean); variances.append(var)
    # save spatial snapshots at requested times
    if n in snapshot_steps:
        snapshots[snapshot_steps[n]] = P.copy()
    if n == nsteps:
        break
    # --- explicit update (one Euler step of the diffusion equation) ---
    lap = np.zeros(N)
    lap[1:-1] = P[2:] - 2*P[1:-1] + P[:-2]   # centered second difference
    P = P + alpha * lap                       # forward-Euler time step
    P[0] = 0.0; P[-1] = 0.0                    # Dirichlet (absorbing) ends

times = np.array(times); totals = np.array(totals)
means = np.array(means); variances = np.array(variances)

# Theoretical variance for free diffusion from a point source
theory_var = 2 * D * times

# ---------------------------------------------------------------
# Numerical results
# ---------------------------------------------------------------
print(f"Grid points N = {N}, X from {X[0]:.1f} to {X[-1]:.1f}, X=0 at index {i0}")
print(f"Initial total probability = {totals[0]:.6f}")
print(f"Initial mean = {means[0]:.6f}")
print(f"Initial variance = {variances[0]:.6f}")

for t_check in [10.0, 50.0, 100.0, 200.0]:
    j = np.argmin(np.abs(times - t_check))
    print(f"t = {times[j]:6.2f} | total P = {totals[j]:.6f} | mean = {means[j]:+.4f} "
          f"| var = {variances[j]:8.4f} | 2Dt = {theory_var[j]:8.4f} "
          f"| rel.err = {abs(variances[j]-theory_var[j])/theory_var[j]*100:6.2f}%")

# Early-time agreement (before absorption at the ends matters)
early_mask = times <= 100.0
early_err = np.mean(np.abs(variances[early_mask][1:] - theory_var[early_mask][1:])
                    / theory_var[early_mask][1:]) * 100
print(f"Mean relative variance error for t <= 100 = {early_err:.3f}%")
print(f"Final total probability at t = {times[-1]:.1f} = {totals[-1]:.6f} (leaked = {1-totals[-1]:.6f})")

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Left: P(X) spreading over time
for t in snapshot_times:
    ax1.plot(X, snapshots[t], label=f"t = {t:.0f}")
ax1.set_xlim(-40, 40)
ax1.set_xlabel("X"); ax1.set_ylabel("P(X)")
ax1.set_title("Point distribution spreading into a Gaussian")
ax1.legend()

# Right: variance vs t against 2*D*t
ax2.plot(times, variances, 'b-', lw=2, label="numerical variance")
ax2.plot(times, theory_var, 'r--', lw=2, label="theory  2*D*t")
ax2.set_xlabel("t"); ax2.set_ylabel("variance $\\sigma^2$")
ax2.set_title("Variance vs time")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.2.1_s5.png", dpi=120)

# ---------------------------------------------------------------
# One-sentence explanation of why this check confirms the result
# ---------------------------------------------------------------
print("Explanation: The numerically integrated point source spreads into a Gaussian "
      "whose measured variance tracks 2*D*t at early times (confirming correct diffusive "
      "dynamics), while the total probability slowly decays as the absorbing Dirichlet "
      "ends remove probability, exactly the behavior expected of the discretized diffusion "
      "equation, so matching both the variance law and the leakage confirms the integrator.")
