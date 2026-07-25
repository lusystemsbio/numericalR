import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 1D explicit finite-difference diffusion integrator (from 7A.1)
#   dP/dt = D * d2P/dX2
# solved with forward Euler in time and centered second difference
# in space, using Dirichlet (absorbing) boundaries P=0 at both ends.
# ---------------------------------------------------------------

# Model parameters
L = 100          # domain length
dX = 1.0         # spatial step
dt = 0.01        # time step
D = 1.0          # diffusion coefficient

# Spatial grid centered on 0 so that "X = 0" is a real node
X = np.arange(-L/2, L/2 + dX, dX)   # -50 ... +50
N = len(X)
i0 = np.argmin(np.abs(X))           # index of X = 0

# Stability (CFL) number for the explicit scheme; must be <= 0.5
alpha = D * dt / dX**2
print(f"Stability number alpha = D*dt/dX^2 = {alpha:.4f} (must be <= 0.5)")

# Initial condition: all probability concentrated at X = 0
P = np.zeros(N)
P[i0] = 1.0 / dX          # unit total probability once multiplied by dX

# Time stepping
T_final = 100.0
n_steps = int(round(T_final / dt))

# Times at which we snapshot P(X) for the spreading plot
snapshot_times = [0.0, 1.0, 5.0, 20.0, 50.0, 100.0]
snapshot_steps = {int(round(t/dt)): t for t in snapshot_times}
snapshots = {}

# Records for diagnostics vs time
rec_t, rec_mean, rec_var, rec_total = [], [], [], []

def moments(P):
    # total probability, mean, and variance of the distribution
    total = np.sum(P) * dX
    if total <= 0:
        return total, 0.0, 0.0
    mean = np.sum(X * P) * dX / total
    var = np.sum((X - mean)**2 * P) * dX / total
    return total, mean, var

# store initial snapshot / diagnostics
if 0 in snapshot_steps:
    snapshots[snapshot_steps[0]] = P.copy()
tot, mu, var = moments(P)
rec_t.append(0.0); rec_mean.append(mu); rec_var.append(var); rec_total.append(tot)

for step in range(1, n_steps + 1):
    # centered second spatial difference (interior nodes only)
    lap = np.zeros(N)
    lap[1:-1] = (P[2:] - 2.0*P[1:-1] + P[:-2]) / dX**2
    # explicit forward-Euler update
    P = P + D * dt * lap
    # Dirichlet (absorbing) boundaries: probability at the ends is removed
    P[0] = 0.0
    P[-1] = 0.0

    t = step * dt
    tot, mu, var = moments(P)
    rec_t.append(t); rec_mean.append(mu); rec_var.append(var); rec_total.append(tot)

    if step in snapshot_steps:
        snapshots[snapshot_steps[step]] = P.copy()

rec_t = np.array(rec_t)
rec_mean = np.array(rec_mean)
rec_var = np.array(rec_var)
rec_total = np.array(rec_total)

# ---------------------------------------------------------------
# Compare variance with theory sigma^2 = 2*D*t
# ---------------------------------------------------------------
theory_var = 2.0 * D * rec_t

# Report a few representative comparisons
for t_check in [1.0, 5.0, 20.0, 50.0, 100.0]:
    k = int(round(t_check/dt))
    print(f"t = {rec_t[k]:6.2f}: mean = {rec_mean[k]:+.4e}, "
          f"variance = {rec_var[k]:8.4f}, theory 2Dt = {theory_var[k]:8.4f}, "
          f"total prob = {rec_total[k]:.6f}")

# Early-time agreement metric (before boundary leakage matters)
early = rec_t <= 20.0
rel_err_early = np.abs(rec_var[early] - theory_var[early]) / np.maximum(theory_var[early], 1e-12)
print(f"Max relative variance error for t <= 20 (early times): "
      f"{np.nanmax(rel_err_early[1:]):.4%}")
print(f"Total probability at t = 0   : {rec_total[0]:.6f}")
print(f"Total probability at t = {T_final:.0f} : {rec_total[-1]:.6f} "
      f"(leaked = {rec_total[0]-rec_total[-1]:.6f})")

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# (a) P(X) spreading over time
axL = axes[0]
for t in sorted(snapshots.keys()):
    axL.plot(X, snapshots[t], label=f"t = {t:g}")
axL.set_xlabel("X")
axL.set_ylabel("P(X)")
axL.set_title("Diffusion of a point source: P(X) spreading over time")
axL.set_xlim(-50, 50)
axL.legend()

# (b) variance vs t against 2*D*t
axR = axes[1]
axR.plot(rec_t, rec_var, "b-", label="measured variance")
axR.plot(rec_t, theory_var, "r--", label=r"theory $2Dt$")
axR.set_xlabel("t")
axR.set_ylabel(r"variance $\sigma^2$")
axR.set_title("Variance vs t compared with 2Dt")
axR.legend()
# secondary annotation of total probability leakage
axR2 = axR.twinx()
axR2.plot(rec_t, rec_total, "g:", label="total probability")
axR2.set_ylabel("total probability", color="g")
axR2.tick_params(axis="y", labelcolor="g")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.2.1_s4.png", dpi=120)

# ---------------------------------------------------------------
# One-sentence explanation of why this check confirms the result
# ---------------------------------------------------------------
print("Explanation: The measured variance tracks 2Dt closely at early times "
      "while the profile stays Gaussian and the total probability is nearly "
      "conserved, so the numerical scheme reproduces the analytic diffusion "
      "law; the later downward deviation of the variance coincides with the "
      "slow drop in total probability as the absorbing Dirichlet ends remove "
      "the spreading tails, exactly the expected boundary effect.")
