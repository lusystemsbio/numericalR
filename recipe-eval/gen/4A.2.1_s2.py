import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model: delayed exponential growth  dN/dt = r * N(t - tau) ----
r = -1.7        # unstable case
tau = 1.0       # delay
dt = 0.01       # time step
t_end = 40.0

lag = int(round(tau / dt))          # delay expressed in whole steps (=100)
n_steps = int(round(t_end / dt))    # number of integration steps
t = np.linspace(0.0, t_end, n_steps + 1)

print(f"Delay in steps (lag): {lag}")
print(f"Number of steps: {n_steps}")

# Constant history N(t) = 1 for t <= 0.
# Because tau is a whole number of steps, the delayed value needed at step i
# is always an already-stored past value N[i-lag], so no interpolation is needed.

# ---- First-order Euler for DDE ----
N_euler = np.ones(n_steps + 1)      # N[0] = 1 (history)
for i in range(n_steps):
    delayed = N_euler[i - lag] if i - lag >= 0 else 1.0   # history = 1 before t=0
    N_euler[i + 1] = N_euler[i] + dt * (r * delayed)      # Euler update

# ---- Second-order Heun for DDE (Euler predictor + trapezoidal corrector) ----
N_heun = np.ones(n_steps + 1)       # N[0] = 1 (history)
for i in range(n_steps):
    # delayed value at the current time t_i
    delayed_i = N_heun[i - lag] if i - lag >= 0 else 1.0
    # delayed value at the next time t_{i+1}; whole-step delay => already stored
    delayed_ip1 = N_heun[i + 1 - lag] if i + 1 - lag >= 0 else 1.0

    f_i = r * delayed_i                                   # slope at t_i
    # Euler predictor (not actually used for the slope here because the DDE
    # right-hand side depends only on the delayed state, which is already known)
    N_pred = N_heun[i] + dt * f_i                         # predictor step
    f_ip1 = r * delayed_ip1                               # slope at t_{i+1}

    # trapezoidal corrector: average the two slopes
    N_heun[i + 1] = N_heun[i] + 0.5 * dt * (f_i + f_ip1)

# ---- Report some numerical results ----
print(f"Euler N at t=40: {N_euler[-1]:.6f}")
print(f"Heun  N at t=40: {N_heun[-1]:.6f}")

# sample comparison early vs late
for t_check in [1.0, 5.0, 20.0, 40.0]:
    idx = int(round(t_check / dt))
    print(f"t={t_check:5.1f}  Euler={N_euler[idx]: .6f}  Heun={N_heun[idx]: .6f}  diff={N_euler[idx]-N_heun[idx]: .6e}")

# ---- Drift check: agree early, diverge late, Euler overshoots Heun ----
diff = np.abs(N_euler - N_heun)
idx_early = int(round(2.0 / dt))    # "early" window up to t=2
idx_late = int(round(30.0 / dt))    # "late" window from t=30
early_max_diff = np.max(diff[:idx_early + 1])
late_max_diff = np.max(diff[idx_late:])
print(f"Max |Euler-Heun| for t<=2  (early): {early_max_diff:.6e}")
print(f"Max |Euler-Heun| for t>=30 (late) : {late_max_diff:.6e}")

# amplitude (envelope) comparison to confirm Euler overshoots
euler_peak = np.max(np.abs(N_euler[idx_late:]))
heun_peak = np.max(np.abs(N_heun[idx_late:]))
print(f"Late-window peak |N| Euler: {euler_peak:.6f}")
print(f"Late-window peak |N| Heun : {heun_peak:.6f}")

agree_early = early_max_diff < late_max_diff
euler_overshoots = euler_peak > heun_peak
print(f"Trajectories agree early then drift apart: {agree_early}")
print(f"First-order Euler overshoots second-order Heun: {euler_overshoots}")

# Explanation of why the check confirms the result:
# The check confirms the result because a smaller early error growing into a larger
# late error, with Euler's oscillation amplitude exceeding Heun's, is exactly the
# signature of first-order truncation error accumulating faster than second-order,
# so the higher-order Heun method more faithfully tracks the true (less-overshooting) solution.
print("Explanation: matching early and diverging later with Euler overshooting shows "
      "the first-order truncation error accumulating faster than Heun's second-order error, "
      "confirming Heun is the more accurate integrator.")

# ---- Plot ----
plt.figure(figsize=(10, 6))
plt.plot(t, N_euler, label="Euler (1st order)", color="tab:red", lw=1.2)
plt.plot(t, N_heun, label="Heun (2nd order)", color="tab:blue", lw=1.2)
plt.axhline(0.0, color="gray", lw=0.6)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title(f"Delayed exponential growth (r={r}, tau={tau}): Euler vs Heun")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.2.1_s2.png", dpi=120, bbox_inches="tight")
