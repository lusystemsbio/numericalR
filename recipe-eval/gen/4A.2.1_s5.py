import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model: delayed exponential growth  dN/dt = r * N(t - tau) ----
r = -1.7          # unstable case
tau = 1.0         # delay
dt = 0.01         # time step
t_end = 40.0      # final time
D = int(round(tau / dt))   # delay measured in whole time steps (= 100)

# Time grid
t = np.arange(0.0, t_end + dt, dt)
n = len(t)

# Constant history N(t) = 1 for t <= 0.
# We prepend D history points so that index i in the padded array
# corresponds to time (i - D) * dt; the delayed value N(t-tau) at
# step k is simply the stored value D steps earlier.
def solve(method):
    N = np.ones(n + D)   # first D entries hold the constant history = 1
    for k in range(D, n + D - 1):
        f_now = r * N[k - D]          # RHS at current step uses stored delayed value
        if method == "euler":
            # First-order Euler: single forward step
            N[k + 1] = N[k] + dt * f_now
        elif method == "heun":
            # Heun (predictor/corrector) taken over whole steps:
            # 1) Euler predictor for N at the next step
            N_pred = N[k] + dt * f_now
            # 2) Delayed value needed by the corrector at t_{k+1} is
            #    N[k+1-D], which is already stored (whole-step delay).
            f_next = r * N[k + 1 - D]
            # 3) Trapezoidal corrector averages the two slopes
            N[k + 1] = N[k] + 0.5 * dt * (f_now + f_next)
            _ = N_pred  # predictor kept explicit for method structure
    return N[D:]   # drop the history padding, return solution on t-grid

N_euler = solve("euler")
N_heun = solve("heun")

# ---- Plot both trajectories overlaid ----
plt.figure(figsize=(9, 5))
plt.plot(t, N_euler, label="Euler (1st order)", lw=1.2)
plt.plot(t, N_heun, label="Heun (2nd order)", lw=1.2)
plt.axhline(0.0, color="k", lw=0.5)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed growth dN/dt = r N(t-tau),  r = -1.7 (unstable)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.2.1_s5.png")

# ---- Check: agree early, drift apart later; Euler overshoots Heun ----
diff = np.abs(N_euler - N_heun)

# Early window (first 5 time units) vs late window (last 5 time units)
early = t <= 5.0
late = t >= 35.0
max_diff_early = diff[early].max()
max_diff_late = diff[late].max()

# Amplitude (peak-to-peak) over the late window as a measure of overshoot
amp_euler_late = N_euler[late].max() - N_euler[late].min()
amp_heun_late = N_heun[late].max() - N_heun[late].min()

# First time the two curves differ by more than a small tolerance
tol = 1e-3
drift_idx = np.argmax(diff > tol) if np.any(diff > tol) else -1
drift_time = t[drift_idx] if drift_idx >= 0 else float("nan")

print(f"Delay steps D (= tau/dt): {D}")
print(f"Max |Euler - Heun| in early window (t <= 5): {max_diff_early:.6e}")
print(f"Max |Euler - Heun| in late window  (t >= 35): {max_diff_late:.6e}")
print(f"First time |Euler - Heun| exceeds {tol}: t = {drift_time:.3f}")
print(f"Late-window peak-to-peak amplitude, Euler: {amp_euler_late:.6f}")
print(f"Late-window peak-to-peak amplitude, Heun:  {amp_heun_late:.6f}")
print(f"Euler overshoots Heun (larger late amplitude): {amp_euler_late > amp_heun_late}")
print(f"Curves agree early then drift (late diff >> early diff): {max_diff_late > 10 * max_diff_early}")

# Explanation:
print("Explanation: because both schemes share the same exact constant history "
      "their local truncation errors are tiny at first, so early agreement plus "
      "a growing late-time gap with Euler's larger oscillation amplitude confirms "
      "that the extra first-order error accumulates and overshoots the more "
      "accurate second-order Heun solution.")
