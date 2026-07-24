import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Delayed exponential growth:  dN/dt = r * N(t - tau)
# Constant history N(t) = 1 for t <= 0
# ----------------------------------------------------------------------
r     = -1.7      # unstable case
tau   = 1.0       # delay
dt    = 0.01      # time step
t_end = 40.0
D     = int(round(tau / dt))   # delay expressed as a whole number of steps (=100)

n_steps = int(round(t_end / dt))
t = np.linspace(0.0, t_end, n_steps + 1)

# f(N_delayed) = r * N(t - tau); the delay is an exact number of steps,
# so the needed delayed values are always ones we have already stored.
def f(N_delayed):
    return r * N_delayed

# ----------------------------------------------------------------------
# First-order Euler for DDE
# ----------------------------------------------------------------------
N_euler = np.ones(n_steps + 1)   # history N=1 covers the first D indices too
for i in range(n_steps):
    N_delay = N_euler[i - D] if i - D >= 0 else 1.0   # constant history before t=0
    N_euler[i + 1] = N_euler[i] + dt * f(N_delay)     # Euler step

# ----------------------------------------------------------------------
# Second-order Heun (Euler predictor + trapezoidal corrector) for DDE
# ----------------------------------------------------------------------
N_heun = np.ones(n_steps + 1)
for i in range(n_steps):
    # delayed values at the current and next step (whole-step delay -> both stored)
    N_delay_i  = N_heun[i - D]     if i - D     >= 0 else 1.0   # N(t_i - tau)
    N_delay_ip = N_heun[i + 1 - D] if i + 1 - D >= 0 else 1.0   # N(t_{i+1} - tau)

    slope_i  = f(N_delay_i)                       # slope at current step
    # (Euler predictor uses slope_i; not needed for the trapezoid RHS here
    #  because the delayed term at i+1 is already known from storage)
    N_pred   = N_heun[i] + dt * slope_i           # Euler predictor
    slope_ip = f(N_delay_ip)                       # slope at next step (delayed term known)

    # trapezoidal corrector: average the two slopes over the whole step
    N_heun[i + 1] = N_heun[i] + 0.5 * dt * (slope_i + slope_ip)

# ----------------------------------------------------------------------
# Check: agree early, drift apart later, Euler overshooting Heun
# ----------------------------------------------------------------------
i_early = int(round(2.0 / dt))    # sample at t = 2
i_late  = n_steps                 # final time t = 40
diff_early = abs(N_euler[i_early] - N_heun[i_early])
diff_late  = abs(N_euler[i_late]  - N_heun[i_late])
euler_overshoots = abs(N_euler[i_late]) > abs(N_heun[i_late])

print(f"Delay in steps D                         : {D}")
print(f"Euler N at t=2                           : {N_euler[i_early]:.6f}")
print(f"Heun  N at t=2                           : {N_heun[i_early]:.6f}")
print(f"|Euler - Heun| at t=2 (early)            : {diff_early:.6e}")
print(f"Euler N at t=40                          : {N_euler[i_late]:.6f}")
print(f"Heun  N at t=40                          : {N_heun[i_late]:.6f}")
print(f"|Euler - Heun| at t=40 (late)            : {diff_late:.6e}")
print(f"Late difference larger than early        : {diff_late > diff_early}")
print(f"Euler magnitude overshoots Heun at t=40  : {euler_overshoots}")

# ----------------------------------------------------------------------
# Plot both trajectories overlaid
# ----------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(t, N_euler, label="Euler (1st order)", color="tab:red", lw=1.2)
plt.plot(t, N_heun,  label="Heun (2nd order)",  color="tab:blue", lw=1.2)
plt.axhline(0.0, color="gray", lw=0.6, ls="--")
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth dN/dt = r N(t-tau), r=-1.7 (unstable)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.2.1_s4.png")

# One-sentence explanation of the check:
# The check confirms the result because both methods share the same true early
# solution (so they agree while local error is tiny), and since Heun's O(dt^2)
# accuracy tracks the oscillatory growth more faithfully than Euler's O(dt),
# the growing late-time gap with Euler on the outside is exactly the accumulated
# first-order overshoot we expect from the lower-order scheme.
print("Explanation: they agree early because the accumulated error is still "
      "negligible, and Euler's larger first-order error makes its amplitude "
      "overshoot the more accurate Heun solution as the unstable oscillation "
      "grows, which is the expected order-of-accuracy signature.")
