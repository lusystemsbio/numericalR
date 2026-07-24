import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Two-node loop with DELAYED repression:
#   dx/dt = g/(1 + y(t-tau)^3) - x
#   dy/dt = h*x^3/(1 + x^3) - y
# X represses itself indirectly through Y, but Y's action on X
# is delayed by tau. We integrate with a generic multi-variable
# Heun (predictor-corrector) scheme adapted for DDEs.
# ---------------------------------------------------------------

# Parameters
g, h = 10.0, 10.0
tau = 2.0
dt = 0.01
t_end = 30.0
n_steps = int(round(t_end / dt))
lag = int(round(tau / dt))          # delay expressed in integer steps

# Generic RHS: f(t, state, delayed_state) -> derivative vector.
# 'state' is (x, y) now; 'delayed_state' is (x, y) evaluated at t-tau.
def f(t, state, delayed):
    x, y = state
    yd = delayed[1]                 # only y(t-tau) enters this model
    dx = g / (1.0 + yd**3) - x
    dy = h * x**3 / (1.0 + x**3) - y
    return np.array([dx, dy])

# Constant history: (x, y) = (1, 1) for all t <= 0
hist = np.array([1.0, 1.0])

# Storage for the full trajectory (index 0 corresponds to t = 0).
X = np.zeros((n_steps + 1, 2))
X[0] = hist.copy()

# Helper: delayed state at step index k (k may be negative -> history).
def delayed_at(k):
    return X[k] if k >= 0 else hist

# ---- Generic multi-variable Heun integration for the DDE ----
for n in range(n_steps):
    t = n * dt
    # delayed states needed at t-tau (predictor) and t+dt-tau (corrector)
    d_now = delayed_at(n - lag)         # state at t - tau
    d_next = delayed_at(n + 1 - lag)    # state at (t+dt) - tau (already known)

    # Predictor (explicit Euler step)
    k1 = f(t, X[n], d_now)
    X_pred = X[n] + dt * k1

    # Corrector (average the slopes)
    k2 = f(t + dt, X_pred, d_next)
    X[n + 1] = X[n] + 0.5 * dt * (k1 + k2)

t_arr = np.linspace(0.0, t_end, n_steps + 1)
x_sol, y_sol = X[:, 0], X[:, 1]

# ---------------------------------------------------------------
# Steady state (fixed point) of the system: identical with or
# without delay, since at a constant state y(t-tau) = y(t).
#   x* = g/(1 + y*^3),  y* = h*x*^3/(1 + x*^3)
# Solve by simple fixed-point iteration.
# ---------------------------------------------------------------
xs, ys = 1.0, 1.0
for _ in range(100000):
    xs_new = g / (1.0 + ys**3)
    ys_new = h * xs_new**3 / (1.0 + xs_new**3)
    if abs(xs_new - xs) < 1e-12 and abs(ys_new - ys) < 1e-12:
        xs, ys = xs_new, ys_new
        break
    xs, ys = xs_new, ys_new

# ---------------------------------------------------------------
# Stability check: compare peak-to-peak amplitude of x in two
# consecutive late windows. A decaying (stable) response would
# shrink toward zero; a sustained oscillation keeps a comparable
# amplitude in both windows.
# ---------------------------------------------------------------
def ptp_in(window):
    lo, hi = window
    mask = (t_arr >= lo) & (t_arr <= hi)
    return x_sol[mask].max() - x_sol[mask].min()

amp_early = ptp_in((20.0, 25.0))
amp_late = ptp_in((25.0, 30.0))

print(f"Delay in steps (tau/dt)         : {lag}")
print(f"Fixed point x*                  : {xs:.6f}")
print(f"Fixed point y*                  : {ys:.6f}")
print(f"x at t=30                       : {x_sol[-1]:.6f}")
print(f"y at t=30                       : {y_sol[-1]:.6f}")
print(f"Peak-to-peak amplitude of x, t in [20,25]: {amp_early:.6f}")
print(f"Peak-to-peak amplitude of x, t in [25,30]: {amp_late:.6f}")
print(f"Amplitude ratio (late/early)    : {amp_late/amp_early:.6f}")
print(f"Sustained oscillation (ratio > 0.8 and amp_late > 0.1): "
      f"{bool(amp_late/amp_early > 0.8 and amp_late > 0.1)}")

# Explanation (one sentence):
print("Explanation: because x* and y* are a genuine fixed point that "
      "the undelayed system settles onto, a non-decaying oscillation of "
      "comparable amplitude in both late windows shows the delay has "
      "turned that previously stable steady state unstable.")

# ---- Time-series plot ----
plt.figure(figsize=(9, 4.5))
plt.plot(t_arr, x_sol, label="x(t)", color="tab:blue")
plt.plot(t_arr, y_sol, label="y(t)", color="tab:red")
plt.axhline(xs, color="tab:blue", ls=":", lw=1, label="x* (steady state)")
plt.axhline(ys, color="tab:red", ls=":", lw=1, label="y* (steady state)")
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title("Delayed two-node loop (g=10, h=10, tau=2): sustained oscillation")
plt.legend(loc="upper right", ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.2.1_s4.png")
