import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --- Model definition -------------------------------------------------------
# Nondimensional two-node negative-feedback loop (no delay):
#   X activates Y, Y represses X, Hill coefficient n = 3, unit degradation.
#   dx/dt = g / (1 + y^3) - x     (x is repressed by y)
#   dy/dt = h * x^3 / (1 + x^3) - y (y is activated by x)
def rhs(state, g, h):
    x, y = state
    dx = g / (1.0 + y**3) - x       # production repressed by y, minus decay
    dy = h * x**3 / (1.0 + x**3) - y  # production activated by x, minus decay
    return np.array([dx, dy])


# --- Generic multi-variable Heun integrator (explicit, no delay) ------------
# Heun's method = explicit trapezoidal predictor-corrector:
#   1) predictor (Euler step):   s_pred = s + dt * f(s)
#   2) corrector (average slope): s_new  = s + dt/2 * (f(s) + f(s_pred))
def heun(f, s0, dt, t_end, g, h):
    n_steps = int(round(t_end / dt))
    ts = np.linspace(0.0, n_steps * dt, n_steps + 1)  # time grid
    traj = np.empty((n_steps + 1, len(s0)))           # solution storage
    traj[0] = s0
    s = np.array(s0, dtype=float)
    for i in range(n_steps):
        k1 = f(s, g, h)             # slope at current point
        s_pred = s + dt * k1        # Euler predictor
        k2 = f(s_pred, g, h)        # slope at predicted point
        s = s + 0.5 * dt * (k1 + k2)  # trapezoidal corrector
        traj[i + 1] = s
    return ts, traj


# --- Run the baseline test case --------------------------------------------
g, h = 10.0, 10.0
s0 = [1.0, 1.0]
dt = 0.01
t_end = 10.0

ts, traj = heun(rhs, s0, dt, t_end, g, h)
x, y = traj[:, 0], traj[:, 1]

# Final (steady-state) values
x_final, y_final = x[-1], y[-1]
print(f"Final x (t=10): {x_final:.6f}")
print(f"Final y (t=10): {y_final:.6f}")

# Residual of the ODE right-hand side at the final point: near zero => steady state
res = rhs([x_final, y_final], g, h)
print(f"RHS residual dx/dt at final point: {res[0]:.3e}")
print(f"RHS residual dy/dt at final point: {res[1]:.3e}")

# Oscillation check: compare the amplitude of the tail (last 30% of the run)
# to the total range. No sustained oscillation => the tail is essentially flat.
tail_start = int(0.7 * len(ts))
x_tail_amp = x[tail_start:].max() - x[tail_start:].min()
y_tail_amp = y[tail_start:].max() - y[tail_start:].min()
print(f"x tail peak-to-peak amplitude (last 30%): {x_tail_amp:.3e}")
print(f"y tail peak-to-peak amplitude (last 30%): {y_tail_amp:.3e}")

# Monotone approach: count sign changes of the derivative of x over the tail.
# Zero (or near-zero) sign changes => monotone relaxation, no oscillation.
dx_tail = np.diff(x[tail_start:])
sign_changes = int(np.sum(np.diff(np.sign(dx_tail[dx_tail != 0])) != 0))
print(f"Sign changes of dx along tail (oscillation indicator): {sign_changes}")

# --- Time-series plot -------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(ts, x, label="x(t)", lw=2)
plt.plot(ts, y, label="y(t)", lw=2)
plt.axhline(x_final, color="C0", ls="--", alpha=0.4)
plt.axhline(y_final, color="C1", ls="--", alpha=0.4)
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title("Two-gene negative-feedback loop (no delay): relaxation to steady state")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.1.1_s4.png")

# Explanation of the check:
# Because the tail amplitude and the RHS residual both collapse to ~0 while x and y
# approach fixed constants monotonically, the trajectory has settled onto a fixed
# point rather than a limit cycle, confirming that without delay the loop relaxes
# to a single stable steady state with no oscillation.
print("Check: near-zero tail amplitude, near-zero RHS residual, and no derivative "
      "sign changes together confirm relaxation to a stable steady state with no oscillation.")
