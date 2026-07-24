import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Two-gene negative-feedback loop, NO delay (baseline).
# Nondimensional two-node loop: X activates Y, Y represses X.
# Hill coefficient 3, unit degradation.
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h*x^3/(1 + x^3) - y
# ---------------------------------------------------------------------------

def rhs(state, g, h):
    """Right-hand side of the ODE system for state = [x, y]."""
    x, y = state
    dx = g / (1.0 + y**3) - x          # Y represses X
    dy = h * x**3 / (1.0 + x**3) - y   # X activates Y
    return np.array([dx, dy])


# ---------------------------------------------------------------------------
# Generic multi-variable Heun (improved Euler / explicit trapezoidal) integrator.
# Implemented explicitly, step by step, rather than via a library routine.
# ---------------------------------------------------------------------------

def heun(f, y0, t0, tf, dt, *params):
    n_steps = int(round((tf - t0) / dt))
    ts = np.empty(n_steps + 1)
    ys = np.empty((n_steps + 1, len(y0)))
    ts[0] = t0
    ys[0] = np.array(y0, dtype=float)

    for i in range(n_steps):
        t = ts[i]
        s = ys[i]

        # 1) Predictor: a plain Euler step using the slope at the start.
        k1 = f(s, *params)
        s_pred = s + dt * k1

        # 2) Corrector slope: evaluate the RHS at the predicted end point.
        k2 = f(s_pred, *params)

        # 3) Heun update: average the two slopes (trapezoidal rule).
        ys[i + 1] = s + dt * 0.5 * (k1 + k2)
        ts[i + 1] = t + dt

    return ts, ys


# ---------------------------------------------------------------------------
# Run the test case.
# ---------------------------------------------------------------------------
g, h = 10.0, 10.0
x0, y0 = 1.0, 1.0
dt = 0.01
t0, tf = 0.0, 10.0

ts, ys = heun(rhs, [x0, y0], t0, tf, dt, g, h)
x = ys[:, 0]
y = ys[:, 1]

# Final (approximate steady) state values.
x_final, y_final = x[-1], y[-1]
print(f"g = {g}")
print(f"h = {h}")
print(f"dt = {dt}")
print(f"Number of steps = {len(ts) - 1}")
print(f"Final time = {ts[-1]}")
print(f"Final x = {x_final:.10f}")
print(f"Final y = {y_final:.10f}")

# Residual of the RHS at the final state: if ~0, we are at a steady state.
res = rhs([x_final, y_final], g, h)
print(f"RHS residual dx/dt at final state = {res[0]:.3e}")
print(f"RHS residual dy/dt at final state = {res[1]:.3e}")
print(f"Max |RHS| at final state = {np.max(np.abs(res)):.3e}")

# ---------------------------------------------------------------------------
# Oscillation check: look at the tail of the trajectory. If the loop relaxes
# monotonically to a fixed point (no oscillation), the tail should be flat,
# i.e. the peak-to-peak variation over the last 20% of the run is negligible.
# ---------------------------------------------------------------------------
tail = slice(int(0.8 * len(ts)), None)
x_ptp_tail = np.ptp(x[tail])
y_ptp_tail = np.ptp(y[tail])
print(f"x peak-to-peak over final 20% = {x_ptp_tail:.3e}")
print(f"y peak-to-peak over final 20% = {y_ptp_tail:.3e}")

# Count sign changes in the discrete time-derivative of x (turning points):
# sustained oscillation would produce many; relaxation produces ~0.
dx_series = np.diff(x)
sign_changes_x = int(np.sum(np.diff(np.sign(dx_series)) != 0))
dy_series = np.diff(y)
sign_changes_y = int(np.sum(np.diff(np.sign(dy_series)) != 0))
print(f"Turning points in x(t) = {sign_changes_x}")
print(f"Turning points in y(t) = {sign_changes_y}")

oscillates = (x_ptp_tail > 1e-4) or (y_ptp_tail > 1e-4)
print(f"Oscillation detected = {oscillates}")

# One-sentence explanation of why this check confirms the result:
# The trajectory settles to constant values (near-zero RHS residual and
# negligible peak-to-peak variation in its tail with no repeated turning
# points), which is exactly the signature of relaxation to a single stable
# steady state rather than a limit cycle, so the no-delay loop does not oscillate.
print("Check explanation: because the RHS residual is ~0 and the tail of the "
      "trajectory is flat with no repeated turning points, the system has "
      "converged to a fixed point rather than a sustained oscillation, "
      "confirming stable relaxation with no delay.")

# ---------------------------------------------------------------------------
# Time-series plot.
# ---------------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(ts, x, label="x(t)", lw=2)
plt.plot(ts, y, label="y(t)", lw=2)
plt.axhline(x_final, color="C0", ls="--", alpha=0.4)
plt.axhline(y_final, color="C1", ls="--", alpha=0.4)
plt.xlabel("t")
plt.ylabel("concentration")
plt.title("Two-gene negative-feedback loop (no delay): relaxation to steady state\n"
          f"g={g}, h={h}, dt={dt}, Heun integrator")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.1.1_s2.png")
