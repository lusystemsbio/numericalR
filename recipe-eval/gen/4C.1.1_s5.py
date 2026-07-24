import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def heun_step(f, state, t, dt):
    """One step of the generic multi-variable Heun (improved Euler) method.

    Heun = predictor (Euler) + corrector (trapezoidal average of slopes).
    Works for any vector-valued RHS f(t, state) -> np.ndarray.
    """
    state = np.asarray(state, dtype=float)
    k1 = np.asarray(f(t, state), dtype=float)          # slope at the start of the step
    predictor = state + dt * k1                         # Euler predictor for the end of the step
    k2 = np.asarray(f(t + dt, predictor), dtype=float)  # slope estimated at the predicted end point
    return state + 0.5 * dt * (k1 + k2)                 # corrector: average the two slopes


def integrate_heun(f, state0, t0, t_end, dt):
    """Integrate f from t0 to t_end using the Heun step above (no delay)."""
    n_steps = int(round((t_end - t0) / dt))            # number of steps to take
    ts = np.empty(n_steps + 1)                          # storage for times
    ys = np.empty((n_steps + 1, len(state0)))           # storage for states
    ts[0] = t0
    ys[0] = np.asarray(state0, dtype=float)
    for i in range(n_steps):                            # march forward step by step
        ys[i + 1] = heun_step(f, ys[i], ts[i], dt)      # advance the state one dt
        ts[i + 1] = ts[i] + dt                          # advance time
    return ts, ys


# Model parameters: nondimensional two-node negative-feedback loop.
# X activates Y and Y represses X, Hill coefficient 3, unit degradation.
g = 10.0
h = 10.0


def rhs(t, state):
    """Right-hand side of the ODE system (no delay: y feeds back instantaneously).

    dx/dt = g / (1 + y^3) - x        (X repressed by Y)
    dy/dt = h * x^3 / (1 + x^3) - y  (Y activated by X)
    """
    x, y = state
    dxdt = g / (1.0 + y**3) - x
    dydt = h * x**3 / (1.0 + x**3) - y
    return np.array([dxdt, dydt])


# Integration settings.
state0 = [1.0, 1.0]   # initial (x, y)
t0 = 0.0
t_end = 10.0
dt = 0.01

# Run the simulation.
ts, ys = integrate_heun(rhs, state0, t0, t_end, dt)
x_traj = ys[:, 0]
y_traj = ys[:, 1]

# Report final (steady-state) values.
x_final = x_traj[-1]
y_final = y_traj[-1]
print(f"Final x at t={t_end}: {x_final:.6f}")
print(f"Final y at t={t_end}: {y_final:.6f}")

# Residual of the RHS at the final point: near zero means we have reached a fixed point.
res = rhs(t_end, [x_final, y_final])
print(f"dx/dt at final point: {res[0]:.6e}")
print(f"dy/dt at final point: {res[1]:.6e}")
print(f"RHS residual norm at final point: {np.linalg.norm(res):.6e}")

# Oscillation check: measure the peak-to-peak variation over the second half of the run.
# If the loop merely relaxes (no oscillation), the late-time signal is essentially flat.
half = len(ts) // 2
x_late_ptp = np.ptp(x_traj[half:])
y_late_ptp = np.ptp(y_traj[half:])
print(f"x peak-to-peak over second half of run: {x_late_ptp:.6e}")
print(f"y peak-to-peak over second half of run: {y_late_ptp:.6e}")
oscillating = (x_late_ptp > 1e-3) or (y_late_ptp > 1e-3)
print(f"Oscillation detected in late-time window: {oscillating}")

# Time-series plot of x and y relaxing to steady state.
plt.figure(figsize=(8, 5))
plt.plot(ts, x_traj, label="x(t)")
plt.plot(ts, y_traj, label="y(t)")
plt.axhline(x_final, color="C0", ls="--", lw=0.8, alpha=0.6)
plt.axhline(y_final, color="C1", ls="--", lw=0.8, alpha=0.6)
plt.xlabel("t")
plt.ylabel("concentration")
plt.title("Two-gene negative-feedback loop (no delay): relaxation to steady state")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.1.1_s5.png")

# Explanation: the near-zero RHS residual together with a flat (essentially zero peak-to-peak)
# late-time signal confirms the trajectory settles onto a single fixed point rather than a limit
# cycle, so without delay the loop relaxes to a stable steady state with no oscillation.
print("Check: a near-zero RHS residual with negligible late-time peak-to-peak variation confirms "
      "the system reaches a single stable fixed point instead of oscillating, so no delay => stable relaxation.")
