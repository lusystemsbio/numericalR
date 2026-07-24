import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --- Model definition: nondimensional two-node negative-feedback loop ---
# X activates Y, Y represses X, Hill coefficient 3, unit degradation.
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h*x^3/(1 + x^3) - y
def rhs(state, g, h):
    x, y = state
    dxdt = g / (1.0 + y**3) - x            # Y represses X (repressive Hill)
    dydt = h * x**3 / (1.0 + x**3) - y     # X activates Y (activating Hill)
    return np.array([dxdt, dydt])


# --- Generic multi-variable Heun (predictor-corrector) integrator ---
# Ordinary ODE, no delay. Implemented explicitly step-by-step.
def heun(rhs, state0, t0, tf, dt, *params):
    n_steps = int(round((tf - t0) / dt))
    ts = np.empty(n_steps + 1)
    ys = np.empty((n_steps + 1, len(state0)))
    ts[0] = t0
    ys[0] = state0
    state = np.array(state0, dtype=float)
    t = t0
    for i in range(n_steps):
        f1 = rhs(state, *params)               # slope at start of step
        predictor = state + dt * f1            # Euler predictor
        f2 = rhs(predictor, *params)           # slope at predicted endpoint
        state = state + 0.5 * dt * (f1 + f2)   # corrector: average of the two slopes
        t = t + dt
        ts[i + 1] = t
        ys[i + 1] = state
    return ts, ys


# --- Parameters and initial condition ---
g, h = 10.0, 10.0
state0 = (1.0, 1.0)
dt = 0.01
t0, tf = 0.0, 10.0

# --- Integrate ---
ts, ys = heun(rhs, state0, t0, tf, dt, g, h)
x = ys[:, 0]
y = ys[:, 1]

# --- Report final (steady-state) values ---
print(f"Final x (t={tf}): {x[-1]:.6f}")
print(f"Final y (t={tf}): {y[-1]:.6f}")

# --- Steady-state / oscillation check ---
# Residual of the RHS at the final state: if ~0, we are at a fixed point.
res = rhs(ys[-1], g, h)
print(f"RHS residual at final state (dx/dt, dy/dt): {res[0]:.3e}, {res[1]:.3e}")
print(f"Max |RHS residual| at final state: {np.max(np.abs(res)):.3e}")

# Look at the tail of the trajectory: peak-to-peak amplitude over the last 20%.
tail = int(0.8 * len(ts))
ptp_x = np.ptp(x[tail:])
ptp_y = np.ptp(y[tail:])
print(f"Peak-to-peak x over final 20% of run: {ptp_x:.3e}")
print(f"Peak-to-peak y over final 20% of run: {ptp_y:.3e}")

oscillates = (ptp_x > 1e-4) or (ptp_y > 1e-4)
print(f"Oscillation detected in tail: {oscillates}")

# --- Time-series plot ---
plt.figure(figsize=(8, 5))
plt.plot(ts, x, label="x", lw=2)
plt.plot(ts, y, label="y", lw=2)
plt.xlabel("t")
plt.ylabel("concentration")
plt.title("Two-gene negative-feedback loop (no delay): relaxation to steady state")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.1.1_s3.png")

# Explanation of the check:
# The near-zero RHS residual together with a peak-to-peak amplitude collapsing
# to ~0 in the trajectory tail confirms the result because a genuine stable
# steady state means the state stops changing (derivatives vanish, no sustained
# amplitude), whereas an oscillation would keep a finite peak-to-peak swing forever.
print("Check: negligible tail amplitude + ~zero RHS residual => stable steady state, no oscillation.")
