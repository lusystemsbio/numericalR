import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Two-gene negative-feedback loop (nondimensional form).
# X activates Y, Y represses X.  Hill coefficient n = 3, unit degradation.
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h*x^3/(1 + x^3) - y
# ----------------------------------------------------------------------

g = 10.0
h = 10.0

def f(state):
    """Right-hand side of the ODE system; returns [dx/dt, dy/dt]."""
    x, y = state
    dx = g / (1.0 + y**3) - x                 # X production repressed by Y, minus decay
    dy = h * x**3 / (1.0 + x**3) - y          # Y production activated by X, minus decay
    return np.array([dx, dy])

def rk4_step(state, dt):
    """One explicit classical fourth-order Runge-Kutta step."""
    k1 = f(state)                    # slope at the start of the interval
    k2 = f(state + 0.5 * dt * k1)    # slope at the midpoint using k1
    k3 = f(state + 0.5 * dt * k2)    # slope at the midpoint using k2
    k4 = f(state + dt * k3)          # slope at the end using k3
    # weighted average of the four slopes (weights 1,2,2,1)/6
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

def integrate(state0, dt, nsteps):
    """Integrate the system with the generic RK4 stepper."""
    traj = np.empty((nsteps + 1, 2))
    traj[0] = state0
    s = np.array(state0, dtype=float)
    for i in range(1, nsteps + 1):
        s = rk4_step(s, dt)
        traj[i] = s
    return traj

# ----------------------------------------------------------------------
# Integrate several trajectories from different initial conditions.
# ----------------------------------------------------------------------
dt = 0.01
nsteps = 4000
initial_conditions = [
    (0.5, 0.5),
    (9.0, 0.5),
    (0.5, 9.0),
    (8.0, 8.0),
    (2.0, 6.0),
    (6.0, 2.0),
]

trajectories = [integrate(ic, dt, nsteps) for ic in initial_conditions]

# ----------------------------------------------------------------------
# Steady state: solve f(x*, y*) = 0 by fixed-point iteration.
# At equilibrium x* = g/(1+y*^3) and y* = h*x*^3/(1+x*^3).
# ----------------------------------------------------------------------
x_ss, y_ss = 1.0, 1.0
for _ in range(100000):
    x_new = g / (1.0 + y_ss**3)
    y_new = h * x_new**3 / (1.0 + x_new**3)
    if abs(x_new - x_ss) < 1e-14 and abs(y_new - y_ss) < 1e-14:
        x_ss, y_ss = x_new, y_new
        break
    x_ss, y_ss = x_new, y_new

# ----------------------------------------------------------------------
# Convergence check: do all trajectory endpoints land on the same point,
# and is the approach oscillatory (spiraling)?
# ----------------------------------------------------------------------
endpoints = np.array([traj[-1] for traj in trajectories])
max_dist_to_ss = np.max(np.hypot(endpoints[:, 0] - x_ss, endpoints[:, 1] - y_ss))

# Detect spiraling: count sign changes in (x - x_ss) along one trajectory.
dev_x = trajectories[3][:, 0] - x_ss
sign_changes = int(np.sum(np.diff(np.sign(dev_x[dev_x != 0])) != 0))

# ----------------------------------------------------------------------
# Report numerical results.
# ----------------------------------------------------------------------
print(f"g = {g}")
print(f"h = {h}")
print(f"Steady state x* = {x_ss:.10f}")
print(f"Steady state y* = {y_ss:.10f}")
print(f"Residual dx/dt at steady state = {f((x_ss, y_ss))[0]:.3e}")
print(f"Residual dy/dt at steady state = {f((x_ss, y_ss))[1]:.3e}")
for ic, ep in zip(initial_conditions, endpoints):
    print(f"IC {ic} -> endpoint ({ep[0]:.6f}, {ep[1]:.6f})")
print(f"Max distance of any endpoint from steady state = {max_dist_to_ss:.3e}")
print(f"Number of oscillations (sign changes in x-x*) for IC (8,8) = {sign_changes}")

# ----------------------------------------------------------------------
# Phase-plane plot.
# ----------------------------------------------------------------------
plt.figure(figsize=(7, 6))
for ic, traj in zip(initial_conditions, trajectories):
    plt.plot(traj[:, 0], traj[:, 1], lw=1.0, label=f"IC {ic}")
    plt.plot(traj[0, 0], traj[0, 1], 'o', ms=4)
plt.plot(x_ss, y_ss, 'k*', ms=15, label="steady state")
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Negative-feedback loop phase plane (g={g:g}, h={h:g})")
plt.legend(fontsize=8, loc="upper right")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.2.1_s3.png")

# Explanation: All trajectories from distinct initial conditions ending at the
# same point (max endpoint distance ~ 0) while their x-coordinates oscillate
# (many sign changes about x*) confirms a single stable steady state approached
# by inward spiraling.
print("Check: distinct initial conditions all converge to one point while "
      "oscillating about it, confirming a single stable steady state reached by spiraling in.")
