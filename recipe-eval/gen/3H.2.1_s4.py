import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Two-gene negative-feedback loop (X activates Y, Y represses X)
# Nondimensional form, Hill coefficient n = 3, unit degradation:
#   dx/dt = g / (1 + y^3) - x
#   dy/dt = h * x^3 / (1 + x^3) - y
# ---------------------------------------------------------------

g = 10.0  # parameter for X production (repressed by Y)
h = 10.0  # parameter for Y production (activated by X)

def f(state):
    """Right-hand side of the ODE system; returns [dx/dt, dy/dt]."""
    x, y = state
    dxdt = g / (1.0 + y**3) - x
    dydt = h * x**3 / (1.0 + x**3) - y
    return np.array([dxdt, dydt])

def rk4_step(state, dt):
    """One explicit classical Runge-Kutta 4th-order step."""
    k1 = f(state)                    # slope at the start
    k2 = f(state + 0.5 * dt * k1)    # slope at midpoint using k1
    k3 = f(state + 0.5 * dt * k2)    # slope at midpoint using k2
    k4 = f(state + dt * k3)          # slope at the end using k3
    # weighted average of the four slopes
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

def integrate(state0, dt, n_steps):
    """Integrate the system from state0 for n_steps of size dt."""
    traj = np.empty((n_steps + 1, 2))
    traj[0] = state0
    s = np.array(state0, dtype=float)
    for i in range(n_steps):
        s = rk4_step(s, dt)   # advance one RK4 step
        traj[i + 1] = s
    return traj

# ---------------------------------------------------------------
# Integrate several trajectories from different initial conditions
# ---------------------------------------------------------------
dt = 0.01
n_steps = 4000  # total time = 40, long enough to settle

initial_conditions = [
    (0.5, 0.5),
    (8.0, 0.5),
    (0.5, 8.0),
    (9.0, 9.0),
    (2.0, 6.0),
]

trajectories = [integrate(ic, dt, n_steps) for ic in initial_conditions]

# ---------------------------------------------------------------
# Phase-plane plot: trajectories spiraling into the steady state
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
for ic, traj in zip(initial_conditions, trajectories):
    plt.plot(traj[:, 0], traj[:, 1], lw=1.0,
             label=f"start ({ic[0]}, {ic[1]})")
    plt.plot(traj[0, 0], traj[0, 1], 'o', ms=4)  # mark start

# Mark the common endpoint (numerical steady state)
ss = trajectories[0][-1]
plt.plot(ss[0], ss[1], 'k*', ms=15, label="steady state")

plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Two-gene negative-feedback loop (g={g:g}, h={h:g})\n"
          "Trajectories spiraling into a single steady state")
plt.legend(fontsize=8)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.2.1_s4.png")

# ---------------------------------------------------------------
# CHECK: confirm convergence to a single stable steady state,
#        spiraling in.
# ---------------------------------------------------------------

# (a) All trajectories end at the same point -> single steady state.
endpoints = np.array([traj[-1] for traj in trajectories])
mean_endpoint = endpoints.mean(axis=0)
max_endpoint_spread = np.max(np.linalg.norm(endpoints - mean_endpoint, axis=1))

# (b) Residual of the RHS at the endpoint ~ 0 -> it is a fixed point.
residual = np.linalg.norm(f(mean_endpoint))

# (c) Spiraling in: the distance to the steady state oscillates
#     (non-monotone decay) while overall decreasing.
#     Count sign changes in the increments of distance-to-steady-state.
ref_traj = trajectories[4]  # the (2, 6) start shows clear spiraling
dist = np.linalg.norm(ref_traj - mean_endpoint, axis=1)
d_dist = np.diff(dist)
sign_changes = int(np.sum(np.diff(np.sign(d_dist)) != 0))
final_distance = dist[-1]

# ---------------------------------------------------------------
# Print all numerical results
# ---------------------------------------------------------------
print(f"Steady state (mean of endpoints) x*: {mean_endpoint[0]:.6f}")
print(f"Steady state (mean of endpoints) y*: {mean_endpoint[1]:.6f}")
for ic, ep in zip(initial_conditions, endpoints):
    print(f"Endpoint from start {ic}: x={ep[0]:.6f}, y={ep[1]:.6f}")
print(f"Max spread of endpoints (single-attractor check): {max_endpoint_spread:.3e}")
print(f"RHS residual ||f(x*,y*)|| at steady state (fixed-point check): {residual:.3e}")
print(f"Number of oscillations in distance-to-steady-state (spiral check): {sign_changes}")
print(f"Final distance to steady state (reference trajectory): {final_distance:.3e}")

# Explanation of why the check confirms the result:
print("Explanation: Because every trajectory from distinct initial conditions "
      "converges to the same point (tiny endpoint spread) where the RHS "
      "vanishes, while its distance to that point oscillates before decaying, "
      "the check confirms a single stable steady state approached by spiraling in.")
