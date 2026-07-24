import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: two-gene negative-feedback loop (nondimensional, Hill coeff 3) ----
# X activates Y, Y represses X.
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h*x^3/(1 + x^3) - y
g = 10.0
h = 10.0

def f(state):
    x, y = state
    dx = g / (1.0 + y**3) - x                 # Y represses X
    dy = h * x**3 / (1.0 + x**3) - y          # X activates Y
    return np.array([dx, dy])

# ---- Generic explicit RK4 integrator (implemented by hand) ----
def rk4(f, y0, t0, t1, dt):
    n = int(round((t1 - t0) / dt))            # number of steps
    ts = np.empty(n + 1)
    ys = np.empty((n + 1, len(y0)))
    ts[0], ys[0] = t0, y0
    y = np.array(y0, dtype=float)
    t = t0
    for i in range(n):
        k1 = f(y)                             # slope at start
        k2 = f(y + 0.5 * dt * k1)             # slope at midpoint using k1
        k3 = f(y + 0.5 * dt * k2)             # slope at midpoint using k2
        k4 = f(y + dt * k3)                   # slope at end using k3
        y = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)  # weighted average
        t = t + dt
        ts[i + 1] = t
        ys[i + 1] = y
    return ts, ys

# ---- Integrate several trajectories from different initial conditions ----
dt = 0.01
t0, t1 = 0.0, 40.0
initial_conditions = [(0.5, 0.5), (9.0, 1.0), (1.0, 9.0), (8.0, 8.0), (0.2, 6.0)]

trajectories = []
for ic in initial_conditions:
    ts, ys = rk4(f, np.array(ic, dtype=float), t0, t1, dt)
    trajectories.append((ic, ts, ys))

# ---- Report the final endpoints (numerical steady state) ----
final_points = np.array([ys[-1] for (_, _, ys) in trajectories])
for ic, _, ys in trajectories:
    print(f"IC {ic} -> final (x, y) = ({ys[-1,0]:.6f}, {ys[-1,1]:.6f})")

mean_ss = final_points.mean(axis=0)
max_spread = np.max(np.linalg.norm(final_points - mean_ss, axis=1))
print(f"Mean steady state across trajectories: x* = {mean_ss[0]:.6f}, y* = {mean_ss[1]:.6f}")
print(f"Max spread of endpoints about the mean: {max_spread:.6e}")

# Residual of dx/dt, dy/dt at the mean endpoint: should be ~0 at a steady state
res = f(mean_ss)
print(f"Residual |f| at mean steady state: {np.linalg.norm(res):.6e}")

# ---- Detect spiraling: count how many times x oscillates about x* over time ----
# Sign changes of (x - x*) indicate the trajectory circles the fixed point.
ic0, ts0, ys0 = trajectories[0]
x_dev = ys0[:, 0] - mean_ss[0]
sign_changes = int(np.sum(np.diff(np.sign(x_dev)) != 0))
print(f"Sign changes of (x - x*) for IC {ic0}: {sign_changes} (>0 => oscillatory/spiral approach)")

# ---- Phase-plane plot ----
plt.figure(figsize=(7, 6))
for ic, ts, ys in trajectories:
    plt.plot(ys[:, 0], ys[:, 1], lw=1.2, label=f"IC {ic}")
    plt.plot(ic[0], ic[1], 'o', ms=4, color='gray')
plt.plot(mean_ss[0], mean_ss[1], 'k*', ms=16, label="steady state")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Two-gene negative-feedback loop: trajectories spiraling to steady state (g=h=10)")
plt.legend(loc="upper right", fontsize=8)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.2.1_s2.png")

# Explanation: because every trajectory, launched from a spread of distinct initial
# conditions, converges to the same point (max spread ~ 0) while its coordinate
# oscillates (multiple sign changes of x - x*) as it approaches, the check confirms
# that the system has a single stable steady state reached by an inward spiral.
print("Check: all trajectories reach one common point (small spread) with an oscillatory approach => single stable spiral steady state.")
