import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Toggle switch (nondimensional), Hill coefficient 3, unit decay:
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h/(1 + x^3) - y
# ---------------------------------------------------------------
g = 5.0
h = 5.0

def f(state):
    """Right-hand side of the ODE system; returns [dx/dt, dy/dt]."""
    x, y = state
    dx = g / (1.0 + y**3) - x
    dy = h / (1.0 + x**3) - y
    return np.array([dx, dy])

# ---------------------------------------------------------------
# Generic RK4 integrator, implemented explicitly step by step.
# ---------------------------------------------------------------
def rk4(rhs, s0, t0, t1, dt):
    n = int(round((t1 - t0) / dt))          # number of steps
    ts = np.empty(n + 1)                    # time samples
    ss = np.empty((n + 1, len(s0)))         # state history
    ts[0] = t0
    ss[0] = s0
    s = np.array(s0, dtype=float)
    t = t0
    for i in range(n):
        k1 = rhs(s)                         # slope at start
        k2 = rhs(s + 0.5 * dt * k1)         # slope at midpoint using k1
        k3 = rhs(s + 0.5 * dt * k2)         # slope at midpoint using k2
        k4 = rhs(s + dt * k3)               # slope at end using k3
        s = s + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)  # weighted average
        t = t + dt
        ts[i+1] = t
        ss[i+1] = s
    return ts, ss

# ---------------------------------------------------------------
# Integrate from several initial conditions.
# ---------------------------------------------------------------
dt = 0.01
T  = 20.0
initial_conditions = [
    (0.5, 4.0),   # biased toward y-high
    (1.0, 3.0),
    (4.0, 0.5),   # biased toward x-high
    (3.0, 1.0),
    (2.0, 0.5),
    (0.5, 2.0),
    (2.5, 2.5),   # near the diagonal
    (2.4, 2.6),
    (2.6, 2.4),
]

trajectories = []
final_states = []
for ic in initial_conditions:
    ts, ss = rk4(f, ic, 0.0, T, dt)
    trajectories.append(ss)
    final_states.append(ss[-1])
    print(f"IC = ({ic[0]:.2f}, {ic[1]:.2f}) -> final (x, y) = ({ss[-1,0]:.6f}, {ss[-1,1]:.6f})")

# ---------------------------------------------------------------
# Separate check: find and classify the steady states.
# A steady state satisfies f(s) = 0.  We locate them by running the
# integrator to convergence from the two representative corners, then
# verify stability via the eigenvalues of the Jacobian.
# ---------------------------------------------------------------
def jacobian(x, y):
    """Analytic Jacobian of the RHS."""
    dfx_dx = -1.0
    dfx_dy = -3.0 * g * y**2 / (1.0 + y**3)**2
    dfy_dx = -3.0 * h * x**2 / (1.0 + x**3)**2
    dfy_dy = -1.0
    return np.array([[dfx_dx, dfx_dy], [dfy_dx, dfy_dy]])

print()
seeds = {"x-high/y-low": (4.0, 0.5), "x-low/y-high": (0.5, 4.0)}
fixed_points = {}
for name, seed in seeds.items():
    _, ss = rk4(f, seed, 0.0, 100.0, dt)   # integrate long enough to settle
    fp = ss[-1]
    fixed_points[name] = fp
    residual = f(fp)                        # should be ~0 at a steady state
    eigvals = np.linalg.eigvals(jacobian(fp[0], fp[1]))
    stable = np.all(eigvals.real < 0)       # stable iff all eigenvalues have negative real part
    print(f"Steady state '{name}': (x, y) = ({fp[0]:.6f}, {fp[1]:.6f})")
    print(f"  RHS residual f(x,y) = ({residual[0]:.2e}, {residual[1]:.2e})")
    print(f"  Jacobian eigenvalues = ({eigvals[0]:.4f}, {eigvals[1]:.4f})")
    print(f"  Stable = {stable}")

# ---------------------------------------------------------------
# Phase-plane plot.
# ---------------------------------------------------------------
plt.figure(figsize=(7, 7))
for ss, ic in zip(trajectories, initial_conditions):
    plt.plot(ss[:, 0], ss[:, 1], '-', lw=1.2, alpha=0.8)
    plt.plot(ic[0], ic[1], 'k.', ms=6)      # start point

for name, fp in fixed_points.items():
    plt.plot(fp[0], fp[1], 'r*', ms=18, zorder=5)
    plt.annotate(name, fp, textcoords="offset points", xytext=(8, 8))

plt.xlabel("x")
plt.ylabel("y")
plt.title("Toggle switch phase plane (g = 5, h = 5): trajectories to two stable states")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.4.1_s1.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# Because integrating from the two corner seeds converges to two
# distinct points that each give a near-zero RHS residual and a
# Jacobian whose eigenvalues are all negative, the check confirms the
# system has two genuinely stable steady states (x-high/y-low and
# x-low/y-high), and the differing final states of the many initial
# conditions show which basin of attraction each trajectory falls into.
# ---------------------------------------------------------------
print()
print("Explanation: The two seeds settle to distinct points with ~0 RHS residual and all-negative Jacobian eigenvalues, confirming two stable steady states, while the varied ICs converging to one or the other confirms the outcome is set by the initial condition (basin of attraction).")
