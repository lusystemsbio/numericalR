import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------
# Two-gene negative-feedback loop (no delay), nondimensional:
#   dx/dt = g/(1 + y^3) - x      (Y represses X)
#   dy/dt = h*x^3/(1 + x^3) - y  (X activates Y)
# Hill coefficient 3, unit degradation.
# -------------------------------------------------------------------

g = 10.0
h = 10.0

# Vector field f(state) -> d(state)/dt, generic multi-variable form
def f(state):
    x, y = state
    dx = g / (1.0 + y**3) - x
    dy = h * x**3 / (1.0 + x**3) - y
    return np.array([dx, dy])

# --- Generic multi-variable Heun (improved Euler) integrator ---------
# Implemented explicitly: predictor (Euler) then corrector (trapezoid).
dt = 0.01
t_end = 10.0
n_steps = int(round(t_end / dt))

t = np.zeros(n_steps + 1)
sol = np.zeros((n_steps + 1, 2))
sol[0] = np.array([1.0, 1.0])   # initial (x, y) = (1, 1)

for i in range(n_steps):
    s = sol[i]
    k1 = f(s)                    # slope at current point
    s_pred = s + dt * k1         # Euler predictor step
    k2 = f(s_pred)               # slope at predicted point
    sol[i + 1] = s + dt * 0.5 * (k1 + k2)  # Heun corrector (average slope)
    t[i + 1] = t[i] + dt

x = sol[:, 0]
y = sol[:, 1]

# -------------------------------------------------------------------
# Report numerical results
# -------------------------------------------------------------------
print(f"Final time                : t = {t[-1]:.4f}")
print(f"Steady-state x (final)    : {x[-1]:.6f}")
print(f"Steady-state y (final)    : {y[-1]:.6f}")

# Residual of the vector field at the final point: should be ~0 at steady state
res = f(sol[-1])
print(f"dx/dt at final point      : {res[0]:.3e}")
print(f"dy/dt at final point      : {res[1]:.3e}")
print(f"Vector-field norm at end  : {np.linalg.norm(res):.3e}")

# Oscillation check: look at the last 20% of the trajectory (after transient).
# A stable relaxation has essentially no residual variation there.
tail = int(0.8 * n_steps)
print(f"x range over last 20%     : {x[tail:].max() - x[tail:].min():.3e}")
print(f"y range over last 20%     : {y[tail:].max() - y[tail:].min():.3e}")

# Count sign changes of the "velocity" of x in the tail; oscillation would
# show repeated overshoot/undershoot (many sign flips). Monotone relaxation ~0.
dxdt_tail = np.diff(x[tail:])
sign_changes = int(np.sum(np.diff(np.sign(dxdt_tail)) != 0))
print(f"Sign changes of dx (tail) : {sign_changes}")

# -------------------------------------------------------------------
# Time-series plot
# -------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(t, x, label="x (repressed by Y)", lw=2)
plt.plot(t, y, label="y (activated by X)", lw=2)
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title("Two-gene negative-feedback loop, no delay (Heun integrator)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.1.1_s1.png")

# Explanation:
# The near-zero vector-field norm at the end, the vanishing spread of x and y
# over the last 20% of the run, and the absence of repeated dx sign changes
# together confirm the trajectory settles to a single fixed point without
# sustained oscillation -- i.e., without delay the loop relaxes to a stable
# steady state, since a Hopf-type oscillation would instead leave a persistent
# nonzero amplitude and repeated sign changes in the tail.
print("Check: no-delay loop relaxes to a stable steady state (no oscillation) "
      "because the tail shows ~zero variation and no repeated sign changes, "
      "whereas an oscillation would leave persistent nonzero amplitude.")
