import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Two-node loop with delayed repression (delay differential eqs):
#   dx/dt = g / (1 + y(t-tau)^3) - x
#   dy/dt = h * x^3 / (1 + x^3) - y
# X is repressed by Y, but Y acts on X after a time delay tau.
# ---------------------------------------------------------------

# ---- parameters ----
g = 10.0
h = 10.0
tau = 2.0
dt = 0.01
t_end = 30.0

# ---- generic multi-variable RHS for the DDE ----
# state = [x, y] (current), yd = y value delayed by tau
def rhs(state, yd):
    x, y = state
    dx = g / (1.0 + yd**3) - x
    dy = h * x**3 / (1.0 + x**3) - y
    return np.array([dx, dy])

# ---- time grid ----
n_steps = int(round(t_end / dt))
t = np.linspace(0.0, t_end, n_steps + 1)

# delay expressed in number of steps
delay_steps = int(round(tau / dt))

# ---- storage; constant history (x, y) = (1, 1) for t <= 0 ----
X = np.empty(n_steps + 1)
Y = np.empty(n_steps + 1)
X[0] = 1.0
Y[0] = 1.0

# helper: delayed y value at step index i (i can be negative -> history)
def delayed_y(i):
    j = i - delay_steps
    if j < 0:
        return 1.0            # constant history value
    return Y[j]

# ---- generic multi-variable Heun integrator, done explicitly ----
for i in range(n_steps):
    s = np.array([X[i], Y[i]])

    # delayed y needed at the start (time t_i) and end (time t_{i+1}) of the step
    yd_now = delayed_y(i)          # y(t_i - tau)
    yd_next = delayed_y(i + 1)     # y(t_{i+1} - tau)

    # predictor: explicit Euler step
    k1 = rhs(s, yd_now)
    s_pred = s + dt * k1

    # corrector: slope at predicted end state, then average the two slopes
    k2 = rhs(s_pred, yd_next)
    s_new = s + 0.5 * dt * (k1 + k2)

    X[i + 1] = s_new[0]
    Y[i + 1] = s_new[1]

# ---- report the final state ----
print(f"Final time t = {t[-1]:.2f}")
print(f"x(t_end) = {X[-1]:.6f}")
print(f"y(t_end) = {Y[-1]:.6f}")

# ---- stability check ----------------------------------------------------
# Locate the steady state of the NON-delayed system (tau = 0), where the
# delayed argument y(t-tau) collapses to y(t). At a fixed point dx=dy=0:
#   x* = g / (1 + y*^3)   and   y* = h x*^3 / (1 + x*^3)
# Solve by simple fixed-point iteration.
xs, ys = 1.0, 1.0
for _ in range(100000):
    xs_new = g / (1.0 + ys**3)
    ys_new = h * xs_new**3 / (1.0 + xs_new**3)
    if abs(xs_new - xs) < 1e-14 and abs(ys_new - ys) < 1e-14:
        xs, ys = xs_new, ys_new
        break
    xs, ys = xs_new, ys_new
print(f"Steady state (fixed point): x* = {xs:.6f}, y* = {ys:.6f}")

# Measure the amplitude of oscillation in the second half of the run
# (after transients), to confirm the delay produced a sustained oscillation.
half = n_steps // 2
x_amp = X[half:].max() - X[half:].min()
y_amp = Y[half:].max() - Y[half:].min()
print(f"x oscillation amplitude (last half of run) = {x_amp:.6f}")
print(f"y oscillation amplitude (last half of run) = {y_amp:.6f}")

sustained = x_amp > 1e-2 and y_amp > 1e-2
print(f"Delay destabilizes the stable steady state (sustained oscillation): {sustained}")

# Explanation:
# The steady state x*, y* is unchanged by the delay (the delayed and undelayed
# systems share the same fixed point), so a non-decaying oscillation about that
# same point in the last half of the run confirms the delay turned a stable
# fixed point into an unstable one with a limit cycle.

# ---- time-series plot ----
plt.figure(figsize=(10, 5))
plt.plot(t, X, label="x(t)")
plt.plot(t, Y, label="y(t)")
plt.axhline(xs, color="C0", ls="--", lw=0.8, alpha=0.6, label="x* (steady state)")
plt.axhline(ys, color="C1", ls="--", lw=0.8, alpha=0.6, label="y* (steady state)")
plt.xlabel("t")
plt.ylabel("concentration")
plt.title(f"Delayed two-node repressor loop (g={g}, h={h}, tau={tau})")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.2.1_s2.png")
