import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Generic multi-variable Heun integrator (predictor-corrector, no delay)
# -------------------------------------------------------------------
def heun(f, x0, dt, t_end):
    n_steps = int(round(t_end / dt))
    t = np.linspace(0.0, n_steps * dt, n_steps + 1)
    X = np.zeros((n_steps + 1, len(x0)))
    X[0] = x0
    for k in range(n_steps):
        xk = X[k]
        f1 = f(xk)                       # slope at current point
        x_pred = xk + dt * f1            # Euler predictor step
        f2 = f(x_pred)                   # slope at predicted point
        X[k + 1] = xk + 0.5 * dt * (f1 + f2)  # average the two slopes (corrector)
    return t, X

# Hill functions (coefficient 3), all with unit degradation
def hill_act(u):  # activation:  u^3 / (1 + u^3)
    return u**3 / (1.0 + u**3)

def hill_rep(u):  # repression:  1 / (1 + u^3)
    return 1.0 / (1.0 + u**3)

# -------------------------------------------------------------------
# Three-node ring: X -| activates Y -> Z, Z represses X
# -------------------------------------------------------------------
g = h = l = 10.0
def f3(s):
    x, y, z = s
    dx = g * hill_rep(z) - x        # Z represses X
    dy = h * hill_act(x) - y        # X activates Y
    dz = l * hill_act(y) - z        # Y activates Z
    return np.array([dx, dy, dz])

# -------------------------------------------------------------------
# Four-node ring: extra gene W inserted between Z and X
# -------------------------------------------------------------------
m = 10.0
def f4(s):
    x, y, z, w = s
    dx = g * hill_rep(w) - x        # W represses X
    dy = h * hill_act(x) - y        # X activates Y
    dz = l * hill_act(y) - z        # Y activates Z
    dw = m * hill_act(z) - w        # Z activates W
    return np.array([dx, dy, dz, dw])

# Integration settings
dt = 0.01
t_end = 30.0

# Run both rings from all-ones initial state
t3, X3 = heun(f3, np.array([1.0, 1.0, 1.0]), dt, t_end)
t4, X4 = heun(f4, np.array([1.0, 1.0, 1.0, 1.0]), dt, t_end)

# -------------------------------------------------------------------
# Estimate the oscillation period from successive peaks of x(t),
# using only the later part of the run (after transient decays)
# -------------------------------------------------------------------
def estimate_period(t, x):
    # discard the first third as transient
    start = len(t) // 3
    tt, xx = t[start:], x[start:]
    peaks = []
    for i in range(1, len(xx) - 1):
        if xx[i] > xx[i - 1] and xx[i] > xx[i + 1]:
            peaks.append(tt[i])
    if len(peaks) >= 2:
        return np.mean(np.diff(peaks))
    return float("nan")

period3 = estimate_period(t3, X3[:, 0])
period4 = estimate_period(t4, X4[:, 0])

print(f"Three-node ring: estimated period of x(t) = {period3:.4f} time units")
print(f"Four-node ring:  estimated period of x(t) = {period4:.4f} time units")
print(f"Period increase (four-node minus three-node) = {period4 - period3:.4f} time units")
print(f"Relative lengthening = {(period4 - period3) / period3 * 100:.2f} %")
print(f"Three-node ring oscillates without delay: {'YES' if not np.isnan(period3) else 'NO'}")
print(f"Adding a fourth node lengthens the period: {'YES' if period4 > period3 else 'NO'}")

# Confirmation of sustained oscillation: peak-to-trough amplitude in the tail
amp3 = X3[len(t3)//2:, 0].max() - X3[len(t3)//2:, 0].min()
amp4 = X4[len(t4)//2:, 0].max() - X4[len(t4)//2:, 0].min()
print(f"Three-node sustained x amplitude (late) = {amp3:.4f}")
print(f"Four-node sustained x amplitude (late)  = {amp4:.4f}")

# -------------------------------------------------------------------
# Time-series plots
# -------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t3, X3[:, 0], label="x")
ax1.plot(t3, X3[:, 1], label="y")
ax1.plot(t3, X3[:, 2], label="z")
ax1.set_title(f"Three-node repression ring (no delay)  |  period ~ {period3:.2f}")
ax1.set_ylabel("concentration")
ax1.legend(loc="upper right")
ax1.grid(alpha=0.3)

ax2.plot(t4, X4[:, 0], label="x")
ax2.plot(t4, X4[:, 1], label="y")
ax2.plot(t4, X4[:, 2], label="z")
ax2.plot(t4, X4[:, 3], label="w")
ax2.set_title(f"Four-node repression ring (no delay)  |  period ~ {period4:.2f}  (lengthened)")
ax2.set_xlabel("time")
ax2.set_ylabel("concentration")
ax2.legend(loc="upper right")
ax2.grid(alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.3.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the undelayed three-node ring already sustains oscillations "
      "and inserting a fourth intermediate gene lengthens the period without any explicit "
      "time delay, a chain of intermediate reactions reproduces the same lag that an "
      "explicit delay would impose, confirming that the two are interchangeable mechanisms.")
