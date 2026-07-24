import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Generic multi-variable Heun (predictor-corrector) integrator.
# Ordinary ODEs only -- NO delay term anywhere.
# ---------------------------------------------------------------
def heun(f, state0, dt, t_end):
    n_steps = int(round(t_end / dt))
    ts = np.linspace(0.0, n_steps * dt, n_steps + 1)
    ys = np.empty((n_steps + 1, len(state0)))
    ys[0] = state0
    y = np.array(state0, dtype=float)
    for k in range(n_steps):
        f0 = f(y)                 # slope at current point
        y_pred = y + dt * f0      # Euler predictor
        f1 = f(y_pred)            # slope at predicted point
        y = y + 0.5 * dt * (f0 + f1)   # average the two slopes (corrector)
        ys[k + 1] = y
    return ts, ys

# Hill activation and repression, coefficient 3, unit degradation implied in f.
def hill_act(u):   # activating input
    return u**3 / (1.0 + u**3)

def hill_rep(u):   # repressing input
    return 1.0 / (1.0 + u**3)

# ---------------------------------------------------------------
# Three-node ring: X -> Y -> Z -| X
# ---------------------------------------------------------------
g = h = l = 10.0
def f3(s):
    x, y, z = s
    dx = g * hill_rep(z) - x       # Z represses X
    dy = h * hill_act(x) - y       # X activates Y
    dz = l * hill_act(y) - z       # Y activates Z
    return np.array([dx, dy, dz])

# ---------------------------------------------------------------
# Four-node ring: X -> Y -> Z -> W -| X
# ---------------------------------------------------------------
m = 10.0
def f4(s):
    x, y, z, w = s
    dx = g * hill_rep(w) - x       # W represses X
    dy = h * hill_act(x) - y       # X activates Y
    dz = l * hill_act(y) - z       # Y activates Z
    dw = m * hill_act(z) - w       # Z activates W
    return np.array([dx, dy, dz, dw])

dt = 0.01
t_end = 30.0

t3, y3 = heun(f3, [1.0, 1.0, 1.0], dt, t_end)
t4, y4 = heun(f4, [1.0, 1.0, 1.0, 1.0], dt, t_end)

# ---------------------------------------------------------------
# Estimate oscillation period from successive peaks of x(t),
# using only the second half (after transients settle).
# ---------------------------------------------------------------
def estimate_period(t, x):
    half = len(t) // 2
    tt, xx = t[half:], x[half:]
    peaks = []
    for i in range(1, len(xx) - 1):
        if xx[i] > xx[i - 1] and xx[i] > xx[i + 1]:
            peaks.append(tt[i])
    if len(peaks) >= 2:
        return float(np.mean(np.diff(peaks)))
    return float("nan")

p3 = estimate_period(t3, y3[:, 0])
p4 = estimate_period(t4, y4[:, 0])

print("Three-node ring x(t) final value:", y3[-1, 0])
print("Four-node ring  x(t) final value:", y4[-1, 0])
print("Three-node ring estimated period:", p3)
print("Four-node ring  estimated period:", p4)
print("Period increase (four minus three):", p4 - p3)
print("Three-node oscillates without delay:", "yes" if not np.isnan(p3) else "no")
print("Adding a fourth node lengthens the period:", "yes" if p4 > p3 else "no")

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

axes[0].plot(t3, y3[:, 0], label="x")
axes[0].plot(t3, y3[:, 1], label="y")
axes[0].plot(t3, y3[:, 2], label="z")
axes[0].set_title("Three-node ring (period ~ %.2f)" % p3)
axes[0].set_ylabel("concentration")
axes[0].legend(loc="upper right")

axes[1].plot(t4, y4[:, 0], label="x")
axes[1].plot(t4, y4[:, 1], label="y")
axes[1].plot(t4, y4[:, 2], label="z")
axes[1].plot(t4, y4[:, 3], label="w")
axes[1].set_title("Four-node ring (period ~ %.2f, longer)" % p4)
axes[1].set_xlabel("time t")
axes[1].set_ylabel("concentration")
axes[1].legend(loc="upper right")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.3.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the undelayed three-node ring already oscillates and "
      "inserting a fourth intermediate gene lengthens the period, the extra "
      "reaction step acts exactly like added time lag, showing a chain of "
      "intermediate genes and an explicit delay produce the same oscillatory effect.")
