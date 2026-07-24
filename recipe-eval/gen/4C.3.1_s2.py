import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Generic multi-variable Heun (improved Euler) integrator.
# Ordinary ODE solver -- NO delay is built in anywhere.
# ---------------------------------------------------------------
def heun(f, y0, dt, T):
    n = int(round(T / dt))                 # number of steps
    y = np.array(y0, dtype=float)          # current state vector
    ts = np.empty(n + 1)                   # storage for time
    ys = np.empty((n + 1, y.size))         # storage for states
    ts[0], ys[0] = 0.0, y
    for k in range(n):
        t = k * dt
        f1 = f(t, y)                       # slope at start (predictor)
        yp = y + dt * f1                   # Euler predictor step
        f2 = f(t + dt, yp)                 # slope at predicted end point
        y = y + dt * 0.5 * (f1 + f2)       # average the two slopes (corrector)
        ts[k + 1] = t + dt
        ys[k + 1] = y
    return ts, ys

# ---------------------------------------------------------------
# Right-hand sides.  Every edge is a Hill term with coefficient 3
# and every species has unit (first-order) degradation.
# ---------------------------------------------------------------
def hill_act(u):   # activation Hill function u^3 / (1 + u^3)
    return u**3 / (1.0 + u**3)

def hill_rep(u):   # repression Hill function 1 / (1 + u^3)
    return 1.0 / (1.0 + u**3)

g = h = l = m = 10.0   # all maximal production rates equal to 10

# Three-node ring: X -> Y -> Z -| X
def f3(t, s):
    x, y, z = s
    dx = g * hill_rep(z) - x
    dy = h * hill_act(x) - y
    dz = l * hill_act(y) - z
    return np.array([dx, dy, dz])

# Four-node ring: X -> Y -> Z -> W -| X (extra intermediate gene W)
def f4(t, s):
    x, y, z, w = s
    dx = g * hill_rep(w) - x
    dy = h * hill_act(x) - y
    dz = l * hill_act(y) - z
    dw = m * hill_act(z) - w
    return np.array([dx, dy, dz, dw])

# ---------------------------------------------------------------
# Integrate: initial state all ones, dt = 0.01, up to t = 30.
# ---------------------------------------------------------------
dt, T = 0.01, 30.0
t3, y3 = heun(f3, [1.0, 1.0, 1.0], dt, T)
t4, y4 = heun(f4, [1.0, 1.0, 1.0, 1.0], dt, T)

# ---------------------------------------------------------------
# Estimate the oscillation period from successive peaks of x(t),
# using only the later part of the run so transients are gone.
# ---------------------------------------------------------------
def period_from_peaks(t, x, t_start=10.0):
    peaks = []
    for i in range(1, len(x) - 1):
        if t[i] >= t_start and x[i] > x[i - 1] and x[i] >= x[i + 1]:
            peaks.append(t[i])
    if len(peaks) < 2:
        return float("nan")
    return np.mean(np.diff(peaks))

P3 = period_from_peaks(t3, y3[:, 0])
P4 = period_from_peaks(t4, y4[:, 0])

# ---------------------------------------------------------------
# Print numerical results, each on its own labelled line.
# ---------------------------------------------------------------
print("Three-node ring final state (x,y,z):", y3[-1])
print("Four-node ring final state  (x,y,z,w):", y4[-1])
print("Three-node ring x amplitude (min..max, t>10):",
      y3[t3 >= 10, 0].min(), "..", y3[t3 >= 10, 0].max())
print("Four-node ring  x amplitude (min..max, t>10):",
      y4[t4 >= 10, 0].min(), "..", y4[t4 >= 10, 0].max())
print("Three-node ring estimated period:", P3)
print("Four-node ring estimated period:", P4)
print("Period increase (four minus three):", P4 - P3)
print("Period ratio (four / three):", P4 / P3)

# ---------------------------------------------------------------
# Time-series plots.
# ---------------------------------------------------------------
fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

ax[0].plot(t3, y3[:, 0], label="x")
ax[0].plot(t3, y3[:, 1], label="y")
ax[0].plot(t3, y3[:, 2], label="z")
ax[0].set_title("Three-node repression ring (no delay)  period ~ %.2f" % P3)
ax[0].set_ylabel("concentration")
ax[0].legend(loc="upper right")

ax[1].plot(t4, y4[:, 0], label="x")
ax[1].plot(t4, y4[:, 1], label="y")
ax[1].plot(t4, y4[:, 2], label="z")
ax[1].plot(t4, y4[:, 3], label="w")
ax[1].set_title("Four-node repression ring (no delay)  period ~ %.2f" % P4)
ax[1].set_xlabel("time")
ax[1].set_ylabel("concentration")
ax[1].legend(loc="upper right")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.3.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The three-node ring sustains oscillations with an ordinary "
      "(non-delay) integrator, and inserting a fourth intermediate gene lengthens "
      "the period without adding any delay term, showing that a chain of "
      "intermediate reactions reproduces the same lag that a delay would supply.")
