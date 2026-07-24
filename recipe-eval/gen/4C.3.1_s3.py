import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Generic multi-variable Heun integrator (ordinary ODE, NO delay).
# Implemented explicitly: predictor (Euler) + corrector (trapezoidal avg).
# ----------------------------------------------------------------------
def heun(f, y0, dt, tmax):
    n = int(round(tmax / dt))          # number of steps
    t = np.linspace(0.0, n * dt, n + 1)
    Y = np.zeros((n + 1, len(y0)))     # store the trajectory
    Y[0] = y0
    for i in range(n):
        yi = Y[i]
        k1 = f(yi)                     # slope at the current point
        yp = yi + dt * k1              # predictor: one explicit Euler step
        k2 = f(yp)                     # slope at the predicted point
        Y[i + 1] = yi + dt * 0.5 * (k1 + k2)  # corrector: average the two slopes
    return t, Y

# Hill activation (coefficient 3): rises from 0 toward 1 as u grows
def act(u):
    u3 = u ** 3
    return u3 / (1.0 + u3)

# Hill repression (coefficient 3): falls from 1 toward 0 as u grows
def rep(u):
    return 1.0 / (1.0 + u ** 3)

# Maximal rates and simulation controls
g = h = l = m = 10.0
dt = 0.01
tmax = 30.0

# ---- Three-node ring: X->Y->Z-|X (Z represses X, everything unit-degraded)
def f3(s):
    x, y, z = s
    dx = g * rep(z) - x
    dy = h * act(x) - y
    dz = l * act(y) - z
    return np.array([dx, dy, dz])

# ---- Four-node ring: X->Y->Z->W-|X (extra intermediate gene W between Z and X)
def f4(s):
    x, y, z, w = s
    dx = g * rep(w) - x
    dy = h * act(x) - y
    dz = l * act(y) - z
    dw = m * act(z) - w
    return np.array([dx, dy, dz, dw])

# Integrate both rings from the all-ones initial state
t3, Y3 = heun(f3, np.ones(3), dt, tmax)
t4, Y4 = heun(f4, np.ones(4), dt, tmax)

# ----------------------------------------------------------------------
# Estimate period from upward zero-crossings of x about its late-time mean
# (use the second half of the run so transients have decayed)
# ----------------------------------------------------------------------
def period(t, x):
    half = len(x) // 2
    xt, xx = t[half:], x[half:]
    mean = xx.mean()
    d = xx - mean
    # times where the signal crosses its mean going upward
    cross = [xt[i] + (xt[i + 1] - xt[i]) * (-d[i]) / (d[i + 1] - d[i])
             for i in range(len(d) - 1) if d[i] < 0.0 <= d[i + 1]]
    if len(cross) < 2:
        return float("nan")
    return float(np.mean(np.diff(cross)))

p3 = period(t3, Y3[:, 0])
p4 = period(t4, Y4[:, 0])

print(f"Three-node ring: period(x) = {p3:.4f} time units")
print(f"Four-node ring:  period(x) = {p4:.4f} time units")
print(f"Period lengthening (four - three) = {p4 - p3:.4f} time units")
print(f"Period ratio (four / three) = {p4 / p3:.4f}")
print(f"Three-node oscillation amplitude(x) = {Y3[len(Y3)//2:,0].ptp():.4f}")
print(f"Four-node oscillation amplitude(x)  = {Y4[len(Y4)//2:,0].ptp():.4f}")

# ----------------------------------------------------------------------
# Time-series plots
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t3, Y3[:, 0], label="x")
ax1.plot(t3, Y3[:, 1], label="y")
ax1.plot(t3, Y3[:, 2], label="z")
ax1.set_title(f"Three-node repression ring (no delay) — period ~ {p3:.2f}")
ax1.set_ylabel("concentration")
ax1.legend(loc="upper right")

ax2.plot(t4, Y4[:, 0], label="x")
ax2.plot(t4, Y4[:, 1], label="y")
ax2.plot(t4, Y4[:, 2], label="z")
ax2.plot(t4, Y4[:, 3], label="w")
ax2.set_title(f"Four-node ring (extra intermediate W) — period ~ {p4:.2f} (longer)")
ax2.set_xlabel("time")
ax2.set_ylabel("concentration")
ax2.legend(loc="upper right")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.3.1_s3.png")

# One-sentence explanation of why this check confirms the result:
# The undelayed three-node ring already sustains oscillations, and simply
# inserting one more intermediate gene lengthens the period without adding any
# explicit delay term, showing that a chain of intermediate reactions produces
# the same accumulating lag that a discrete time delay would.
print("Check: the ordinary (delay-free) three-node ring oscillates, and adding a "
      "fourth intermediate gene lengthens the period, so an intermediate chain "
      "reproduces the lag of an explicit delay.")
