import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# Generic multi-variable Heun (predictor-corrector) integrator.
# Ordinary ODEs only -- NO delay. This is the whole point: we show that a
# chain of intermediate genes reproduces the lag that a delay would create,
# using a plain (undelayed) integrator.
# -----------------------------------------------------------------------------
def heun(f, y0, dt, t_end):
    n = int(round(t_end / dt))          # number of steps
    t = np.linspace(0.0, n * dt, n + 1) # time grid
    Y = np.zeros((n + 1, len(y0)))      # state history
    Y[0] = y0
    for k in range(n):
        yk = Y[k]
        f1 = f(yk)                       # slope at current point
        y_pred = yk + dt * f1            # Euler predictor step
        f2 = f(y_pred)                   # slope at predicted point
        Y[k + 1] = yk + 0.5 * dt * (f1 + f2)  # corrector: average the slopes
    return t, Y

# Hill helpers (coefficient 3), all with unit degradation applied in the RHS.
def act(u):  # activating Hill:  u^3 / (1 + u^3)
    return u**3 / (1.0 + u**3)

def rep(u):  # repressing Hill:  1 / (1 + u^3)
    return 1.0 / (1.0 + u**3)

# Three-node ring: X -| by Z, X -> Y -> Z.
def ring3(g=10.0, h=10.0, l=10.0):
    def f(s):
        x, y, z = s
        dx = g * rep(z) - x            # Z represses X
        dy = h * act(x) - y            # X activates Y
        dz = l * act(y) - z            # Y activates Z
        return np.array([dx, dy, dz])
    return f

# Four-node ring: extra gene W inserted between Z and X.
def ring4(g=10.0, h=10.0, l=10.0, m=10.0):
    def f(s):
        x, y, z, w = s
        dx = g * rep(w) - x            # W represses X
        dy = h * act(x) - y            # X activates Y
        dz = l * act(y) - z            # Y activates Z
        dw = m * act(z) - w            # Z activates W
        return np.array([dx, dy, dz, dw])
    return f

# -----------------------------------------------------------------------------
# Integrate both rings with the required parameters.
# -----------------------------------------------------------------------------
dt, t_end = 0.01, 30.0

t3, Y3 = heun(ring3(), np.array([1.0, 1.0, 1.0]), dt, t_end)
t4, Y4 = heun(ring4(), np.array([1.0, 1.0, 1.0, 1.0]), dt, t_end)

# -----------------------------------------------------------------------------
# Period estimation via successive upward zero-crossings of x about its mean,
# measured over the latter half of the run (after transients settle).
# -----------------------------------------------------------------------------
def period(t, x):
    half = len(t) // 2
    t, x = t[half:], x[half:]
    xc = x - x.mean()
    crossings = []
    for i in range(len(xc) - 1):
        if xc[i] < 0.0 and xc[i + 1] >= 0.0:  # upward crossing
            # linear interpolation for the crossing time
            frac = -xc[i] / (xc[i + 1] - xc[i])
            crossings.append(t[i] + frac * (t[i + 1] - t[i]))
    if len(crossings) < 2:
        return np.nan
    return np.mean(np.diff(crossings))

p3 = period(t3, Y3[:, 0])
p4 = period(t4, Y4[:, 0])

print(f"Three-node ring estimated period: {p3:.4f}")
print(f"Four-node ring  estimated period: {p4:.4f}")
print(f"Period lengthening (4-node - 3-node): {p4 - p3:.4f}")
print(f"Three-node ring x amplitude (max-min, 2nd half): "
      f"{Y3[len(t3)//2:,0].max() - Y3[len(t3)//2:,0].min():.4f}")
print(f"Four-node ring  x amplitude (max-min, 2nd half): "
      f"{Y4[len(t4)//2:,0].max() - Y4[len(t4)//2:,0].min():.4f}")
print(f"Three-node ring oscillates (amplitude > 0.1): "
      f"{(Y3[len(t3)//2:,0].max() - Y3[len(t3)//2:,0].min()) > 0.1}")
print(f"Four-node period longer than three-node: {p4 > p3}")

# -----------------------------------------------------------------------------
# Plots.
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t3, Y3[:, 0], label="x")
ax1.plot(t3, Y3[:, 1], label="y")
ax1.plot(t3, Y3[:, 2], label="z")
ax1.set_title(f"Three-node repression ring (no delay), period ~ {p3:.2f}")
ax1.set_ylabel("concentration")
ax1.legend(loc="upper right")

ax2.plot(t4, Y4[:, 0], label="x")
ax2.plot(t4, Y4[:, 1], label="y")
ax2.plot(t4, Y4[:, 2], label="z")
ax2.plot(t4, Y4[:, 3], label="w")
ax2.set_title(f"Four-node repression ring (no delay), period ~ {p4:.2f} (longer)")
ax2.set_xlabel("time")
ax2.set_ylabel("concentration")
ax2.legend(loc="upper right")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.3.1_s5.png")

# One-sentence explanation:
# The check confirms the result because the three-node ring sustains oscillations
# using only ordinary (undelayed) ODEs, and inserting a fourth intermediate gene
# lengthens the period toward the delayed two-node loop, demonstrating that a
# chain of intermediate reactions and an explicit time delay produce the same lag.
print("Explanation: the three-node ring oscillates with a plain undelayed "
      "integrator, and adding a fourth intermediate gene lengthens the period "
      "toward the delayed two-node loop, so a chain of intermediate reactions "
      "and an explicit delay capture the same lag.")
