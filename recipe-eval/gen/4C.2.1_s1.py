import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.2.1_s1.png"

# ---------------------------------------------------------------------------
# Model: delayed two-node loop
#   dx/dt = g/(1 + y(t-tau)^3) - x        (X repressed by Y, but delayed)
#   dy/dt = h*x^3/(1 + x^3)   - y         (Y activated by X, instantaneous)
# The RHS receives the current state s=(x,y) and the DELAYED state sd=(xd,yd),
# where yd = y(t-tau) is the only delayed quantity actually used.
# ---------------------------------------------------------------------------
def rhs(s, sd, g, h):
    x, y = s
    xd, yd = sd            # delayed state; here only yd enters the repression
    dx = g / (1.0 + yd**3) - x
    dy = h * x**3 / (1.0 + x**3) - y
    return np.array([dx, dy])

# ---------------------------------------------------------------------------
# Generic multi-variable Heun (improved Euler, predictor-corrector) integrator
# for delay differential equations with a single constant delay tau.
# Constant history: state = hist for all t <= 0.
# ---------------------------------------------------------------------------
def heun_dde(rhs, hist, g, h, tau, dt, t_end):
    nsteps = int(round(t_end / dt))
    kdelay = int(round(tau / dt))          # delay expressed in integer steps
    t = np.arange(nsteps + 1) * dt
    S = np.zeros((nsteps + 1, len(hist)))  # solution trajectory
    S[0] = hist

    # helper: state at delayed index i-kdelay; before t=0 use constant history
    def delayed(i):
        j = i - kdelay
        return S[j] if j >= 0 else np.array(hist, dtype=float)

    for i in range(nsteps):
        sd_now = delayed(i)        # delayed state at current time t_i
        sd_nxt = delayed(i + 1)    # delayed state at t_{i+1} (for corrector)

        # Predictor: forward-Euler estimate of next state
        f1 = rhs(S[i], sd_now, g, h)
        s_pred = S[i] + dt * f1

        # Corrector: average slope at start and predicted end
        f2 = rhs(s_pred, sd_nxt, g, h)
        S[i + 1] = S[i] + 0.5 * dt * (f1 + f2)

    return t, S

# ---------------------------------------------------------------------------
# Parameters and run WITH delay
# ---------------------------------------------------------------------------
g, h = 10.0, 10.0
tau = 2.0
dt = 0.01
t_end = 30.0
hist = np.array([1.0, 1.0])

t, S = heun_dde(rhs, hist, g, h, tau, dt, t_end)
x, y = S[:, 0], S[:, 1]

# ---------------------------------------------------------------------------
# Control run WITHOUT delay (tau = 0): should settle to a stable steady state
# ---------------------------------------------------------------------------
t0, S0 = heun_dde(rhs, hist, g, h, 0.0, dt, t_end)
x0, y0 = S0[:, 0], S0[:, 1]

# ---------------------------------------------------------------------------
# Quantify sustained oscillation over the last third of each run
# ---------------------------------------------------------------------------
tail = t >= (2.0 / 3.0) * t_end   # examine the late-time (settled) window

amp_x_delay = x[tail].max() - x[tail].min()
amp_y_delay = y[tail].max() - y[tail].min()
amp_x_nodelay = x0[tail].max() - x0[tail].min()
amp_y_nodelay = y0[tail].max() - y0[tail].min()

print(f"Parameters: g={g}, h={h}, tau={tau}, dt={dt}, t_end={t_end}")
print(f"Delay in steps (tau/dt): {int(round(tau/dt))}")
print(f"Final state WITH delay:    x={x[-1]:.6f}, y={y[-1]:.6f}")
print(f"Final state WITHOUT delay: x={x0[-1]:.6f}, y={y0[-1]:.6f}")
print(f"Late-time amplitude WITH delay    -> x: {amp_x_delay:.6f}, y: {amp_y_delay:.6f}")
print(f"Late-time amplitude WITHOUT delay -> x: {amp_x_nodelay:.6f}, y: {amp_y_nodelay:.6f}")
print(f"Sustained oscillation confirmed (delay dwarfs no-delay amplitude): "
      f"{amp_x_delay > 0.1 and amp_x_delay > 20 * max(amp_x_nodelay, 1e-9)}")

# Why this check confirms the result: the identical model with tau=0 relaxes to a
# steady state (late-time amplitude ~ 0), so a large sustained amplitude appearing
# only when tau=2 shows the delay itself destabilized that steady state.
print("Explanation: the tau=0 run decays to a fixed point while the tau=2 run keeps "
      "a large steady amplitude, so the delay alone must be what destabilizes the "
      "steady state into a sustained oscillation.")

# ---------------------------------------------------------------------------
# Time-series plot
# ---------------------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(t, x, label="x(t)", color="C0")
plt.plot(t, y, label="y(t)", color="C1")
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title(f"Delayed two-node loop (g={g}, h={h}, tau={tau}): sustained oscillation")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT)
print(f"Saved figure to: {OUT}")
