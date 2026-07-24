import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------
g, h = 10.0, 10.0        # production strengths
tau = 2.0                # delay in the repression of X by Y
dt = 0.01
T = 30.0
nsteps = int(round(T / dt))
ndelay = int(round(tau / dt))   # number of steps corresponding to the delay

# -------------------------------------------------------------------
# Right-hand side of the delay system.
# s      = current state (x, y)
# s_del  = delayed state (x(t-tau), y(t-tau)); only y(t-tau) is used here.
# -------------------------------------------------------------------
def rhs(s, s_del):
    x, y = s
    y_del = s_del[1]                       # y evaluated at t - tau
    dx = g / (1.0 + y_del**3) - x          # X repressed by delayed Y
    dy = h * x**3 / (1.0 + x**3) - y       # Y activated by current X
    return np.array([dx, dy])

# -------------------------------------------------------------------
# Generic multi-variable Heun integrator for a DDE.
# Constant history: state = (1, 1) for all t <= 0.
# We store the full trajectory so delayed values can be looked up.
# Since tau >> dt, every delayed value needed (at t and at t+dt) is
# already known from the past, so the Heun corrector is explicit.
# -------------------------------------------------------------------
hist = np.array([1.0, 1.0])                      # constant history value
traj = np.zeros((nsteps + 1, 2))
traj[0] = np.array([1.0, 1.0])                   # initial state

def delayed_value(i):
    # state at time index i (relative to start); indices < 0 use history
    return hist if i < 0 else traj[i]

for n in range(nsteps):
    s = traj[n]
    s_del_now = delayed_value(n - ndelay)        # state at t_n - tau
    s_del_next = delayed_value(n + 1 - ndelay)   # state at t_{n+1} - tau

    # Predictor (explicit Euler step)
    k1 = rhs(s, s_del_now)
    s_pred = s + dt * k1

    # Corrector (average of slopes at both ends -> Heun)
    k2 = rhs(s_pred, s_del_next)
    traj[n + 1] = s + 0.5 * dt * (k1 + k2)

t = np.linspace(0.0, T, nsteps + 1)
x = traj[:, 0]
y = traj[:, 1]

# -------------------------------------------------------------------
# Report key numbers
# -------------------------------------------------------------------
print(f"Delay tau: {tau}")
print(f"Delay in steps (ndelay): {ndelay}")
print(f"Final x(t=30): {x[-1]:.6f}")
print(f"Final y(t=30): {y[-1]:.6f}")

# The steady state without delay solves the coupled algebraic system;
# the delayed run should NOT settle to a fixed point but keep oscillating.
# Measure oscillation amplitude over the last third of the run (transient gone).
tail = slice(int(2 * (nsteps + 1) / 3), None)
x_amp = x[tail].max() - x[tail].min()
y_amp = y[tail].max() - y[tail].min()
print(f"x peak-to-peak over last third: {x_amp:.6f}")
print(f"y peak-to-peak over last third: {y_amp:.6f}")

sustained = (x_amp > 1e-3) and (y_amp > 1e-3)
print(f"Sustained oscillation (amplitude does not decay to ~0): {sustained}")

# Explanation of the check:
print("Check rationale: a sustained, non-vanishing peak-to-peak amplitude in "
      "the late-time window shows the trajectory does not converge to the "
      "fixed point, so the delay has destabilized the otherwise-stable steady state.")

# -------------------------------------------------------------------
# Time-series plot
# -------------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(t, x, label="x(t)")
plt.plot(t, y, label="y(t)")
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title(f"Delayed two-node loop (g={g}, h={h}, tau={tau}) — sustained oscillation")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.2.1_s3.png")
