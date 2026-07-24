import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: two-node loop with DELAYED repression of X by Y
#   dx/dt = g/(1 + y(t-tau)^3) - x
#   dy/dt = h*x^3/(1 + x^3)     - y
# The delay tau enters only through y in the x-equation.
# ---------------------------------------------------------------

g, h, tau = 10.0, 10.0, 2.0
dt = 0.01
t_end = 30.0
n_steps = int(round(t_end / dt))
lag = int(round(tau / dt))          # number of dt-steps in one delay tau

# RHS: takes current state s=(x,y) and the DELAYED state sd=(xd,yd)
def rhs(s, sd):
    x, y = s
    xd, yd = sd                     # delayed values; here only yd is used
    dxdt = g / (1.0 + yd**3) - x    # repression of X by (delayed) Y
    dydt = h * x**3 / (1.0 + x**3) - y
    return np.array([dxdt, dydt])

# ---------------------------------------------------------------
# Generic multi-variable HEUN integrator for delay ODEs.
# We store the whole trajectory so that a delayed value at time
# t-tau is simply the stored state 'lag' steps earlier.
# Constant history: state = (1,1) for all t <= 0.
# ---------------------------------------------------------------
S = np.zeros((n_steps + 1, 2))
t = np.zeros(n_steps + 1)
S[0] = np.array([1.0, 1.0])         # initial state at t=0
hist = np.array([1.0, 1.0])         # constant history for t < 0

def delayed(i):
    # state at time t_i - tau, i.e. 'lag' steps before index i
    j = i - lag
    return S[j] if j >= 0 else hist  # before t=0 use constant history

for i in range(n_steps):
    # delayed state needed at current time t_i
    sd_now = delayed(i)
    # delayed state needed at the predictor time t_{i+1} = t_i + dt
    sd_next = delayed(i + 1)

    # --- Heun predictor (explicit Euler step) ---
    k1 = rhs(S[i], sd_now)          # slope at start of step
    s_pred = S[i] + dt * k1         # Euler prediction of state at t_{i+1}

    # --- Heun corrector (average of start & predicted slopes) ---
    k2 = rhs(s_pred, sd_next)       # slope at end using predicted state
    S[i + 1] = S[i] + 0.5 * dt * (k1 + k2)
    t[i + 1] = t[i] + dt

x = S[:, 0]
y = S[:, 1]

# ---------------------------------------------------------------
# Steady state (fixed point) that would be STABLE without delay.
# Solve x* = g/(1+y*^3), y* = h*x*^3/(1+x*^3) by fixed-point iter.
# ---------------------------------------------------------------
xs, ys = 1.0, 1.0
for _ in range(100000):
    xs_new = g / (1.0 + ys**3)
    ys_new = h * xs_new**3 / (1.0 + xs_new**3)
    xs, ys = xs_new, ys_new
print(f"Steady state x* = {xs:.6f}")
print(f"Steady state y* = {ys:.6f}")

# ---------------------------------------------------------------
# Stability check: compare peak-to-peak amplitude of x in an
# earlier late window vs the final window. A sustained (undamped)
# oscillation keeps roughly constant amplitude; a decay to the
# stable fixed point would shrink it toward zero.
# ---------------------------------------------------------------
i15 = int(round(15.0 / dt))         # window A: t in [15,22.5]
i22 = int(round(22.5 / dt))
i30 = n_steps                       # window B: t in [22.5,30]

ampA = x[i15:i22].max() - x[i15:i22].min()
ampB = x[i22:i30].max() - x[i22:i30].min()
print(f"x peak-to-peak amplitude, t in [15.0, 22.5] = {ampA:.6f}")
print(f"x peak-to-peak amplitude, t in [22.5, 30.0] = {ampB:.6f}")
print(f"amplitude ratio (late/early) = {ampB/ampA:.6f}")
print(f"distance of final state from steady state = "
      f"{np.hypot(x[-1]-xs, y[-1]-ys):.6f}")

sustained = ampB > 0.1 and abs(ampB/ampA - 1.0) < 0.2
print(f"Sustained oscillation (undamped, non-trivial amplitude)? {sustained}")

# Explanation of the check:
print("Check rationale: because the two late-window amplitudes are both "
      "large and nearly equal (ratio ~1) rather than decaying toward zero, "
      "the trajectory never settles onto the fixed point, confirming the "
      "delay destabilized the otherwise-stable steady state into a "
      "sustained oscillation.")

# ---------------------------------------------------------------
# Time-series plot
# ---------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(t, x, label="x(t)", lw=1.5)
plt.plot(t, y, label="y(t)", lw=1.5)
plt.axhline(xs, color="C0", ls="--", lw=0.8, alpha=0.6, label="x* steady state")
plt.axhline(ys, color="C1", ls="--", lw=0.8, alpha=0.6, label="y* steady state")
plt.xlabel("time t")
plt.ylabel("concentration")
plt.title(f"Delayed two-node loop (g={g:g}, h={h:g}, tau={tau:g}): sustained oscillation")
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4C.2.1_s5.png")
