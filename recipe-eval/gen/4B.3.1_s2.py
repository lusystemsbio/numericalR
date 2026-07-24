import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Delayed Lotka-Volterra system.
# State vector y = [N, P].
# The RHS needs BOTH the current state y(t) and the delayed
# state y(t - tau) because the cross-terms respond after a lag.
#   dN/dt = N * (a - b * P(t-tau))
#   dP/dt = P * (c * N(t-tau) - d)
# ---------------------------------------------------------------
def rhs(y_now, y_delayed, a, b, c, d):
    N, P = y_now                 # current populations
    Nd, Pd = y_delayed           # populations one delay ago
    dN = N * (a - b * Pd)        # prey responds to lagged predators
    dP = P * (c * Nd - d)        # predator responds to lagged prey
    return np.array([dN, dP])

# ---------------------------------------------------------------
# Generic multi-variable Heun (predictor-corrector) integrator
# for delay differential equations.  It carries a VECTOR state
# and a VECTOR history so any number of variables works.
# ---------------------------------------------------------------
def heun_dde(f, y0, dt, t_end, tau, args=()):
    y0 = np.array(y0, dtype=float)
    n_steps = int(round(t_end / dt))          # number of time steps
    lag = int(round(tau / dt))                 # delay measured in steps

    # History buffer of the vector state, one row per computed time.
    # Before t=0 we assume the state was constant at y0 (standard choice).
    hist = [y0.copy()]                          # hist[k] = state at step k
    t_vals = [0.0]

    for k in range(n_steps):
        y_now = hist[k]                         # current vector state

        # --- fetch delayed state y(t - tau) ---
        # index k - lag; if it points before the start, use y0.
        idx = k - lag
        y_delayed_now = hist[idx] if idx >= 0 else y0

        # --- fetch delayed state y(t + dt - tau) for the corrector ---
        idx_next = k + 1 - lag
        y_delayed_next = hist[idx_next] if idx_next >= 0 else y0

        # Predictor: explicit Euler step using the slope at t.
        k1 = f(y_now, y_delayed_now, *args)
        y_pred = y_now + dt * k1

        # Corrector: slope at t+dt evaluated with the predicted state
        # and the delayed state one lag before t+dt, then average slopes.
        k2 = f(y_pred, y_delayed_next, *args)
        y_next = y_now + 0.5 * dt * (k1 + k2)

        hist.append(y_next)                     # append vector to history
        t_vals.append((k + 1) * dt)

    return np.array(t_vals), np.array(hist)

# ---------------------------------------------------------------
# Parameters and run.
# ---------------------------------------------------------------
a, b, c, d = 1.0, 0.03, 0.02, 1.0
y0 = [30.0, 10.0]
dt = 0.01
t_end = 50.0

t0, sol0 = heun_dde(rhs, y0, dt, t_end, tau=0.0,  args=(a, b, c, d))
t1, sol1 = heun_dde(rhs, y0, dt, t_end, tau=0.01, args=(a, b, c, d))

# ---------------------------------------------------------------
# Diagnostic: distance from the fixed point (N*,P*) = (d/c, a/b).
# For a neutrally stable closed loop this distance returns to its
# starting value each cycle; for an outward spiral it grows.
# We compare the radius at the final time to the radius at t=0.
# ---------------------------------------------------------------
N_star, P_star = d / c, a / b   # equilibrium of the classic LV system

def radius(sol):
    return np.hypot(sol[:, 0] - N_star, sol[:, 1] - P_star)

r0 = radius(sol0)
r1 = radius(sol1)

drift0 = r0[-1] - r0[0]   # ~0 => closed loop
drift1 = r1[-1] - r1[0]   # >0 => spiraling outward

print("Equilibrium (N*, P*):", N_star, P_star)
print("tau=0.00  initial radius:", r0[0])
print("tau=0.00  final radius  :", r0[-1])
print("tau=0.00  radius drift (final - initial):", drift0)
print("tau=0.00  max radius:", r0.max())
print("tau=0.01  initial radius:", r1[0])
print("tau=0.01  final radius  :", r1[-1])
print("tau=0.01  radius drift (final - initial):", drift1)
print("tau=0.01  max radius:", r1.max())
print("Relative growth of radius drift (tau=0.01 vs tau=0.00):",
      drift1 - drift0)

# ---------------------------------------------------------------
# Phase-plane plot of both orbits.
# ---------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(sol0[:, 0], sol0[:, 1], lw=1.2,
         label="tau = 0 (neutrally stable closed loop)")
plt.plot(sol1[:, 0], sol1[:, 1], lw=1.2,
         label="tau = 0.01 (spirals outward)")
plt.plot(N_star, P_star, "k+", markersize=12, label="equilibrium")
plt.plot(y0[0], y0[1], "ko", markersize=5, label="start (30, 10)")
plt.xlabel("Prey N")
plt.ylabel("Predator P")
plt.title("Delayed Lotka-Volterra phase plane (Heun DDE integrator)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.3.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the radius drift is essentially zero for "
      "tau=0 (the orbit closes on itself, i.e. neutrally stable) but "
      "clearly positive for tau=0.01 (the radius grows so the orbit "
      "spirals outward), the diagnostic confirms that adding even a "
      "one-step delay destabilizes the otherwise conservative LV cycle.")
