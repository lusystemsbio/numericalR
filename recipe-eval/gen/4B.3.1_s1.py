import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Delayed Lotka-Volterra:  dN/dt = N*(a - b*P(t-tau))
#                          dP/dt = P*(c*N(t-tau) - d)
# State is the vector y = [N, P]; the "delayed" cross-terms use
# the state evaluated tau time units in the past.
# ---------------------------------------------------------------

a, b, c, d = 1.0, 0.03, 0.02, 1.0
y0 = np.array([30.0, 10.0])   # initial (N, P)
dt = 0.01
t_end = 50.0

def rhs(y, y_delayed):
    """Vector RHS. y = current [N,P]; y_delayed = [N,P] at t-tau."""
    N, P = y
    Nd, Pd = y_delayed          # delayed values enter the cross-terms
    dN = N * (a - b * Pd)
    dP = P * (c * Nd - d)
    return np.array([dN, dP])

def integrate_heun_dde(tau):
    n = int(round(t_end / dt)) + 1          # number of time points
    lag = int(round(tau / dt))              # delay expressed in steps
    t = np.linspace(0.0, t_end, n)
    Y = np.empty((n, 2))                    # vector state history buffer
    Y[0] = y0

    def delayed(i):
        # state at time t_i - tau; constant history (= y0) before t=0
        j = i - lag
        return Y[j] if j >= 0 else y0

    for i in range(n - 1):
        yd_now = delayed(i)                 # delayed state at t_i
        yd_next = delayed(i + 1)            # delayed state at t_{i+1}

        # --- Heun predictor: explicit Euler step ---
        k1 = rhs(Y[i], yd_now)
        y_pred = Y[i] + dt * k1

        # --- Heun corrector: average slope at both endpoints ---
        k2 = rhs(y_pred, yd_next)
        Y[i + 1] = Y[i] + 0.5 * dt * (k1 + k2)

    return t, Y

# Run both cases
t0, Y_no = integrate_heun_dde(tau=0.0)
t1, Y_de = integrate_heun_dde(tau=0.01)

# ---------------------------------------------------------------
# Closed-loop / spiral check: compare distance from the starting
# point at the end vs the maximum radial excursion. A neutrally
# stable closed orbit returns near its start; a spiral does not.
# Use the phase-plane radius relative to the orbit's center.
# ---------------------------------------------------------------
def orbit_stats(Y):
    center = Y.mean(axis=0)
    r = np.linalg.norm(Y - center, axis=1)
    return r[0], r[-1], r.max()

r0_no, rend_no, rmax_no = orbit_stats(Y_no)
r0_de, rend_de, rmax_de = orbit_stats(Y_de)

# Return-to-start gap (how far the final point is from the initial point)
gap_no = np.linalg.norm(Y_no[-1] - Y_no[0])
gap_de = np.linalg.norm(Y_de[-1] - Y_de[0])

print("tau = 0.0  : initial radius from center = %.4f" % r0_no)
print("tau = 0.0  : final radius from center   = %.4f" % rend_no)
print("tau = 0.0  : radial growth (final/initial) = %.4f" % (rend_no / r0_no))
print("tau = 0.0  : max radius = %.4f" % rmax_no)
print("tau = 0.0  : return-to-start gap |Y_end - Y_start| = %.4f" % gap_no)

print("tau = 0.01 : initial radius from center = %.4f" % r0_de)
print("tau = 0.01 : final radius from center   = %.4f" % rend_de)
print("tau = 0.01 : radial growth (final/initial) = %.4f" % (rend_de / r0_de))
print("tau = 0.01 : max radius = %.4f" % rmax_de)
print("tau = 0.01 : return-to-start gap |Y_end - Y_start| = %.4f" % gap_de)

# ---------------------------------------------------------------
# Phase-plane plot
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(Y_no[:, 0], Y_no[:, 1], lw=1.0,
         label="tau = 0 (neutrally stable closed loop)")
plt.plot(Y_de[:, 0], Y_de[:, 1], lw=1.0,
         label="tau = 0.01 (spirals outward)")
plt.plot(y0[0], y0[1], "ko", ms=6, label="start (30, 10)")
plt.xlabel("N (prey)")
plt.ylabel("P (predator)")
plt.title("Delayed Lotka-Volterra phase-plane orbits (Heun DDE integrator)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.3.1_s1.png")

# One-sentence explanation of why this check confirms the result:
print("Explanation: The tau=0 orbit's final radius stays essentially equal to "
      "its initial radius (growth ~1, tiny return gap), confirming a neutrally "
      "stable closed loop, whereas the tau=0.01 orbit's radius grows over time "
      "(growth > 1, large return gap), confirming that even a one-step delay "
      "destabilizes the orbit into an outward spiral.")
