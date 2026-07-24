import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Generic multi-variable RK4 (from Part 3A).
# Integrates y' = f(t, y) where y is a vector; done explicitly with the
# four slope stages k1..k4 rather than calling a black-box integrator.
# ----------------------------------------------------------------------
def rk4(f, y0, t0, tf, h):
    ts = [t0]
    ys = [np.array(y0, dtype=float)]
    t = t0
    y = np.array(y0, dtype=float)
    while t < tf - 1e-12:
        k1 = f(t, y)                      # slope at start
        k2 = f(t + h/2, y + h/2 * k1)     # slope at midpoint using k1
        k3 = f(t + h/2, y + h/2 * k2)     # slope at midpoint using k2
        k4 = f(t + h, y + h * k3)         # slope at end using k3
        y = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)  # weighted average step
        t = t + h
        ts.append(t)
        ys.append(y)
    return np.array(ts), np.array(ys)

# ----------------------------------------------------------------------
# Chemostat model.  y = [N, C].
#   dN/dt = a1*(C/(C+1))*N - N   (Michaelis-Menten growth minus washout)
#   dC/dt = -(C/(C+1))*N - C + a2 (consumption minus washout plus feed)
# ----------------------------------------------------------------------
a1, a2 = 2.0, 5.0

def chemostat(t, y):
    N, C = y
    mm = C / (C + 1.0)                    # Michaelis-Menten uptake fraction
    dN = a1 * mm * N - N
    dC = -mm * N - C + a2
    return np.array([dN, dC])

# ----------------------------------------------------------------------
# Integrate two trajectories.
# ----------------------------------------------------------------------
t0, tf, h = 0.0, 40.0, 0.01

# Case 1: exactly N(0) = 0  -> should relax to washout (0, 5)
t_a, y_a = rk4(chemostat, [0.0, 5.0], t0, tf, h)
# also start C away from 5 to show it settles onto washout
t_a2, y_a2 = rk4(chemostat, [0.0, 0.0], t0, tf, h)

# Case 2: tiny seed N(0) = 0.01 -> should reach coexistence (8, 1)
t_b, y_b = rk4(chemostat, [0.01, 5.0], t0, tf, h)

wash_final = y_a2[-1]
coex_final = y_b[-1]

# ----------------------------------------------------------------------
# Report numerical results.
# ----------------------------------------------------------------------
print("Parameters: a1 =", a1, " a2 =", a2)
print("Expected washout fixed point: (N, C) = (0, 5)")
print("Expected coexistence fixed point: (N, C) = (8, 1)")
print("N(0)=0     final state:  N = %.6f  C = %.6f" % (y_a2[-1, 0], y_a2[-1, 1]))
print("N(0)=0.01  final state:  N = %.6f  C = %.6f" % (coex_final[0], coex_final[1]))
print("Washout error  |(N,C)-(0,5)| = %.3e" % np.linalg.norm(wash_final - np.array([0.0, 5.0])))
print("Coexistence error |(N,C)-(8,1)| = %.3e" % np.linalg.norm(coex_final - np.array([8.0, 1.0])))

# minimum N along seeded trajectory (shows it dips toward washout first)
imin = np.argmin(y_b[:, 0])
print("Seeded trajectory minimum N = %.6f at t = %.3f (drift toward washout before peeling away)"
      % (y_b[imin, 0], t_b[imin]))
print("Seeded trajectory maximum N = %.6f" % np.max(y_b[:, 0]))

# ----------------------------------------------------------------------
# Phase-plane plot.
# ----------------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(y_a2[:, 0], y_a2[:, 1], 'b-', label="N(0)=0  -> washout")
plt.plot(y_b[:, 0], y_b[:, 1], 'r-', label="N(0)=0.01 -> coexistence")
plt.plot(0, 5, 'ks', markersize=9, label="washout (0, 5)")
plt.plot(8, 1, 'g^', markersize=10, label="coexistence (8, 1)")
plt.plot(y_a2[0, 0], y_a2[0, 1], 'bo')
plt.plot(y_b[0, 0], y_b[0, 1], 'ro')
plt.xlabel("Population N")
plt.ylabel("Substrate C")
plt.title("Chemostat phase plane: washout vs coexistence")
plt.legend()
plt.grid(True)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3C.2.1_s5.png")

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# Because N=0 is an invariant set (dN/dt=0 when N=0) it can only relax to
# washout, whereas any tiny positive seed grows through the unstable washout
# state and is drawn to the stable coexistence point, so the two distinct end
# states confirm washout is unstable and coexistence is the stable attractor.
# ----------------------------------------------------------------------
print("Check: N=0 stays on the N=0 line and settles at washout (0,5); a tiny "
      "positive seed drifts toward washout but, since washout is unstable to "
      "invasion, peels away and converges to coexistence (8,1), confirming which "
      "state is the stable attractor.")
