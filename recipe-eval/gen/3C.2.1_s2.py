import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------------------------------
# Generic multi-variable RK4 (from Part 3A): advances a vector state
# y by one step h for the system dy/dt = f(t, y).
# Implemented explicitly (the four slopes) rather than via a library.
# ------------------------------------------------------------------
def rk4_step(f, t, y, h):
    k1 = f(t, y)                      # slope at the start
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)  # slope at midpoint using k1
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)  # slope at midpoint using k2
    k4 = f(t + h, y + h * k3)         # slope at the end using k3
    return y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted average


def rk4_integrate(f, y0, t0, tf, h):
    n = int(round((tf - t0) / h))     # number of steps
    ts = np.empty(n + 1)
    ys = np.empty((n + 1, len(y0)))
    ts[0], ys[0] = t0, np.asarray(y0, dtype=float)
    for i in range(n):                # march forward one RK4 step at a time
        ts[i + 1] = ts[i] + h
        ys[i + 1] = rk4_step(f, ts[i], ys[i], h)
    return ts, ys


# ------------------------------------------------------------------
# Chemostat model: y = [N, C]
# dN/dt = a1*(C/(C+1))*N - N   (Michaelis-Menten growth minus dilution)
# dC/dt = -(C/(C+1))*N - C + a2 (uptake, washout, scaled feed a2)
# ------------------------------------------------------------------
def make_chemostat(a1, a2):
    def f(t, y):
        N, C = y
        mm = C / (C + 1.0)            # Michaelis-Menten saturation factor
        dN = a1 * mm * N - N
        dC = -mm * N - C + a2
        return np.array([dN, dC])
    return f


# Parameters
a1, a2 = 2.0, 5.0
f = make_chemostat(a1, a2)
t0, tf, h = 0.0, 30.0, 0.001

# Trajectory 1: N(0) = 0 exactly -> should relax to washout (0, 5)
_, y_washout = rk4_integrate(f, [0.0, 5.0], t0, tf, h)

# Trajectory 2: tiny seed N(0) = 0.01 -> should end at coexistence (8, 1)
_, y_coexist = rk4_integrate(f, [0.01, 5.0], t0, tf, h)

end_washout = y_washout[-1]
end_coexist = y_coexist[-1]

# Report final states
print(f"Parameters: a1 = {a1}, a2 = {a2}")
print(f"Final state from N(0)=0    (expect washout    (0, 5)): N = {end_washout[0]:.6f}, C = {end_washout[1]:.6f}")
print(f"Final state from N(0)=0.01 (expect coexistence (8, 1)): N = {end_coexist[0]:.6f}, C = {end_coexist[1]:.6f}")

# Track the "drift toward washout then peel away" for the seeded run:
# find the minimum N reached before it climbs to coexistence.
N_seed = y_coexist[:, 0]
i_min = int(np.argmin(N_seed))
print(f"Seeded run minimum N (drift toward washout): N_min = {N_seed[i_min]:.6f} at t = {i_min*h:.3f}")
print(f"Seeded run peak N (coexistence approach):    N_max = {N_seed.max():.6f}")

# Phase-plane plot
plt.figure(figsize=(7, 6))
plt.plot(y_washout[:, 0], y_washout[:, 1], 'b-', label='N(0)=0  -> washout (0,5)')
plt.plot(y_coexist[:, 0], y_coexist[:, 1], 'r-', label='N(0)=0.01 -> coexistence (8,1)')
plt.plot(0, 5, 'bo', markersize=10, label='washout eq (0,5)')
plt.plot(8, 1, 'r^', markersize=10, label='coexistence eq (8,1)')
plt.plot(y_coexist[0, 0], y_coexist[0, 1], 'k.', markersize=8)
plt.xlabel('Population N')
plt.ylabel('Substrate C')
plt.title('Chemostat phase plane: washout vs. coexistence (a1=2, a2=5)')
plt.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3C.2.1_s2.png")

# Explanation
print("Explanation: Because N=0 is invariant (dN/dt=0 when N=0) the population can")
print("only grow from a nonzero seed, so the seeded run peeling off toward (8,1)")
print("while the exact-zero run stays at (0,5) confirms washout is unstable and")
print("coexistence is the attractor for any positive initial population.")
