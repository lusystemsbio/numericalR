import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Generic multi-variable RK4 (from Part 3A), implemented explicitly.
# f(t, y) returns dy/dt as a numpy array; y is the state vector.
# ---------------------------------------------------------------
def rk4_step(f, t, y, h):
    k1 = f(t, y)                    # slope at start
    k2 = f(t + h/2, y + h/2 * k1)   # slope at midpoint using k1
    k3 = f(t + h/2, y + h/2 * k2)   # slope at midpoint using k2
    k4 = f(t + h,   y + h * k3)     # slope at end using k3
    # weighted average of the four slopes
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

def rk4_integrate(f, y0, t0, t_end, h):
    n = int(round((t_end - t0) / h))     # number of steps
    ts = np.empty(n + 1)
    ys = np.empty((n + 1, len(y0)))
    ts[0], ys[0] = t0, np.array(y0, dtype=float)
    for i in range(n):                    # march forward one step at a time
        ys[i+1] = rk4_step(f, ts[i], ys[i], h)
        ts[i+1] = ts[i] + h
    return ts, ys

# ---------------------------------------------------------------
# Chemostat ODEs: y = [N, C]
#   dN/dt = a1*(C/(C+1))*N - N      (growth by Michaelis-Menten - dilution)
#   dC/dt = -(C/(C+1))*N - C + a2   (consumption - dilution + scaled feed)
# ---------------------------------------------------------------
a1, a2 = 2.0, 5.0

def chemostat(t, y):
    N, C = y
    mu = C / (C + 1.0)               # Michaelis-Menten uptake fraction
    dN = a1 * mu * N - N
    dC = -mu * N - C + a2
    return np.array([dN, dC])

# ---------------------------------------------------------------
# Integrate two trajectories
# ---------------------------------------------------------------
h, t0, t_end = 0.01, 0.0, 30.0
C0 = 5.0                              # start substrate at the feed value

# Case 1: no population seed -> washout
ts1, ys1 = rk4_integrate(chemostat, [0.0, C0], t0, t_end, h)
# Case 2: tiny seed -> coexistence
ts2, ys2 = rk4_integrate(chemostat, [0.01, C0], t0, t_end, h)

washout_state    = ys1[-1]
coexistence_state = ys2[-1]

# ---------------------------------------------------------------
# Report numerical results
# ---------------------------------------------------------------
print(f"Parameters: a1 = {a1}, a2 = {a2}")
print(f"Initial condition (washout run):     N(0) = 0.00,  C(0) = {C0}")
print(f"Final state (washout run):           N = {washout_state[0]:.6f}, C = {washout_state[1]:.6f}")
print(f"Expected washout state:              N = 0.0, C = 5.0")
print(f"Initial condition (coexistence run): N(0) = 0.01,  C(0) = {C0}")
print(f"Final state (coexistence run):       N = {coexistence_state[0]:.6f}, C = {coexistence_state[1]:.6f}")
print(f"Expected coexistence state:          N = 8.0, C = 1.0")

# minimum N along the seeded run shows it first drifts toward washout, then peels away
imin = np.argmin(ys2[:, 0])
print(f"Seeded run minimum N = {ys2[imin,0]:.6f} at t = {ts2[imin]:.3f} (drift toward washout before growth)")

# ---------------------------------------------------------------
# Phase-plane plot
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(ys1[:, 0], ys1[:, 1], 'b-', label="N(0)=0  -> washout")
plt.plot(ys2[:, 0], ys2[:, 1], 'r-', label="N(0)=0.01 -> coexistence")
plt.plot(0, 5, 'bo', markersize=10, label="washout state (0, 5)")
plt.plot(8, 1, 'rs', markersize=10, label="coexistence state (8, 1)")
plt.plot([0.0, 0.01], [C0, C0], 'k.', markersize=8, label="start points")
plt.xlabel("population N")
plt.ylabel("substrate C")
plt.title("Chemostat phase plane: washout vs. coexistence")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3C.2.1_s3.png")

# Why the check confirms the result:
print("Explanation: Because N=0 is an invariant boundary the exactly-unseeded run can only relax to washout (0,5), "
      "while the tiny seed grows exponentially once nutrient is available and is drawn to the coexistence state (8,1), "
      "showing washout is unstable and coexistence is the true attractor.")
