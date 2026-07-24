import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Generic multi-variable RK4 (from Part 3A), implemented explicitly ---
def rk4_system(f, y0, t0, tf, h):
    """Integrate dy/dt = f(t, y) for a vector state y using classic RK4.
    Returns arrays of times and states."""
    n_steps = int(round((tf - t0) / h))
    ts = np.empty(n_steps + 1)
    ys = np.empty((n_steps + 1, len(y0)))
    ts[0] = t0
    ys[0] = np.array(y0, dtype=float)
    for i in range(n_steps):
        t = ts[i]
        y = ys[i]
        # four slope estimates
        k1 = f(t, y)                       # slope at start
        k2 = f(t + h/2, y + h/2 * k1)      # slope at midpoint using k1
        k3 = f(t + h/2, y + h/2 * k2)      # slope at midpoint using k2
        k4 = f(t + h,   y + h   * k3)      # slope at end using k3
        # weighted average step
        ys[i+1] = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        ts[i+1] = t + h
    return ts, ys

# --- Chemostat model ---
# dN/dt = a1*(C/(C+1))*N - N      (growth by Michaelis-Menten uptake minus dilution)
# dC/dt = -(C/(C+1))*N - C + a2   (consumption minus dilution plus scaled feed)
a1, a2 = 2.0, 5.0

def chemostat(t, y):
    N, C = y
    mm = C / (C + 1.0)             # Michaelis-Menten saturation factor
    dN = a1 * mm * N - N
    dC = -mm * N - C + a2
    return np.array([dN, dC])

# --- Integrate two trajectories ---
h, t0, tf = 0.001, 0.0, 60.0

# Trajectory 1: no seed population -> should relax to washout (0, 5)
t1, y1 = rk4_system(chemostat, [0.0, a2], t0, tf, h)
# Trajectory 2: tiny seed -> should reach coexistence (8, 1)
t2, y2 = rk4_system(chemostat, [0.01, a2], t0, tf, h)

washout_end = y1[-1]
coexist_end = y2[-1]

print(f"Parameters: a1 = {a1}, a2 = {a2}")
print(f"Expected washout state:     (0, 5)")
print(f"Expected coexistence state: (8, 1)")
print(f"Trajectory from N(0)=0.00 final (N, C):  ({washout_end[0]:.6f}, {washout_end[1]:.6f})")
print(f"Trajectory from N(0)=0.01 final (N, C):  ({coexist_end[0]:.6f}, {coexist_end[1]:.6f})")

# Evidence of "drift toward washout then peel away": minimum N along the seeded run
n_min_idx = np.argmin(y2[:, 0])
print(f"Seeded trajectory minimum N = {y2[n_min_idx,0]:.6f} at t = {t2[n_min_idx]:.3f} (dips low near washout before rising)")
print(f"Seeded trajectory maximum N = {np.max(y2[:,0]):.6f}")

# --- Phase-plane plot ---
plt.figure(figsize=(7, 6))
plt.plot(y1[:, 0], y1[:, 1], 'b-', lw=2, label="N(0)=0  -> washout")
plt.plot(y2[:, 0], y2[:, 1], 'r-', lw=2, label="N(0)=0.01 -> coexistence")
plt.plot(0, 5, 'ks', ms=10, label="washout (0, 5)")
plt.plot(8, 1, 'k^', ms=10, label="coexistence (8, 1)")
plt.plot(y1[0, 0], y1[0, 1], 'bo', ms=6)
plt.plot(y2[0, 0], y2[0, 1], 'ro', ms=6)
plt.xlabel("Population N")
plt.ylabel("Substrate C")
plt.title("Chemostat phase plane: washout vs. coexistence")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3C.2.1_s4.png")

# Explanation: Because N=0 is invariant (dN/dt=0 when N=0), starting with exactly
# zero population can only relax to washout (0,5), whereas an arbitrarily tiny seed
# lets the growth term amplify N away from washout to the coexistence state (8,1) --
# confirming that washout is unstable and coexistence is the stable attractor.
print("Why the check confirms it: N=0 is an invariant set forcing washout, so only a nonzero seed can grow away toward the stable coexistence state, showing washout is unstable and coexistence is the attractor.")
