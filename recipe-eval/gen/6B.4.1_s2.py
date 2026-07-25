import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Setup ---
rng = np.random.default_rng(12)          # seed 12 for reproducibility
n_steps = 10000                          # number of time steps
dx = np.array([1.0, 1.0])                # step-size scale per axis (x, y)
dt = 1.0                                 # time increment per step

# Arrays to hold the path; index 0 is the starting point (0, 0).
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)

# --- Explicit Brownian motion via independent per-axis Gaussian steps ---
# At each time step, add an independent N(0,1) draw (scaled by dx) to x and to y.
for i in range(n_steps):
    x[i + 1] = x[i] + rng.normal(0.0, 1.0) * dx[0]   # x_next = x + N(0,1)*dx_x
    y[i + 1] = y[i] + rng.normal(0.0, 1.0) * dx[1]   # y_next = y + N(0,1)*dy_y

# --- Plot the 2D walk path ---
plt.figure(figsize=(8, 8))
plt.plot(x, y, lw=0.5, color="steelblue")
plt.plot(0, 0, "go", label="start (0,0)")
plt.plot(x[-1], y[-1], "ro", label="end")
plt.xlabel("x")
plt.ylabel("y")
plt.title("2D Brownian motion path (10000 steps, seed 12)")
plt.legend()
plt.axis("equal")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.4.1_s2.png")

# --- Basic path results ---
final_dist = np.hypot(x[-1], y[-1])      # straight-line distance of endpoint from origin
print("Final position x:", x[-1])
print("Final position y:", y[-1])
print("Final distance from origin:", final_dist)

# --- Check: typical distance grows like sqrt(t) ---
# For 2D Brownian motion with unit-variance steps per axis, E[r^2] = 2*t,
# so the RMS distance from origin should track sqrt(2*t).
t = np.arange(n_steps + 1)                       # elapsed time at each step (dt = 1)
r2 = x**2 + y**2                                  # squared distance from origin over time
# Compare measured distance to the theoretical sqrt(2*t) at a few checkpoints.
print("\nDistance vs sqrt(2*t) check:")
for tc in [100, 1000, 5000, 10000]:
    measured = np.sqrt(r2[tc])
    theory = np.sqrt(2.0 * tc)
    print(f"  t={tc:5d}: measured r={measured:10.3f}, sqrt(2t)={theory:10.3f}")

# Aggregate check: ratio of mean(r^2) to t should be ~2 (dimension * step variance).
ratio = np.mean(r2[1:] / t[1:])
print("\nMean of r^2 / t over the walk (expected ~2.0):", ratio)

# Explanation: the check confirms the result because a genuine 2D Brownian walk
# has mean-squared displacement growing linearly in time (E[r^2] = 2t), so the
# typical distance rising as sqrt(t) — verified by r^2/t staying near 2 — is the
# diffusive signature that distinguishes a tangled random walk from directed motion.
print("\nExplanation: r grows like sqrt(t) because E[r^2]=2t for independent unit-variance"
      " Gaussian steps per axis, the defining diffusive signature of a 2D random walk.")
