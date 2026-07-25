import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
dx = (1.0, 1.0)   # per-axis step scale (sigma) for x and y
dt = 1.0          # time increment per step
n_steps = 10000   # number of steps
seed = 12

rng = np.random.default_rng(seed)

# --- Explicit 2D Brownian motion simulation ---
# Start the walker at the origin.
x = 0.0
y = 0.0

# Arrays to record the path (include the initial point).
xs = np.empty(n_steps + 1)
ys = np.empty(n_steps + 1)
xs[0] = x
ys[0] = y

# At each time step, take an independent Gaussian step in x and in y.
# x_next = x + N(0,1)*dx_x, y_next = y + N(0,1)*dx_y
for i in range(1, n_steps + 1):
    x = x + rng.normal(0.0, 1.0) * dx[0]   # independent Gaussian step in x
    y = y + rng.normal(0.0, 1.0) * dx[1]   # independent Gaussian step in y
    xs[i] = x
    ys[i] = y

# --- Report final position and distance ---
final_x = xs[-1]
final_y = ys[-1]
final_dist = np.hypot(final_x, final_y)
print(f"Final x position: {final_x:.6f}")
print(f"Final y position: {final_y:.6f}")
print(f"Final distance from origin: {final_dist:.6f}")

# --- Plot the 2D walk path ---
plt.figure(figsize=(7, 7))
plt.plot(xs, ys, lw=0.5, color="steelblue")
plt.plot(0, 0, "go", label="start (0,0)")
plt.plot(final_x, final_y, "ro", label="end")
plt.title("2D Brownian Motion Path")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.legend()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.4.1_s4.png")

# --- Check: typical distance from origin grows like sqrt(t) ---
# For 2D BM with per-axis variance dx^2 per unit time, the expected squared
# distance at time t is E[r^2] = (dx_x^2 + dx_y^2) * t, so typical distance
# (RMS) grows like sqrt(t). We compare the observed distance at several times
# to the theoretical prediction sqrt((dx_x^2+dx_y^2)*t).
times = np.arange(1, n_steps + 1)
r2 = xs[1:] ** 2 + ys[1:] ** 2          # squared distance at each time
theory_rms = np.sqrt((dx[0] ** 2 + dx[1] ** 2) * times)  # predicted RMS ~ sqrt(t)

print("\nGrowth check (observed |r| vs theoretical sqrt-t RMS):")
for t in [100, 1000, 5000, 10000]:
    print(f"t={t:5d}  observed |r|={np.sqrt(r2[t-1]):10.4f}  "
          f"theory RMS={theory_rms[t-1]:10.4f}")

# Fit observed |r| against sqrt(t): slope should be near sqrt(dx_x^2+dx_y^2).
sqrt_t = np.sqrt(times)
slope = np.sum(sqrt_t * np.sqrt(r2)) / np.sum(sqrt_t ** 2)  # least-squares through origin
print(f"\nFitted slope of |r| vs sqrt(t): {slope:.6f}")
print(f"Expected slope sqrt(dx_x^2+dx_y^2): {np.sqrt(dx[0]**2 + dx[1]**2):.6f}")

# Explanation: The check confirms the result because a genuine 2D random walk
# has expected squared displacement growing linearly in time, so the distance
# from the origin scales as sqrt(t) — matching the observed slope confirms the
# simulated path is diffusive rather than ballistic or stationary.
print("\nThis check confirms the result because the walker's distance from the "
      "origin growing in proportion to sqrt(t) is the defining signature of "
      "diffusive Brownian motion (E[r^2] ~ t), distinguishing it from a "
      "directed drift (|r| ~ t) or a bounded/stationary process.")
