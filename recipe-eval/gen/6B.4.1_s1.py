import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Setup ---
np.random.seed(12)          # reproducible seed
n_steps = 10000             # number of steps
dt = 1.0                    # time increment per step
dx = (1.0, 1.0)             # step-size scale in (x, y)

# Preallocate arrays for the path; index 0 is the starting point (0, 0)
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)

# --- Explicit 2D Brownian motion ---
# At each time step, take an independent Gaussian step in x and in y.
# x_next = x + N(0,1)*dx_x, y_next = y + N(0,1)*dx_y
for i in range(n_steps):
    x[i + 1] = x[i] + np.random.normal(0.0, 1.0) * dx[0]   # independent x-step
    y[i + 1] = y[i] + np.random.normal(0.0, 1.0) * dx[1]   # independent y-step

# --- Plot the 2D walk path in the plane ---
plt.figure(figsize=(8, 8))
plt.plot(x, y, lw=0.5, color="steelblue")
plt.plot(0, 0, "go", label="start (0,0)")          # start point
plt.plot(x[-1], y[-1], "rs", label="end")          # end point
plt.xlabel("x")
plt.ylabel("y")
plt.title("2D Brownian motion (10000 steps)")
plt.legend()
plt.axis("equal")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.4.1_s1.png")

# --- Check: typical distance from origin grows like sqrt(t) ---
# Distance from origin at each time
t = np.arange(n_steps + 1)
dist = np.sqrt(x**2 + y**2)

# For 2D BM with per-axis variance dt*dx^2, E[dist^2] = (dx_x^2 + dx_y^2) * t.
# So the theoretical RMS distance is sqrt((dx_x^2+dx_y^2) * t) ~ sqrt(t).
final_dist = dist[-1]
theoretical_rms_final = np.sqrt((dx[0]**2 + dx[1]**2) * n_steps)

# Ratio of observed distance to sqrt(t) should be order 1 (not growing/shrinking systematically).
ratio_final = final_dist / np.sqrt(n_steps)

# --- Print numerical results ---
print("Number of steps:", n_steps)
print("Final position x:", x[-1])
print("Final position y:", y[-1])
print("Final distance from origin:", final_dist)
print("Theoretical RMS distance at final step (sqrt((dx_x^2+dx_y^2)*t)):", theoretical_rms_final)
print("Ratio of final distance to sqrt(t):", ratio_final)
print("Max distance from origin over the walk:", dist.max())

# Explanation: The distance from the origin scaling as sqrt(t) (final distance
# close to the theoretical sqrt((dx_x^2+dx_y^2)*t)) confirms diffusive spreading,
# the defining signature of Brownian motion rather than ballistic or static motion.
print("Explanation: distance growing like sqrt(t) confirms diffusive (Brownian) spreading, since a random walk's expected squared displacement is proportional to elapsed time.")
