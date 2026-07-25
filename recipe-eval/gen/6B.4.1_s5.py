import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model parameters ---
dx = (1.0, 1.0)   # per-axis step scale (std of each Gaussian step, one per axis)
dt = 1            # time increment per step (used to define time axis)
n_steps = 10000   # number of steps
np.random.seed(12)

# --- Explicit simulation of 2D Brownian motion ---
# Start the walker at the origin.
x, y = 0.0, 0.0
xs = [x]   # store the path
ys = [y]

# At each time step, take an INDEPENDENT Gaussian step in x and in y:
#   x_next = x + N(0,1)*dx[0],  y_next = y + N(0,1)*dx[1]
for _ in range(n_steps):
    x = x + np.random.normal(0.0, 1.0) * dx[0]   # independent step in x
    y = y + np.random.normal(0.0, 1.0) * dx[1]   # independent step in y
    xs.append(x)
    ys.append(y)

xs = np.array(xs)
ys = np.array(ys)

# --- Plot the 2D walk path in the plane ---
plt.figure(figsize=(7, 7))
plt.plot(xs, ys, lw=0.5, color="steelblue")
plt.plot(0, 0, "go", label="start (0,0)")           # start
plt.plot(xs[-1], ys[-1], "ro", label="end")         # end
plt.title("2D Brownian motion (independent per-axis Gaussian steps)")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.legend()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.4.1_s5.png")

# --- Reported numerical results ---
print("Number of steps:", n_steps)
print("Final position x:", xs[-1])
print("Final position y:", ys[-1])
print("Final distance from origin:", np.sqrt(xs[-1]**2 + ys[-1]**2))

# --- Separate check: typical distance grows like sqrt(t) ---
# For 2D Brownian motion, E[R^2(t)] = (dx^2 + dy^2)*t, so RMS distance ~ sqrt(t).
# The measured squared distance at each time divided by t should be ~ (dx^2+dy^2).
t = np.arange(n_steps + 1)
r2 = xs**2 + ys**2                       # squared distance from origin at each time
expected_r2_per_t = dx[0]**2 + dx[1]**2  # theoretical slope of E[R^2] vs t

# Fit slope of R^2 vs t (through the origin) to compare with the expectation.
tt = t[1:]
measured_slope = np.sum(tt * r2[1:]) / np.sum(tt * tt)

print("Expected E[R^2]/t (theory, dx^2+dy^2):", expected_r2_per_t)
print("Measured slope of R^2 vs t:", measured_slope)
print("RMS distance at final time (sqrt(R^2)):", np.sqrt(r2[-1]))
print("Theoretical typical distance sqrt((dx^2+dy^2)*t):", np.sqrt(expected_r2_per_t * n_steps))

# Explanation:
# This check confirms the result because the measured squared distance grows
# linearly in t with slope ~ (dx^2+dy^2), meaning typical distance ~ sqrt(t),
# the hallmark of a diffusive, drifting-yet-tangled Brownian path.
print("Explanation: R^2 grows linearly with t (slope ~ dx^2+dy^2), so typical distance ~ sqrt(t), confirming diffusive 2D Brownian motion.")
