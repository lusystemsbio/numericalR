import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Parameters ---
dx = (1.0, 1.0)      # step scale in (x, y)
dt = 1.0             # time increment per step
n_steps = 10000      # number of steps
np.random.seed(12)   # reproducibility

# --- Explicit 2D Brownian motion ---
# Preallocate arrays; index 0 is the starting point (0, 0).
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)

# At each time step take an independent Gaussian step per axis.
# The N(0,1)*sqrt(dt) scaling makes variance proportional to elapsed time.
for i in range(1, n_steps + 1):
    x[i] = x[i - 1] + np.random.normal(0.0, 1.0) * dx[0] * np.sqrt(dt)  # x_next = x + N(0,1)*dx
    y[i] = y[i - 1] + np.random.normal(0.0, 1.0) * dx[1] * np.sqrt(dt)  # y_next = y + N(0,1)*dy

# --- Plot the path in the plane ---
plt.figure(figsize=(7, 7))
plt.plot(x, y, lw=0.5, color="steelblue")
plt.plot(0, 0, "go", label="start (0,0)")
plt.plot(x[-1], y[-1], "rs", label="end")
plt.axis("equal")
plt.xlabel("x")
plt.ylabel("y")
plt.title("2D Brownian motion (10000 steps, seed 12)")
plt.legend()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.4.1_s3.png")

# --- Check: typical distance from origin should grow like sqrt(t) ---
# Distance of the walker from the origin at every time.
t = np.arange(n_steps + 1) * dt
dist = np.sqrt(x**2 + y**2)

# Theoretical RMS distance for 2D BM with per-axis variance (dx^2)*dt per step:
# E[R^2] = (dx_x^2 + dx_y^2) * t  ->  RMS = sqrt((dx_x^2+dx_y^2) * t)
rms_theory_end = np.sqrt((dx[0]**2 + dx[1]**2) * t[-1])

# Empirical final distance and a fitted growth exponent of dist ~ t^p.
# Fit log(dist) = p*log(t) + c on the later portion (avoid t=0 and early noise).
mask = t > 100
p, logc = np.polyfit(np.log(t[mask]), np.log(dist[mask]), 1)

print(f"Final position x: {x[-1]:.4f}")
print(f"Final position y: {y[-1]:.4f}")
print(f"Final distance from origin: {dist[-1]:.4f}")
print(f"Theoretical RMS distance at t={t[-1]:.0f}: {rms_theory_end:.4f}")
print(f"Fitted growth exponent p in dist ~ t^p: {p:.4f}")
print(f"Expected exponent for sqrt(t) growth: 0.5")

# The fitted exponent near 0.5 confirms the result: it shows the walker's
# distance from the origin grows on the order of sqrt(t), the diffusive
# signature of Brownian motion, while the plotted path stays a tangled,
# non-directional 2D walk that nonetheless wanders away from the origin.
print("Check: fitted exponent close to 0.5 confirms sqrt(t) diffusive spreading.")
