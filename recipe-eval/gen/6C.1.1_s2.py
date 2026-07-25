import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------
# General SDE:  dX = f(X)*dt + sqrt(2*D(X))*dW,  dW ~ N(0, dt)
# Solved explicitly with the Euler-Maruyama scheme.
# Test case: free Brownian motion -> f(X)=0, D(X)=D constant.
# -----------------------------------------------------------

# Model definitions for the pure-diffusion test
def f(x):
    return 0.0          # no drift -> free diffusion

D = 10.0                # constant diffusion coefficient
def Dfun(x):
    return D

# Simulation parameters
X0 = 0.0
t_end = 200.0
dt = 0.01
n_steps = int(round(t_end / dt))
n_traj = 10

np.random.seed(1)       # reproducibility

# Time grid
t = np.linspace(0.0, t_end, n_steps + 1)

# Storage: rows = trajectories, cols = time points
X = np.zeros((n_traj, n_steps + 1))
X[:, 0] = X0

# Explicit Euler-Maruyama integration
for k in range(n_steps):
    x = X[:, k]
    # Wiener increments dW ~ N(0, dt) for every trajectory
    dW = np.random.normal(0.0, np.sqrt(dt), size=n_traj)
    # Update:  X_next = X + f(X)*dt + sqrt(2 D(X)) * dW
    X[:, k + 1] = x + f(x) * dt + np.sqrt(2.0 * Dfun(x)) * dW

# -----------------------------------------------------------
# Check against theory of free Brownian motion.
# For dX = sqrt(2D) dW the position at time t is Gaussian with
# mean 0 and variance <X^2> = 2 D t  (i.e. MSD = 2 D t).
# We compare the empirical ensemble variance at t_end with 2 D t_end.
# -----------------------------------------------------------
emp_mean_end = np.mean(X[:, -1])
emp_var_end = np.var(X[:, -1])
theory_var_end = 2.0 * D * t_end

print(f"Number of trajectories: {n_traj}")
print(f"Number of Euler-Maruyama steps per trajectory: {n_steps}")
print(f"Empirical mean of X at t={t_end}: {emp_mean_end}")
print(f"Empirical variance of X at t={t_end}: {emp_var_end}")
print(f"Theoretical variance 2*D*t at t={t_end}: {theory_var_end}")
print(f"Final X value of each trajectory: {X[:, -1]}")

# One-sentence explanation:
# The check confirms the result because the Euler-Maruyama trajectories,
# having zero drift and increments of variance 2*D*dt, accumulate to a
# Gaussian spread whose variance grows as 2*D*t -- exactly the free
# Brownian motion (mean-squared displacement = 2*D*t) of the previous chapter.
print("Why the check works: with f=0 the scheme sums independent N(0, 2*D*dt) "
      "increments, so X(t) is Gaussian with variance 2*D*t, which is precisely "
      "the free Brownian motion of the previous chapter.")

# Plot the ten Brownian trajectories versus time
plt.figure(figsize=(9, 6))
for i in range(n_traj):
    plt.plot(t, X[i, :], lw=0.8)
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Free Brownian motion via Euler-Maruyama (f=0, D=10)")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.1.1_s2.png")
