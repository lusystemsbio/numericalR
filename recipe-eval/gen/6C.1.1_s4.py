import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model definition: dX = f(X) dt + sqrt(2 D(X)) dW ----
# Test case: pure (free) diffusion => drift f = 0, constant diffusion D = 10.
def f(x):
    return 0.0          # no drift

def D(x):
    return 10.0         # constant diffusion coefficient

# ---- Integration parameters ----
X0 = 0.0
t_end = 200.0
dt = 0.01
n_traj = 10
Dconst = 10.0

np.random.seed(1)                       # reproducibility (seed 1)
n_steps = int(round(t_end / dt))        # number of Euler-Maruyama steps
t = np.linspace(0.0, t_end, n_steps + 1)  # time grid

# Storage for all trajectories: rows = time, cols = trajectory
X = np.zeros((n_steps + 1, n_traj))
X[0, :] = X0

# ---- Euler-Maruyama integration (implemented explicitly) ----
for i in range(n_steps):
    x = X[i, :]                                  # current states
    drift = f(x)                                 # deterministic part: f(X)*dt
    # Wiener increments dW ~ N(0, dt), one per trajectory
    dW = np.random.normal(0.0, np.sqrt(dt), size=n_traj)
    diffusion = np.sqrt(2.0 * D(x))              # noise amplitude sqrt(2 D(X))
    X[i + 1, :] = x + drift * dt + diffusion * dW  # X_next = X + f dt + sqrt(2D) dW

# ---- Plot the ten Brownian trajectories vs time ----
plt.figure(figsize=(9, 5))
for j in range(n_traj):
    plt.plot(t, X[:, j], lw=0.8)
plt.xlabel("t")
plt.ylabel("X(t)")
plt.title("Ten free Brownian trajectories (Euler-Maruyama, f=0, D=10)")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.1.1_s4.png")

# ---- Check: reproduce free Brownian motion of the previous chapter ----
# For free BM, X(t) ~ N(0, 2 D t): mean 0 and variance growing linearly as 2 D t.
final = X[-1, :]                        # X at t = t_end across trajectories
emp_mean = np.mean(final)
emp_var = np.var(final)
theo_mean = 0.0
theo_var = 2.0 * Dconst * t_end        # variance of free BM at t_end

print(f"Number of steps: {n_steps}")
print(f"Final time: {t[-1]}")
print(f"Empirical mean of X at t={t_end}: {emp_mean}")
print(f"Theoretical mean of X at t={t_end}: {theo_mean}")
print(f"Empirical variance of X at t={t_end}: {emp_var}")
print(f"Theoretical variance (2*D*t) at t={t_end}: {theo_var}")
print("Explanation: With zero drift and constant D, Euler-Maruyama sums independent "
      "N(0,dt) increments scaled by sqrt(2D), so X(t) is Gaussian with mean 0 and "
      "variance 2*D*t — exactly the free Brownian motion of the previous chapter, "
      "confirming the integrator reproduces it.")
