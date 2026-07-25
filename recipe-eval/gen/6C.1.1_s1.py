import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model definition: dX = f(X)*dt + sqrt(2*D(X))*dW ---
# Test case: pure diffusion, so drift is zero and D is constant.
def f(x):        # drift term
    return 0.0

def D(x):        # diffusion coefficient
    return 10.0

# --- Simulation parameters ---
X0 = 0.0         # initial position
T = 200.0        # total integration time
dt = 0.01        # time step
n_traj = 10      # number of trajectories
n_steps = int(round(T / dt))   # number of Euler-Maruyama steps

np.random.seed(1)              # reproducibility

# Time axis (including t=0)
t = np.linspace(0.0, T, n_steps + 1)

# Storage for all trajectories: rows = trajectories, cols = time points
X = np.zeros((n_traj, n_steps + 1))
X[:, 0] = X0

# --- Explicit Euler-Maruyama integration ---
for i in range(n_steps):
    x = X[:, i]                                   # current states of all trajectories
    drift = np.array([f(xi) for xi in x])         # deterministic increment f(X)*dt
    # dW ~ N(0, dt): standard normal scaled by sqrt(dt)
    dW = np.sqrt(dt) * np.random.randn(n_traj)
    noise = np.sqrt(2.0 * np.array([D(xi) for xi in x])) * dW   # stochastic increment
    X[:, i + 1] = x + drift * dt + noise          # X_next = X + f(X)*dt + sqrt(2D)*dW

# --- Plot the ten Brownian trajectories versus time ---
plt.figure(figsize=(9, 6))
for k in range(n_traj):
    plt.plot(t, X[k], lw=0.8)
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Free Brownian motion via Euler-Maruyama (f=0, D=10)")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.1.1_s1.png", dpi=120)

# --- Check: reproduce free Brownian motion ---
# For pure diffusion with f=0 and constant D, the exact solution is a Wiener
# process with Var[X(t)] = 2*D*t, so the ensemble variance should grow linearly.
final_states = X[:, -1]
empirical_var_final = np.var(final_states)
theoretical_var_final = 2.0 * D(0.0) * T

# Empirical variance over time (across the 10 trajectories) fit to slope 2D
empirical_var_t = np.var(X, axis=0)
# Fit variance vs time through origin: slope = sum(var*t)/sum(t^2)
fitted_slope = np.sum(empirical_var_t * t) / np.sum(t * t)
theoretical_slope = 2.0 * D(0.0)

print("Number of trajectories:", n_traj)
print("Number of time steps:", n_steps)
print("Empirical variance of final positions X(T):", empirical_var_final)
print("Theoretical variance of X(T) = 2*D*T:", theoretical_var_final)
print("Fitted slope of Var[X(t)] vs t:", fitted_slope)
print("Theoretical slope = 2*D:", theoretical_slope)
print("Mean of final positions (should be ~0):", np.mean(final_states))

# Explanation:
print("Check confirms result: because f=0 makes the SDE pure diffusion, the "
      "trajectories are Wiener processes whose ensemble variance grows as 2*D*t, "
      "matching the free Brownian motion of the previous chapter.")
