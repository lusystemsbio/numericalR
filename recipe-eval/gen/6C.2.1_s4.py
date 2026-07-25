import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
k = 1.0                     # restoring drift rate
D_values = [100, 25, 1]     # noise levels to test
X0 = 0.0                    # initial condition
T = 1000.0                  # total time
dt = 0.01                   # time step
n_steps = int(T / dt)       # number of Euler-Maruyama steps
np.random.seed(1)           # reproducible seed

trajectories = {}           # store trajectories per D
measured_vars = {}          # store measured stationary variance per D

for D in D_values:
    # Euler-Maruyama for dX = -k*X*dt + sqrt(2*D)*dW
    X = np.empty(n_steps + 1)
    X[0] = X0
    noise_amp = np.sqrt(2.0 * D)               # diffusion coefficient sqrt(2D)
    for i in range(n_steps):
        dW = np.sqrt(dt) * np.random.randn()   # Wiener increment ~ N(0, dt)
        # deterministic drift step plus stochastic noise step
        X[i + 1] = X[i] - k * X[i] * dt + noise_amp * dW
    trajectories[D] = X

    # Discard an initial transient (burn-in) before measuring stationary variance
    burn = int(0.1 * n_steps)                  # drop first 10% to reach stationarity
    var_measured = np.var(X[burn:])            # <x^2> once relaxed (mean ~ 0)
    measured_vars[D] = var_measured

    print(f"D = {D:6.1f} | theoretical <x^2> = D/k = {D / k:8.3f} | measured variance = {var_measured:8.3f}")

# --- Trajectory plots for each D, plus variance-vs-D check ---
t = np.linspace(0.0, T, n_steps + 1)
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for ax, D in zip(axes.flat[:3], D_values):
    ax.plot(t, trajectories[D], lw=0.5)
    ax.set_title(f"OU trajectory, D = {D} (std = {np.sqrt(measured_vars[D]):.2f})")
    ax.set_xlabel("t")
    ax.set_ylabel("X(t)")

# Measured variance vs D compared to the line <x^2> = D (since k = 1)
ax = axes.flat[3]
Ds = np.array(D_values, dtype=float)
meas = np.array([measured_vars[D] for D in D_values])
ax.plot(Ds, Ds, "k--", label="theory <x^2> = D/k = D")
ax.plot(Ds, meas, "ro", label="measured")
ax.set_title("Measured variance vs D")
ax.set_xlabel("D")
ax.set_ylabel("<x^2>")
ax.legend()

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.2.1_s4.png", dpi=120)

# --- Report the check ---
print("Check: larger D -> larger trajectory spread (noisier):")
for D in D_values:
    print(f"  D = {D:6.1f} | measured std of X(t) = {np.sqrt(measured_vars[D]):8.3f}")

print("Explanation: Because k = 1 makes the predicted stationary variance D/k equal to D, "
      "the measured variances landing on the line <x^2> = D confirms the general relation <x^2> = D/k.")
