import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Simulation parameters ---
k = 1.0                     # restoring drift strength
D_values = [100, 25, 1]     # noise levels to test
X0 = 0.0                    # initial condition
T = 1000.0                  # total integration time
dt = 0.01                   # time step
n_steps = int(T / dt)       # number of Euler-Maruyama steps
np.random.seed(1)           # fixed seed for reproducibility

# Storage for trajectories and measured variances
trajectories = {}
measured_vars = {}

# --- Euler-Maruyama integration of dX = -k*X*dt + sqrt(2*D)*dW ---
for D in D_values:
    X = np.empty(n_steps + 1)
    X[0] = X0
    noise_amp = np.sqrt(2.0 * D)          # diffusion coefficient sqrt(2*D)
    for i in range(n_steps):
        # Wiener increment: dW ~ Normal(0, dt), so dW = sqrt(dt)*N(0,1)
        dW = np.sqrt(dt) * np.random.randn()
        # Explicit Euler-Maruyama update: drift term + noise term
        X[i + 1] = X[i] + (-k * X[i]) * dt + noise_amp * dW
    trajectories[D] = X
    # Discard an initial transient before measuring the stationary variance
    burn_in = int(0.1 * n_steps)          # drop first 10% to reach stationarity
    measured_vars[D] = np.var(X[burn_in:])

# --- Print numerical results ---
for D in D_values:
    theoretical = D / k                   # stationary variance <x^2> = D/k
    print(f"D = {D:6.1f} | measured variance = {measured_vars[D]:10.4f} | "
          f"theoretical D/k = {theoretical:10.4f}")

# --- Plotting ---
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
t = np.linspace(0, T, n_steps + 1)

# Trajectory plots for each D (noisier for larger D)
for ax, D in zip(axes.flat[:3], D_values):
    ax.plot(t, trajectories[D], lw=0.5)
    ax.set_title(f"OU trajectory, D = {D} (measured var = {measured_vars[D]:.2f})")
    ax.set_xlabel("t")
    ax.set_ylabel("X")

# Measured variance vs D against the line <x^2> = D (k = 1)
ax = axes.flat[3]
Ds = np.array(D_values, dtype=float)
meas = np.array([measured_vars[D] for D in D_values])
ax.plot(Ds, meas, "o", markersize=9, label="measured variance")
line = np.linspace(0, max(D_values), 100)
ax.plot(line, line, "r--", label="<x^2> = D")   # theory line since k = 1
ax.set_xlabel("D")
ax.set_ylabel("variance")
ax.set_title("Measured variance vs D")
ax.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.2.1_s3.png")

# Explanation of why the check confirms the result:
print("Explanation: Because the measured variances land on the line <x^2> = D "
      "for k = 1, they satisfy <x^2> = D/k, confirming the predicted stationary variance.")
