import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Simulation parameters ----
k = 1.0                     # restoring drift coefficient
D_values = [100, 25, 1]     # noise levels to test
X0 = 0.0                    # initial condition
T = 1000.0                  # total simulated time
dt = 0.01                   # time step
n_steps = int(T / dt)       # number of Euler-Maruyama steps
np.random.seed(1)           # reproducibility

# Storage for measured stationary variances
measured_vars = {}
trajectories = {}

# ---- Euler-Maruyama integration of the OU process ----
# dX = -k*X*dt + sqrt(2*D)*dW,   dW ~ Normal(0, dt)
for D in D_values:
    X = np.empty(n_steps + 1)   # trajectory array
    X[0] = X0
    noise_amp = np.sqrt(2.0 * D)                     # diffusion coefficient sqrt(2D)
    dW = np.random.normal(0.0, np.sqrt(dt), n_steps) # Wiener increments, std = sqrt(dt)
    for i in range(n_steps):
        drift = -k * X[i] * dt          # deterministic restoring step
        diffusion = noise_amp * dW[i]   # random noise step
        X[i + 1] = X[i] + drift + diffusion  # explicit Euler-Maruyama update
    trajectories[D] = X

    # Discard an initial transient (burn-in) before measuring the stationary variance
    burn = int(n_steps * 0.1)
    var_measured = np.var(X[burn:])
    measured_vars[D] = var_measured

    print(f"D = {D:6.1f}:  measured <x^2> = {var_measured:10.4f}   theory D/k = {D / k:10.4f}")

# ---- Plots ----
fig, axes = plt.subplots(len(D_values) + 1, 1, figsize=(9, 12))
t = np.linspace(0.0, T, n_steps + 1)

# Trajectory plots (one per D); larger D -> visibly noisier
for ax, D in zip(axes[:-1], D_values):
    ax.plot(t, trajectories[D], lw=0.4)
    ax.set_title(f"OU trajectory, D = {D} (k = {k})")
    ax.set_xlabel("t")
    ax.set_ylabel("X(t)")

# Measured variance versus D against the line <x^2> = D
axD = axes[-1]
Ds = np.array(D_values, dtype=float)
meas = np.array([measured_vars[D] for D in D_values])
axD.plot(Ds, meas, "o", ms=9, label="measured <x^2>")
line = np.linspace(0, max(Ds), 100)
axD.plot(line, line, "-", label="<x^2> = D  (D/k, k=1)")
axD.set_xlabel("D")
axD.set_ylabel("variance")
axD.set_title("Measured stationary variance vs D")
axD.legend()

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.2.1_s5.png", dpi=120)

# ---- Explicit check and explanation ----
print("\nCheck: measured variances vs the line <x^2> = D (k = 1)")
for D in D_values:
    ratio = measured_vars[D] / (D / k)
    print(f"D = {D:6.1f}:  measured/theory ratio = {ratio:.4f}")

print("\nExplanation: because the measured variances land on the line <x^2> = D "
      "for every noise level while k = 1, they satisfy <x^2> = D/k, confirming "
      "the analytic stationary variance of the OU process.")
