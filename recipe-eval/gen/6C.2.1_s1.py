import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
k = 1.0                      # restoring drift strength
D_values = [100, 25, 1]      # noise levels to test
X0 = 0.0                     # initial condition
T = 1000.0                   # total simulation time
dt = 0.01                    # time step
N = int(T / dt)              # number of steps
seed = 1

t = np.linspace(0.0, T, N + 1)  # time grid

# storage for measured stationary variances
measured_vars = []
trajectories = {}

for D in D_values:
    # reset the RNG for each D so runs are reproducible and comparable
    rng = np.random.default_rng(seed)

    X = np.empty(N + 1)
    X[0] = X0

    # noise amplitude for dW increments: sqrt(2*D)*dW, with dW ~ N(0, dt)
    noise_coeff = np.sqrt(2.0 * D)

    # --- Euler-Maruyama integration, done explicitly step by step ---
    for i in range(N):
        drift = -k * X[i]                       # deterministic restoring term
        dW = rng.normal(0.0, np.sqrt(dt))       # Wiener increment, variance dt
        # X_{n+1} = X_n + drift*dt + sqrt(2D)*dW
        X[i + 1] = X[i] + drift * dt + noise_coeff * dW

    trajectories[D] = X

    # measure stationary variance, discarding an initial transient (first 10%)
    burn_in = N // 10
    var_measured = np.var(X[burn_in:])
    measured_vars.append(var_measured)

    print(f"D = {D:6.1f}  theoretical <x^2> = D/k = {D / k:8.3f}  measured <x^2> = {var_measured:8.3f}")

measured_vars = np.array(measured_vars)

# --- Plotting ---
fig, axes = plt.subplots(1, 4, figsize=(20, 5))

# trajectory plots for each D
for ax, D in zip(axes[:3], D_values):
    ax.plot(t, trajectories[D], lw=0.5)
    ax.set_title(f"OU trajectory, D = {D}")
    ax.set_xlabel("t")
    ax.set_ylabel("X(t)")

# measured variance vs D against the line <x^2> = D (since k = 1)
D_arr = np.array(D_values, dtype=float)
line_D = np.linspace(0, max(D_values) * 1.05, 100)
axes[3].plot(line_D, line_D, "k--", label="<x^2> = D  (D/k, k=1)")
axes[3].plot(D_arr, measured_vars, "ro", markersize=8, label="measured variance")
axes[3].set_title("Measured variance vs D")
axes[3].set_xlabel("D")
axes[3].set_ylabel("<x^2>")
axes[3].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.2.1_s1.png")

# --- Report the check ---
for D, v in zip(D_values, measured_vars):
    print(f"check: D = {D:6.1f}  measured/theoretical ratio = {v / (D / k):6.4f}")

# Explanation: Because larger D produces visibly noisier trajectories whose measured
# variances lie on the line <x^2> = D (with k = 1), the simulation reproduces the
# analytic stationary result <x^2> = D/k.
print("Explanation: the measured variances landing on the line <x^2> = D confirms <x^2> = D/k since k = 1.")
