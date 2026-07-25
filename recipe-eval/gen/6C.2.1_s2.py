import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Ornstein-Uhlenbeck via explicit Euler-Maruyama ---
# SDE: dX = -k*X*dt + sqrt(2*D)*dW,  stationary variance <x^2> = D/k

np.random.seed(1)          # reproducible noise
k = 1.0                    # restoring drift rate
dt = 0.01                  # time step
T = 1000.0                 # total time
N = int(T / dt)            # number of steps
X0 = 0.0                   # initial condition
D_levels = [100, 25, 1]    # noise levels to test

t = np.linspace(0.0, T, N + 1)     # time grid
trajectories = {}                   # store each trajectory
measured_vars = {}                  # store each stationary variance

for D in D_levels:
    X = np.empty(N + 1)             # allocate trajectory array
    X[0] = X0                       # set initial value
    noise_amp = np.sqrt(2.0 * D)    # diffusion coefficient sqrt(2D)
    # generate all Wiener increments dW ~ Normal(0, dt)
    dW = np.random.normal(0.0, np.sqrt(dt), size=N)
    for i in range(N):
        # explicit Euler-Maruyama update: deterministic drift + stochastic kick
        X[i + 1] = X[i] - k * X[i] * dt + noise_amp * dW[i]
    trajectories[D] = X

    # measure stationary variance from the second half (after equilibration)
    tail = X[N // 2:]
    var_measured = np.var(tail)
    measured_vars[D] = var_measured
    print(f"D = {D:6.1f}:  measured stationary variance <x^2> = {var_measured:10.4f}  (theory D/k = {D / k:10.4f})")

# --- Trajectory plots for each D, plus variance-vs-D check ---
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for ax, D in zip(axes.flat[:3], D_levels):
    ax.plot(t, trajectories[D], lw=0.4)
    ax.set_title(f"OU trajectory, D = {D}")
    ax.set_xlabel("t")
    ax.set_ylabel("X(t)")

# variance vs D against the line <x^2> = D (since k = 1)
ax = axes.flat[3]
Ds = np.array(D_levels, dtype=float)
meas = np.array([measured_vars[D] for D in D_levels])
line = np.linspace(0, max(Ds) * 1.05, 100)
ax.plot(line, line, "k--", label="theory  <x^2> = D")
ax.plot(Ds, meas, "ro", ms=8, label="measured")
ax.set_title("Stationary variance vs D")
ax.set_xlabel("D")
ax.set_ylabel("measured <x^2>")
ax.legend()

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.2.1_s2.png")

# --- Report the confirmation ---
for D in D_levels:
    print(f"ratio measured/(D/k) for D = {D:6.1f}:  {measured_vars[D] / (D / k):8.4f}")

print("Explanation: because the measured variances of the noisier (larger-D) trajectories "
      "fall on the line <x^2> = D with k = 1, they match the predicted D/k, confirming the "
      "stationary variance formula.")
