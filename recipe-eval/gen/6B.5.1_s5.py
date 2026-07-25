import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters ----
D = 0.5           # diffusion coefficient
dt = 10.0         # deliberately large time step
t_final = 1000.0  # total simulated time
n_walks = 1000    # number of independent walkers
seed = 12

# Standard deviation of each Gaussian increment: sqrt(2*D*dt).
# This scaling is what keeps the variance correct regardless of dt,
# because per step Var = (sqrt(2*D*dt))**2 = 2*D*dt, and variances add.
step_sd = np.sqrt(2.0 * D * dt)
print(f"Step standard deviation step_sd = sqrt(2*D*dt) = {step_sd}")

# Number of steps to reach t_final with step size dt.
n_steps = int(round(t_final / dt))
print(f"Number of steps = {n_steps}")

rng = np.random.default_rng(seed)

# Time grid (t=0 plus one point after each step).
times = np.arange(n_steps + 1) * dt

# Positions array: rows = walkers, columns = time points.
positions = np.zeros((n_walks, n_steps + 1))

# ---- Explicit Wiener-process simulation via Gaussian increments ----
# Start all walkers at x = 0 (already set), then step forward one dt at a time.
for k in range(n_steps):
    x = positions[:, k]                       # current positions of all walkers
    increment = rng.normal(0.0, step_sd, n_walks)  # N(0, step_sd) increment
    positions[:, k + 1] = x + increment       # x_next = x + N(0, step_sd)

# ---- Variance check: empirical vs theoretical 2*D*t ----
emp_var = positions.var(axis=0, ddof=1)       # sample variance across walkers at each time
theo_var = 2.0 * D * times                     # theoretical variance 2*D*t

print("t\tempirical_var\ttheoretical_var(2*D*t)")
for t, ev, tv in zip(times, emp_var, theo_var):
    print(f"{t:.0f}\t{ev:.4f}\t{tv:.4f}")

# Focus on the final time as a summary check.
print(f"Final time t = {times[-1]:.0f}")
print(f"Empirical variance at final time = {emp_var[-1]:.4f}")
print(f"Theoretical variance 2*D*t at final time = {theo_var[-1]:.4f}")
print(f"Ratio empirical/theoretical at final time = {emp_var[-1] / theo_var[-1]:.4f}")

# ---- Box plot of walker positions over time ----
# To keep the plot readable, show a subset of time points.
col_idx = np.linspace(0, n_steps, 11).astype(int)
box_data = [positions[:, c] for c in col_idx]
box_labels = [f"{int(times[c])}" for c in col_idx]

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(box_data, labels=box_labels, showfliers=True)
# Overlay theoretical +/- 1 std envelope (sqrt(2*D*t)) to show variance growth.
env = np.sqrt(theo_var[col_idx])
ax.plot(range(1, len(col_idx) + 1), env, "r--", label="+1 std = sqrt(2*D*t)")
ax.plot(range(1, len(col_idx) + 1), -env, "r--")
ax.set_xlabel("time t")
ax.set_ylabel("walker position x")
ax.set_title("Wiener process: positions over time (D=0.5, dt=10), variance ~ 2*D*t")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.5.1_s5.png")

# ---- One-sentence explanation ----
# The check confirms the result because the empirical spread of the 1000 walkers
# matches the theoretical 2*D*t at every time even with the large dt=10, showing
# that scaling each Gaussian increment by sqrt(2*D*dt) makes the accumulated
# variance depend only on total elapsed time t and not on the step size.
print("Explanation: The empirical variance tracks 2*D*t at all times despite dt=10, "
      "confirming that the sqrt(2*D*dt) increment scaling makes total variance depend "
      "only on elapsed time t, not on step size.")
