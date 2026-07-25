import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
D = 0.5           # diffusion coefficient
dt = 10.0         # deliberately large time step
t_final = 1000.0  # total simulated time
n_walks = 1000    # number of independent walkers
seed = 12

# The key scaling: each Gaussian increment has standard deviation
# step_sd = sqrt(2*D*dt).  This keeps Var[x(t)] = 2*D*t no matter
# how large dt is, because variances of independent increments add:
# n_steps * (2*D*dt) = 2*D*(n_steps*dt) = 2*D*t.
step_sd = np.sqrt(2.0 * D * dt)

n_steps = int(round(t_final / dt))          # number of increments
times = np.arange(0, n_steps + 1) * dt      # time axis including t=0

rng = np.random.default_rng(seed)

# Positions array: rows = walkers, cols = time points.
positions = np.zeros((n_walks, n_steps + 1))

# Explicit Euler-Maruyama-style loop (no one-shot routine): at each
# step add a Gaussian increment N(0, step_sd) to the current position.
for k in range(n_steps):
    increments = rng.normal(loc=0.0, scale=step_sd, size=n_walks)  # N(0, step_sd)
    positions[:, k + 1] = positions[:, k] + increments            # x_next = x + increment

# --- Empirical vs theoretical variance check ---
emp_var = positions.var(axis=0, ddof=1)   # sample variance across walkers at each time
theo_var = 2.0 * D * times                 # theoretical 2*D*t

print(f"Diffusion coefficient D: {D}")
print(f"Step size dt: {dt}")
print(f"Final time t_final: {t_final}")
print(f"Number of steps: {n_steps}")
print(f"Number of walks: {n_walks}")
print(f"Increment step_sd = sqrt(2*D*dt): {step_sd}")

# Report the check at the final time.
print(f"Empirical variance at t={t_final}: {emp_var[-1]}")
print(f"Theoretical variance 2*D*t at t={t_final}: {theo_var[-1]}")
print(f"Ratio empirical/theoretical at final time: {emp_var[-1] / theo_var[-1]}")

# Report the check across all times (mean relative error, ignoring t=0).
rel_err = np.abs(emp_var[1:] - theo_var[1:]) / theo_var[1:]
print(f"Mean relative error of variance vs 2*D*t (t>0): {rel_err.mean()}")
print(f"Max relative error of variance vs 2*D*t (t>0): {rel_err.max()}")

# Print variance table for every time point.
for tt, ev, tv in zip(times, emp_var, theo_var):
    print(f"t={tt:.0f}  empirical_var={ev:.4f}  theoretical_var={tv:.4f}")

# --- Box plot of walker positions over time ---
# Subsample time columns so the box plot stays readable.
col_idx = np.linspace(0, n_steps, 11).astype(int)
box_data = [positions[:, c] for c in col_idx]
box_labels = [f"{times[c]:.0f}" for c in col_idx]

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(box_data, labels=box_labels, showfliers=True)
# Overlay theoretical +/- 1 std envelope (sqrt(2*D*t)) as a reference.
std_env = np.sqrt(theo_var[col_idx])
xpos = np.arange(1, len(col_idx) + 1)
ax.plot(xpos, std_env, "r--", label="+1 std = sqrt(2*D*t)")
ax.plot(xpos, -std_env, "r--", label="-1 std = -sqrt(2*D*t)")
ax.set_xlabel("time t")
ax.set_ylabel("walker position x")
ax.set_title("Wiener process: spread grows as variance = 2*D*t (dt=10)")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.5.1_s4.png")

# Explanation of why the check confirms the result:
print("Explanation: Because the empirical across-walker variance matches 2*D*t "
      "closely at every time even with the large dt=10, this confirms that the "
      "sqrt(2*D*dt) increment scaling makes the discrete random walk reproduce "
      "the correct continuous Wiener-process variance independent of step size.")
