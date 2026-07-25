import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
D = 0.5            # diffusion coefficient
dt = 10.0          # deliberately large time step
t_max = 1000.0     # total simulated time
n_walks = 1000     # number of independent walkers
np.random.seed(12) # reproducibility

# Number of steps and the per-step standard deviation.
# The key scaling: step_sd = sqrt(2*D*dt). This keeps the variance
# growth correct (2*D*t) no matter how large dt is chosen.
n_steps = int(round(t_max / dt))
step_sd = np.sqrt(2.0 * D * dt)
print(f"Number of steps: {n_steps}")
print(f"Per-step standard deviation step_sd = sqrt(2*D*dt) = {step_sd}")

# --- Explicit Wiener-process simulation via Gaussian increments ---
# positions[i] holds the current position of every walker at step i.
# We build the walk one step at a time (no one-shot cumsum routine),
# so the mechanics are visible.
times = np.arange(0, n_steps + 1) * dt          # time at each recorded step
positions = np.zeros((n_steps + 1, n_walks))    # rows = time, cols = walkers

x = np.zeros(n_walks)  # all walkers start at the origin
for i in range(1, n_steps + 1):
    # Draw one Gaussian increment N(0, step_sd) per walker...
    increments = np.random.normal(loc=0.0, scale=step_sd, size=n_walks)
    # ...and add it to the running position: x_next = x + N(0, step_sd).
    x = x + increments
    positions[i] = x

# --- Variance check: compare empirical variance to the theoretical 2*D*t ---
print("\nVariance check (time : empirical_var : theoretical 2*D*t):")
for i in range(0, n_steps + 1):
    t = times[i]
    emp_var = np.var(positions[i], ddof=1) if i > 0 else 0.0
    theo_var = 2.0 * D * t
    print(f"t = {t:7.1f} : empirical = {emp_var:12.4f} : theoretical = {theo_var:12.4f}")

# Focused check at the final time.
final_emp_var = np.var(positions[-1], ddof=1)
final_theo_var = 2.0 * D * t_max
print(f"\nFinal empirical variance at t={t_max}: {final_emp_var}")
print(f"Final theoretical variance 2*D*t at t={t_max}: {final_theo_var}")
print(f"Ratio empirical/theoretical at final time: {final_emp_var / final_theo_var}")

# --- Box plot of walker positions over time ---
fig, ax = plt.subplots(figsize=(12, 6))
ax.boxplot(positions[1:].T, positions=times[1:], widths=dt * 0.6,
           manage_ticks=False, showfliers=True)
# Overlay the theoretical +/- 1 sd envelope: sd = sqrt(2*D*t).
ax.plot(times, np.sqrt(2.0 * D * times), 'r--', label='+/- sqrt(2*D*t)')
ax.plot(times, -np.sqrt(2.0 * D * times), 'r--')
ax.set_xlabel("time t")
ax.set_ylabel("walker position x")
ax.set_title("Wiener process: position spread grows as variance = 2*D*t (dt=10)")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.5.1_s2.png")

# Explanation: The check confirms the result because the empirical variance of
# the 1000 walkers matches 2*D*t at every recorded time even with the large
# dt=10, showing the sqrt(2*D*dt) increment scaling makes the simulated
# variance step-size independent and equal to the true diffusion law.
print("\nExplanation: The empirical variance tracking 2*D*t at each time confirms "
      "that the sqrt(2*D*dt) increment scaling keeps the variance correct and "
      "independent of the (large) step size dt.")
