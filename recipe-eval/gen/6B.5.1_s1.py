import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters ----
D = 0.5            # diffusion coefficient
dt = 10.0          # deliberately large time step
t_final = 1000.0   # total time
n_walks = 1000     # number of independent walkers
seed = 12

# ---- Wiener-process scaling ----
# The key idea: increment standard deviation scales as sqrt(2*D*dt),
# so variance per step is 2*D*dt and total variance after time t is 2*D*t
# regardless of how large dt is.
step_sd = np.sqrt(2 * D * dt)
print(f"Step standard deviation step_sd = sqrt(2*D*dt) = {step_sd}")

n_steps = int(round(t_final / dt))  # number of steps to reach t_final
print(f"Number of steps = {n_steps}")

# Time axis (include t = 0 at the start)
times = np.arange(0, n_steps + 1) * dt

# ---- Explicit simulation of the Gaussian random walk ----
rng = np.random.default_rng(seed)

# positions[step, walk]: store every walker's position at every time
positions = np.zeros((n_steps + 1, n_walks))

# Walk forward one step at a time, explicitly applying x_next = x + N(0, step_sd)
for step in range(1, n_steps + 1):
    # draw one Gaussian increment per walker, scaled by step_sd
    increments = rng.normal(loc=0.0, scale=step_sd, size=n_walks)
    # add increment to previous position (accumulate the random walk)
    positions[step] = positions[step - 1] + increments

# ---- Variance check: empirical vs theoretical 2*D*t ----
empirical_var = positions.var(axis=1, ddof=1)  # sample variance across walkers at each time
theoretical_var = 2 * D * times                # expected variance 2*D*t

print("Variance check (time, empirical_var, theoretical 2*D*t):")
for t, ev, tv in zip(times, empirical_var, theoretical_var):
    print(f"  t = {t:7.1f}  empirical = {ev:12.4f}  theoretical = {tv:12.4f}")

# Summarize the agreement at the final time
print(f"Final time t = {times[-1]}")
print(f"Empirical variance at final time = {empirical_var[-1]}")
print(f"Theoretical variance 2*D*t at final time = {theoretical_var[-1]}")
print(f"Ratio empirical/theoretical at final time = {empirical_var[-1] / theoretical_var[-1]}")

# ---- Box plot of walker positions over time ----
fig, ax = plt.subplots(figsize=(12, 6))
# box plot the distribution of positions at each time step
ax.boxplot(positions.T, positions=times, widths=dt * 0.6,
           manage_ticks=False, showfliers=True)
# overlay theoretical +/- sqrt(2*D*t) spread to show growing variance
ax.plot(times, np.sqrt(theoretical_var), 'r-', lw=2, label=r'$+\sqrt{2Dt}$')
ax.plot(times, -np.sqrt(theoretical_var), 'r-', lw=2, label=r'$-\sqrt{2Dt}$')
ax.set_xlabel("time t")
ax.set_ylabel("walker position x")
ax.set_title("Wiener process: position distribution over time (variance grows as 2*D*t)")
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.5.1_s1.png")

# ---- One-sentence explanation ----
# The check confirms the result because the empirical variance across walkers
# tracks the straight line 2*D*t even with the large dt=10, showing that the
# sqrt(2*D*dt) increment scaling makes the total variance independent of step size.
print("Explanation: The empirical variance matches 2*D*t even at dt=10, which confirms that scaling increments by sqrt(2*D*dt) keeps the accumulated variance correct regardless of step size.")
