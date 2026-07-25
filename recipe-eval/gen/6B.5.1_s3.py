import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
D = 0.5            # diffusion coefficient
dt = 10.0          # deliberately large time step
t_max = 1000.0     # total time
n_walks = 1000     # number of independent walkers
seed = 12

# Gaussian increment scaling: step_sd = sqrt(2*D*dt)
# This keeps Var(x(t)) = 2*D*t regardless of how large dt is,
# because each increment has variance 2*D*dt and increments add up.
step_sd = np.sqrt(2.0 * D * dt)

# Number of steps to reach t_max
n_steps = int(round(t_max / dt))

rng = np.random.default_rng(seed)

# Storage: positions[walk, step_index], step_index 0..n_steps (0 = start)
positions = np.zeros((n_walks, n_steps + 1))
times = np.arange(n_steps + 1) * dt

# --- Explicit Wiener process simulation via Gaussian increments ---
# Start all walkers at x = 0
x = np.zeros(n_walks)
for k in range(1, n_steps + 1):
    # Draw one Gaussian increment per walker: N(0, step_sd)
    increments = rng.normal(0.0, step_sd, size=n_walks)
    # Advance the walk: x_next = x + N(0, step_sd)
    x = x + increments
    # Record positions at this time
    positions[:, k] = x

# --- Empirical vs theoretical variance check ---
emp_var = positions.var(axis=0, ddof=1)   # sample variance at each time
theo_var = 2.0 * D * times                 # theoretical variance 2*D*t

print("D =", D)
print("dt =", dt)
print("t_max =", t_max)
print("n_walks =", n_walks)
print("n_steps =", n_steps)
print("step_sd = sqrt(2*D*dt) =", step_sd)

# Print variance comparison at a selection of times
print("\nVariance check (time : empirical : theoretical 2*D*t):")
for idx in np.linspace(0, n_steps, 6).astype(int):
    print(f"  t={times[idx]:8.1f} : empirical={emp_var[idx]:12.4f} : theoretical={theo_var[idx]:12.4f}")

# Final-time detailed check
print("\nAt final time t =", times[-1])
print("Empirical variance   =", emp_var[-1])
print("Theoretical variance =", theo_var[-1])
print("Ratio empirical/theoretical =", emp_var[-1] / theo_var[-1])

# Overall fit: slope of empirical variance vs t (should be ~2*D)
slope = np.polyfit(times[1:], emp_var[1:], 1)[0]
print("\nFitted slope of empirical variance vs t =", slope)
print("Expected slope 2*D =", 2.0 * D)

# --- Box plot of walker positions over time ---
# Subsample times so the box plot stays readable
plot_idx = np.linspace(0, n_steps, 11).astype(int)
box_data = [positions[:, i] for i in plot_idx]
box_labels = [f"{int(times[i])}" for i in plot_idx]

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(box_data, tick_labels=box_labels, showfliers=True)
# Overlay theoretical +/- 1 std (sqrt(2*D*t)) to show variance growth
std_theo = np.sqrt(2.0 * D * times[plot_idx])
xpos = np.arange(1, len(plot_idx) + 1)
ax.plot(xpos, std_theo, "r--", label=r"$+\sqrt{2Dt}$")
ax.plot(xpos, -std_theo, "r--", label=r"$-\sqrt{2Dt}$")
ax.set_xlabel("time t")
ax.set_ylabel("walker position x")
ax.set_title(f"Wiener process positions (D={D}, dt={dt}, {n_walks} walks)\nvariance grows as 2*D*t")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.5.1_s3.png")

# One-sentence explanation:
print("\nExplanation: The check confirms the result because the empirical variance")
print("of the 1000 walkers matches 2*D*t (slope ~2*D) even at dt=10, showing that")
print("scaling increments by sqrt(2*D*dt) makes the variance step-size independent.")
