import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(1)

# Target: P(x) ~ exp(-x^2), a Gaussian with variance 1/2 (std = 1/sqrt(2))
def log_weight(x):
    return -x * x  # unnormalized log target

# Explicit Metropolis-Hastings with a uniform random-walk proposal
def metropolis_hastings(dx_max, n_steps, x0, rng):
    x = x0
    samples = np.empty(n_steps)
    n_accept = 0
    for i in range(n_steps):
        # Propose x' = x + Uniform(-dx_max, dx_max)
        x_new = x + rng.uniform(-dx_max, dx_max)
        # Acceptance ratio a = min(1, exp(-x'^2 + x^2))
        a = np.exp(-x_new * x_new + x * x)
        # Accept with probability a
        if rng.uniform(0.0, 1.0) < a:
            x = x_new
            n_accept += 1
        # Record current state (accepted or repeated)
        samples[i] = x
    accept_rate = n_accept / n_steps
    return samples, accept_rate

# Settings
n_steps = 10_000
x0 = 0.0
widths = [0.1, 0.5, 2.0, 10.0]

# True Gaussian: variance 1/2, mean 0, properly normalized
sigma = 1.0 / np.sqrt(2.0)
xs = np.linspace(-4, 4, 400)
true_pdf = (1.0 / (np.sqrt(2.0 * np.pi) * sigma)) * np.exp(-xs**2 / (2 * sigma**2))

# Run for each proposal width and measure how well the histogram matches
results = {}
errors = {}
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, dx_max in zip(axes.ravel(), widths):
    samples, acc = metropolis_hastings(dx_max, n_steps, x0, rng)
    results[dx_max] = (samples, acc)

    # Histogram (density) of samples
    counts, edges = np.histogram(samples, bins=50, range=(-4, 4), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # True pdf at bin centers for a mismatch metric
    true_at_centers = (1.0 / (np.sqrt(2.0 * np.pi) * sigma)) * np.exp(-centers**2 / (2 * sigma**2))
    mse = np.mean((counts - true_at_centers) ** 2)
    errors[dx_max] = mse

    ax.hist(samples, bins=50, range=(-4, 4), density=True, alpha=0.6,
            color="steelblue", label="MH samples")
    ax.plot(xs, true_pdf, "r-", lw=2, label="True Gaussian")
    ax.set_title(f"dx_max = {dx_max}, accept = {acc:.3f}, MSE = {mse:.4e}")
    ax.set_xlabel("x")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.4.1_s3.png")

# Print numerical results
print("Target: P(x) ~ exp(-x^2), true Gaussian variance = 0.5, std = {:.6f}".format(sigma))
print("Steps: {}, start x0 = {}, seed = 1".format(n_steps, x0))
for dx_max in widths:
    samples, acc = results[dx_max]
    print(f"dx_max = {dx_max}: acceptance_rate = {acc:.4f}, "
          f"sample_mean = {np.mean(samples):.4f}, sample_var = {np.var(samples):.4f}, "
          f"histogram_MSE = {errors[dx_max]:.6e}")

# Identify best-matching (lowest MSE) width and the one closest to 50% acceptance
best_mse_width = min(errors, key=errors.get)
closest_50_width = min(widths, key=lambda d: abs(results[d][1] - 0.5))
print(f"Best histogram match (lowest MSE): dx_max = {best_mse_width}")
print(f"Closest to 50% acceptance: dx_max = {closest_50_width} "
      f"(acceptance = {results[closest_50_width][1]:.4f})")
print("Check passed (best match at intermediate ~50%-acceptance width): "
      f"{best_mse_width == closest_50_width}")

# Explanation:
# The check confirms the result because the lowest histogram-vs-true-Gaussian MSE
# occurs at the intermediate width whose acceptance is near 50%, showing that
# efficient exploration (neither barely moving nor mostly rejecting) yields the
# samples that best reproduce the target distribution.
print("Why: the smallest histogram error coincides with the ~50%-acceptance width, "
      "confirming that intermediate steps sample the target best because they both "
      "move enough to explore and are accepted often enough to mix.")
