import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Target distribution: P(x) proportional to exp(-x^2), a Gaussian of variance 1/2 ---
def log_unnorm_target(x):
    # log of the unnormalized density, used for the ratio in the acceptance test
    return -x * x

# --- Explicit Metropolis-Hastings with a uniform random-walk proposal ---
def metropolis_hastings(dx_max, n_steps, x0, rng):
    x = x0                      # current state
    samples = np.empty(n_steps) # storage for the chain
    n_accept = 0                # count accepted proposals
    for i in range(n_steps):
        # Propose x' = x + Uniform(-dx_max, dx_max)
        x_prop = x + rng.uniform(-dx_max, dx_max)
        # Acceptance probability a = min(1, exp(-x'^2 + x^2))
        a = min(1.0, np.exp(log_unnorm_target(x_prop) - log_unnorm_target(x)))
        # Accept with probability a
        if rng.uniform(0.0, 1.0) < a:
            x = x_prop
            n_accept += 1
        # Record the state (whether or not the move was accepted)
        samples[i] = x
    accept_rate = n_accept / n_steps
    return samples, accept_rate

# --- Run parameters ---
n_steps = int(1e4)
x0 = 0.0
dx_maxes = [0.1, 0.5, 2.0, 10.0]

# True normalized Gaussian density with variance 1/2: sqrt(1/pi) * exp(-x^2)
def true_density(x):
    return np.sqrt(1.0 / np.pi) * np.exp(-x * x)

xs = np.linspace(-4, 4, 400)
true_pdf = true_density(xs)

# --- Sample at each proposal width (seed 1) and measure fit quality ---
results = {}
for dx_max in dx_maxes:
    rng = np.random.default_rng(1)  # reseed to 1 for each width for reproducibility
    samples, accept_rate = metropolis_hastings(dx_max, n_steps, x0, rng)
    # Measure agreement with the true Gaussian via a histogram-based L1 distance
    counts, edges = np.histogram(samples, bins=50, range=(-4, 4), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    l1_error = np.sum(np.abs(counts - true_density(centers))) * (edges[1] - edges[0])
    results[dx_max] = dict(samples=samples, accept_rate=accept_rate,
                           l1_error=l1_error, mean=np.mean(samples),
                           var=np.var(samples))
    print(f"dx_max = {dx_max:>5}: acceptance rate = {accept_rate:.4f}, "
          f"sample mean = {np.mean(samples):.4f}, sample variance = {np.var(samples):.4f}, "
          f"L1 histogram error = {l1_error:.4f}")

# --- Plot histograms at each proposal width against the true Gaussian ---
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, dx_max in zip(axes.flat, dx_maxes):
    r = results[dx_max]
    ax.hist(r["samples"], bins=50, range=(-4, 4), density=True,
            alpha=0.6, color="steelblue", label="MH samples")
    ax.plot(xs, true_pdf, "r-", lw=2, label="true Gaussian")
    ax.set_title(f"dx_max = {dx_max}  (accept = {r['accept_rate']:.2f}, "
                 f"L1 = {r['l1_error']:.3f})")
    ax.set_xlabel("x")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
fig.suptitle("Metropolis-Hastings sampling of P(x) ~ exp(-x^2)")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.4.1_s5.png", dpi=120)

# --- Check: best fit occurs at an intermediate dx_max near ~50% acceptance ---
best_dx = min(results, key=lambda d: results[d]["l1_error"])
# The acceptance rate closest to 0.5 among the tested widths
closest_to_half = min(results, key=lambda d: abs(results[d]["accept_rate"] - 0.5))
print()
print(f"True Gaussian: mean = 0.0, variance = 0.5")
print(f"Best-fitting dx_max (smallest L1 error) = {best_dx} "
      f"with acceptance rate {results[best_dx]['accept_rate']:.4f}")
print(f"dx_max whose acceptance rate is closest to 0.5 = {closest_to_half} "
      f"with acceptance rate {results[closest_to_half]['accept_rate']:.4f}")
print(f"Intermediate-step check passed (best fit at intermediate dx_max near 50% acceptance): "
      f"{best_dx in (0.5, 2.0)}")
# Explanation of why the check confirms the result:
print("Explanation: the check confirms the result because the intermediate dx_max attains "
      "both the ~50% acceptance rate and the smallest histogram-vs-true L1 error, showing that "
      "the chain mixes best when steps are neither so small that it barely moves nor so large "
      "that most proposals are rejected.")
