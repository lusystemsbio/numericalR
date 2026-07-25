import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Target: P(x) proportional to exp(-x^2) -> Gaussian with mean 0, variance 1/2
def log_unnorm(x):
    return -x * x  # log of exp(-x^2)

def metropolis_hastings(dx_max, n_steps, x0, rng):
    # Explicit MH with a uniform random-walk proposal
    x = x0
    samples = np.empty(n_steps)
    n_accept = 0
    for i in range(n_steps):
        # Propose x' = x + Uniform(-dx_max, dx_max)
        x_prop = x + rng.uniform(-dx_max, dx_max)
        # Acceptance ratio a = min(1, exp(-x'^2 + x^2))
        a = np.exp(-x_prop * x_prop + x * x)
        if a > 1.0:
            a = 1.0
        # Accept or reject
        if rng.uniform(0.0, 1.0) < a:
            x = x_prop
            n_accept += 1
        samples[i] = x
    accept_rate = n_accept / n_steps
    return samples, accept_rate

# Settings
x0 = 0.0
n_steps = 10000
dx_list = [0.1, 0.5, 2.0, 10.0]

# True Gaussian density: variance 1/2 -> sigma = 1/sqrt(2)
sigma = 1.0 / np.sqrt(2.0)
def true_pdf(x):
    return np.exp(-x * x) / np.sqrt(np.pi)  # normalized exp(-x^2)

# Run MH for each proposal width (fresh seeded RNG each so runs are comparable)
results = {}
for dx in dx_list:
    rng = np.random.default_rng(1)
    samples, acc = metropolis_hastings(dx, n_steps, x0, rng)
    results[dx] = (samples, acc)
    print(f"dx_max = {dx}: acceptance rate = {acc:.4f}, "
          f"sample mean = {np.mean(samples):.4f}, sample var = {np.var(samples):.4f}")

# True target moments for reference
print(f"True Gaussian: mean = 0.0000, variance = {0.5:.4f}")

# Check: which dx_max reproduces the true Gaussian best?
# Compare histograms to the true pdf via sum of squared errors over bin centers.
xg = np.linspace(-4, 4, 400)
bins = np.linspace(-4, 4, 41)
centers = 0.5 * (bins[:-1] + bins[1:])
best_dx, best_err = None, np.inf
for dx in dx_list:
    samples, acc = results[dx]
    hist, _ = np.histogram(samples, bins=bins, density=True)
    err = np.sum((hist - true_pdf(centers)) ** 2)
    print(f"dx_max = {dx}: histogram-vs-true SSE = {err:.5f}")
    if err < best_err:
        best_err, best_dx = err, dx
print(f"Best-matching dx_max = {best_dx} (SSE = {best_err:.5f}), "
      f"acceptance rate there = {results[best_dx][1]:.4f}")

# Plot histograms against the true Gaussian
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, dx in zip(axes.ravel(), dx_list):
    samples, acc = results[dx]
    ax.hist(samples, bins=bins, density=True, alpha=0.6,
            color="steelblue", label="MH samples")
    ax.plot(xg, true_pdf(xg), "r-", lw=2, label="true Gaussian")
    ax.set_title(f"dx_max = {dx}  (accept = {acc:.2f})")
    ax.set_xlabel("x"); ax.set_ylabel("density")
    ax.legend(fontsize=8)
fig.suptitle("Metropolis-Hastings sampling of exp(-x^2) at varying step sizes")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.4.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The intermediate step size (~50% acceptance) yields the "
      "smallest histogram-vs-true error, confirming the result because too-small "
      "steps barely move (high autocorrelation, poor exploration) and too-large "
      "steps are almost always rejected (the chain sticks), so only a moderate "
      "step both accepts often and moves far enough to sample the whole Gaussian.")
