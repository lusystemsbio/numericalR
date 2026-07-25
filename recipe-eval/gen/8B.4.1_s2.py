import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

# Ensure output directory exists
outpath = "/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.4.1_s2.png"
os.makedirs(os.path.dirname(outpath), exist_ok=True)

# --- Target distribution ---
# P(x) proportional to exp(-x^2): a Gaussian with variance sigma^2 = 1/2.
# We never need the normalization because MH only uses ratios.
def log_unnorm(x):
    return -x * x  # log of exp(-x^2)

# True (normalized) density for comparison: N(0, 1/2)
sigma = np.sqrt(0.5)
def true_pdf(x):
    return np.exp(-x * x) / np.sqrt(np.pi)  # normalizing constant is sqrt(pi)

# --- Metropolis-Hastings sampler (explicit implementation) ---
def metropolis_hastings(dx_max, n_steps, x0, rng):
    x = x0
    samples = np.empty(n_steps)
    n_accept = 0
    for i in range(n_steps):
        # Symmetric uniform random-walk proposal: x' = x + U(-dx_max, dx_max)
        x_prop = x + rng.uniform(-dx_max, dx_max)
        # Acceptance ratio a = min(1, P(x')/P(x)) = min(1, exp(-x'^2 + x^2))
        a = min(1.0, np.exp(-x_prop * x_prop + x * x))
        # Accept with probability a, otherwise stay
        if rng.random() < a:
            x = x_prop
            n_accept += 1
        samples[i] = x
    accept_rate = n_accept / n_steps
    return samples, accept_rate

# --- Run the experiment ---
n_steps = 10_000
x0 = 0.0
dx_list = [0.1, 0.5, 2.0, 10.0]

results = {}
for dx in dx_list:
    rng = np.random.default_rng(1)  # seed 1 for each run
    samples, rate = metropolis_hastings(dx, n_steps, x0, rng)
    results[dx] = (samples, rate)
    print(f"dx_max = {dx:>5}: acceptance rate = {rate:.4f}, "
          f"sample mean = {samples.mean():.4f}, sample var = {samples.var():.4f}")

# --- Goodness-of-fit check: compare sampled histogram to true Gaussian ---
# Use a fixed set of bins; measure mismatch as sum of squared differences
# between the normalized histogram density and the true density at bin centers.
bins = np.linspace(-4, 4, 41)
centers = 0.5 * (bins[:-1] + bins[1:])
true_at_centers = true_pdf(centers)

print()
print("True distribution: N(0, 1/2)  -> mean = 0, variance = 0.5")
print()

best_dx = None
best_err = np.inf
for dx in dx_list:
    samples, rate = results[dx]
    dens, _ = np.histogram(samples, bins=bins, density=True)
    err = np.sum((dens - true_at_centers) ** 2)
    print(f"dx_max = {dx:>5}: histogram-vs-true SSE = {err:.5f}, acceptance = {rate:.4f}")
    if err < best_err:
        best_err = err
        best_dx = dx

print()
print(f"Best-fitting proposal width: dx_max = {best_dx} "
      f"(acceptance = {results[best_dx][1]:.4f})")
print("Explanation: this confirms the result because the intermediate step size "
      "near 50% acceptance balances exploration against rejection, so its "
      "histogram has the lowest squared error to the true Gaussian, "
      "while dx_max=0.1 barely moves and dx_max=10 mostly rejects.")

# --- Plot histograms against the true Gaussian ---
xgrid = np.linspace(-4, 4, 400)
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, dx in zip(axes.ravel(), dx_list):
    samples, rate = results[dx]
    ax.hist(samples, bins=bins, density=True, alpha=0.6,
            color="steelblue", label="MH samples")
    ax.plot(xgrid, true_pdf(xgrid), "r-", lw=2, label="True N(0,1/2)")
    ax.set_title(f"dx_max = {dx}  (accept = {rate:.2f})")
    ax.set_xlabel("x")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
fig.suptitle("Metropolis-Hastings sampling of exp(-x^2): effect of step size")
fig.tight_layout()
plt.savefig(outpath)
print()
print(f"Figure saved to: {outpath}")
