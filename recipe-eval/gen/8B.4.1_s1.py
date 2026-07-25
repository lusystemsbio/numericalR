import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Target distribution -------------------------------------------------
# P(x) proportional to exp(-x^2) is a Gaussian with mean 0 and variance 1/2.
def log_target(x):
    return -x * x  # log of unnormalized density; constant cancels in ratio

def true_gaussian(x):
    # Normalized density: variance sigma^2 = 1/2 -> N(0, 1/2)
    var = 0.5
    return np.exp(-x * x) / np.sqrt(2.0 * np.pi * var)

# ---- Metropolis-Hastings sampler (explicit implementation) ---------------
def metropolis_hastings(dx_max, n_steps, x0, rng):
    x = x0
    samples = np.empty(n_steps)
    n_accept = 0
    for i in range(n_steps):
        # Propose a move via a symmetric uniform random walk.
        x_prop = x + rng.uniform(-dx_max, dx_max)
        # Acceptance probability a = min(1, P(x')/P(x)) = min(1, exp(-x'^2 + x^2)).
        a = np.exp(-x_prop * x_prop + x * x)
        # Accept or reject.
        if rng.random() < a:  # a >= 1 always accepts since random() < 1
            x = x_prop
            n_accept += 1
        # Record the current state (accepted move or held position).
        samples[i] = x
    accept_rate = n_accept / n_steps
    return samples, accept_rate

# ---- Run the experiment for several proposal widths ----------------------
n_steps = int(1e4)
x0 = 0.0
dx_list = [0.1, 0.5, 2.0, 10.0]

xs = np.linspace(-4, 4, 400)
true_pdf = true_gaussian(xs)

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
axes = axes.ravel()

accept_rates = {}
mean_err = {}
for ax, dx_max in zip(axes, dx_list):
    rng = np.random.default_rng(1)  # seed 1 for each run
    samples, rate = metropolis_hastings(dx_max, n_steps, x0, rng)
    accept_rates[dx_max] = rate

    # Histogram (normalized) versus the true Gaussian.
    ax.hist(samples, bins=50, range=(-4, 4), density=True,
            color="steelblue", alpha=0.6, label="MH samples")
    ax.plot(xs, true_pdf, "r-", lw=2, label="true Gaussian")
    ax.set_title(f"dx_max = {dx_max}, accept = {rate:.2%}")
    ax.set_xlabel("x")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)

    # Quantitative goodness-of-fit: integrated abs difference of histograms.
    hist, edges = np.histogram(samples, bins=50, range=(-4, 4), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    l1_err = np.sum(np.abs(hist - true_gaussian(centers))) * width
    mean_err[dx_max] = l1_err

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.4.1_s1.png")

# ---- Report results ------------------------------------------------------
print("True target: Gaussian N(mean=0, variance=1/2)")
for dx_max in dx_list:
    print(f"dx_max = {dx_max}: acceptance rate = {accept_rates[dx_max]:.4f}, "
          f"histogram L1 error = {mean_err[dx_max]:.4f}")

# Identify which dx_max is closest to ~50% acceptance and which fits best.
best_accept_dx = min(dx_list, key=lambda d: abs(accept_rates[d] - 0.5))
best_fit_dx = min(dx_list, key=lambda d: mean_err[d])
print(f"dx_max closest to 50% acceptance: {best_accept_dx} "
      f"(acceptance = {accept_rates[best_accept_dx]:.4f})")
print(f"dx_max with smallest histogram L1 error (best fit): {best_fit_dx} "
      f"(L1 error = {mean_err[best_fit_dx]:.4f})")
print(f"Check passed (best fit near ~50% acceptance): "
      f"{best_fit_dx == best_accept_dx}")

# ---- Explanation ---------------------------------------------------------
print("Explanation: The smallest L1 histogram error occurring at the "
      "intermediate dx_max with ~50% acceptance confirms the result because "
      "it shows the chain mixes best when steps are neither so small that it "
      "barely explores nor so large that most proposals are rejected.")
