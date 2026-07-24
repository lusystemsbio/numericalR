import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(0)

# ---- Polar (Marsaglia) Box-Muller: uniforms -> standard normals ----
def polar_box_muller(n, rng):
    """Return two independent N(0,1) streams, each of length n."""
    out1 = np.empty(n)  # first normal from each accepted point
    out2 = np.empty(n)  # second normal from each accepted point
    filled = 0
    attempts = 0
    while filled < n:
        # draw a point uniform in the square [-1, 1]^2
        x = rng.uniform(-1.0, 1.0)
        y = rng.uniform(-1.0, 1.0)
        attempts += 1
        r2 = x * x + y * y      # squared radius
        # keep only points strictly inside the unit disk (rejection step)
        if 0.0 < r2 < 1.0:
            factor = np.sqrt(-2.0 * np.log(r2) / r2)  # common polar factor
            out1[filled] = x * factor  # one standard normal
            out2[filled] = y * factor  # second, independent standard normal
            filled += 1
    return out1, out2, attempts

n = 10000
z1, z2, attempts = polar_box_muller(n, rng)

# ---- Basic numerical checks ----
print("Accepted draws (n):", n)
print("Total square draws (attempts):", attempts)
print("Acceptance rate:", n / attempts)
print("Theoretical acceptance rate (pi/4):", np.pi / 4.0)
print("Stream 1 mean:", np.mean(z1))
print("Stream 1 std:", np.std(z1))
print("Stream 2 mean:", np.mean(z2))
print("Stream 2 std:", np.std(z2))
combined = np.concatenate([z1, z2])
print("Combined mean:", np.mean(combined))
print("Combined std:", np.std(combined))
print("Correlation between streams:", np.corrcoef(z1, z2)[0, 1])

# ---- Plot: density histograms vs standard normal density ----
grid = np.linspace(-4, 4, 400)
pdf = np.exp(-0.5 * grid**2) / np.sqrt(2.0 * np.pi)  # N(0,1) density

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

axes[0].hist(z1, bins=50, density=True, alpha=0.6, color="steelblue",
             label="samples (stream 1)")
axes[0].plot(grid, pdf, "r-", lw=2, label="N(0,1) pdf")
axes[0].set_title("Stream 1 (n=10000)")
axes[0].legend()

axes[1].hist(z2, bins=50, density=True, alpha=0.6, color="seagreen",
             label="samples (stream 2)")
axes[1].plot(grid, pdf, "r-", lw=2, label="N(0,1) pdf")
axes[1].set_title("Stream 2 (n=10000)")
axes[1].legend()

axes[2].hist(combined, bins=60, density=True, alpha=0.6, color="darkorange",
             label="samples (both streams)")
axes[2].plot(grid, pdf, "r-", lw=2, label="N(0,1) pdf")
axes[2].set_title("Both streams combined")
axes[2].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.4.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because each output stream's density histogram closely "
      "overlays the standard normal curve (and their means/stds are ~0 and ~1), "
      "the polar Box-Muller transform is confirmed to convert uniforms into "
      "independent N(0,1) draws.")
