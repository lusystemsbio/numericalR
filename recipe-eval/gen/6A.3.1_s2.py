import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(42)

# --- Inverse-transform sampling of the exponential distribution ---
# Target: exponential density f(y) = exp(-y) on y >= 0 (rate = 1).
# The CDF is F(y) = 1 - exp(-y); inverting gives y = -ln(1 - x).
# Since x is uniform on (0,1), (1 - x) is also uniform on (0,1),
# so we may use the equivalent simplified rule y = -ln(x).

n = 10000

# Step 1: draw n uniform samples on (0, 1)
x = rng.uniform(0.0, 1.0, size=n)

# Step 2: apply the inverse-CDF transform explicitly
y = -np.log(x)   # transformed samples ~ Exponential(1)

# Report basic numerical checks against the exponential(1) distribution.
# For Exponential(rate=1): mean = 1, variance = 1.
print(f"Number of samples n: {n}")
print(f"Sample mean of y (expected 1.0): {np.mean(y):.6f}")
print(f"Sample variance of y (expected 1.0): {np.var(y):.6f}")
print(f"Sample min of y (>= 0): {np.min(y):.6f}")
print(f"Sample max of y: {np.max(y):.6f}")
print(f"Uniform input mean (expected 0.5): {np.mean(x):.6f}")
print(f"Uniform input variance (expected {1/12:.6f}): {np.var(x):.6f}")

# --- Plot: density histogram of samples vs. the exponential density ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: the flat uniform input as a density histogram
ax1.hist(x, bins=40, density=True, color="steelblue",
         edgecolor="white", alpha=0.8, label="uniform input")
ax1.axhline(1.0, color="red", lw=2, label="uniform density = 1")
ax1.set_title("Input: flat uniform on (0, 1)")
ax1.set_xlabel("x")
ax1.set_ylabel("density")
ax1.legend()

# Right panel: transformed samples vs. the exponential(1) density
counts, edges, _ = ax2.hist(y, bins=50, density=True, color="darkorange",
                            edgecolor="white", alpha=0.8,
                            label="transformed samples")
grid = np.linspace(0, np.max(y), 400)
ax2.plot(grid, np.exp(-grid), color="black", lw=2,
         label=r"exponential density $e^{-y}$")
ax2.set_title("Output: decaying exponential after y = -ln(x)")
ax2.set_xlabel("y")
ax2.set_ylabel("density")
ax2.legend()

# Report how well the histogram matches the density at bin centers.
centers = 0.5 * (edges[:-1] + edges[1:])
max_abs_dev = np.max(np.abs(counts - np.exp(-centers)))
print(f"Max abs deviation (histogram vs exp density at bin centers): {max_abs_dev:.6f}")

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6A.3.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the flat, constant-density uniform input is turned "
      "into the characteristic monotonically decaying exponential density e^{-y} "
      "that closely matches the overlaid theoretical curve, the transform y = -ln(x) "
      "is confirmed to correctly generate exponentially distributed samples.")
