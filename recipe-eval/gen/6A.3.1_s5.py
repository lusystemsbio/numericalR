import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reproducibility
rng = np.random.default_rng(0)

# --- Step 1: draw n uniform random numbers on (0, 1) ---
n = 10000
x = rng.uniform(0.0, 1.0, size=n)  # flat uniform input

# --- Step 2: inverse-transform sampling of the exponential ---
# For a uniform x on (0,1), y = -ln(x) is exponentially distributed (rate 1).
# (This works because if X ~ U(0,1) then -ln(X) ~ Exp(1); we do it explicitly.)
y = -np.log(x)  # transform each uniform draw into an exponential draw

# --- Step 3: summary statistics of the transformed samples ---
# For Exp(1) the theoretical mean and variance are both 1.
print("Number of samples n:", n)
print("Sample mean of y (theory = 1):", np.mean(y))
print("Sample variance of y (theory = 1):", np.var(y))
print("Sample min of y:", np.min(y))
print("Sample max of y:", np.max(y))

# --- Step 4: compare density histogram of y against the exponential density ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: flat uniform input (the "before" picture)
ax1.hist(x, bins=40, range=(0, 1), density=True,
         color="lightsteelblue", edgecolor="white")
ax1.axhline(1.0, color="red", lw=2, label="uniform density = 1")
ax1.set_title("Uniform input x ~ U(0,1) (flat)")
ax1.set_xlabel("x")
ax1.set_ylabel("density")
ax1.legend()

# Right panel: transformed samples vs. exponential density (the "after" picture)
counts, edges, _ = ax2.hist(y, bins=50, density=True,
                            color="lightsteelblue", edgecolor="white",
                            label="samples y = -ln(x)")
grid = np.linspace(0, np.max(y), 400)
pdf = np.exp(-grid)  # Exp(1) density: f(y) = e^{-y} for y >= 0
ax2.plot(grid, pdf, color="red", lw=2, label="exponential density e^{-y}")
ax2.set_title("Transformed output vs. exponential density")
ax2.set_xlabel("y")
ax2.set_ylabel("density")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6A.3.1_s5.png")

# --- Step 5: quantitative check that the transform reshapes the density ---
# Compare the histogram bin heights of y to the exponential density at bin centers.
centers = 0.5 * (edges[:-1] + edges[1:])
model = np.exp(-centers)
max_abs_dev = np.max(np.abs(counts - model))
print("Max absolute deviation (histogram vs. e^{-y}):", max_abs_dev)

# The check confirms the result because a flat uniform input (constant density 1)
# is turned by y = -ln(x) into a monotonically decaying density that overlays the
# analytic exponential curve e^{-y}, which is exactly what correct sampling requires.
print("Explanation: the flat uniform density becomes the decaying e^{-y} curve, "
      "matching the exponential density, so the transform samples Exp(1) correctly.")
