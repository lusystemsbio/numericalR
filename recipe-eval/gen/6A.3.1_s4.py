import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
np.random.seed(0)

# Number of samples
n = 10000

# Step 1: draw n uniform random numbers on (0, 1)
# (avoid exactly 0 so that -ln(x) stays finite)
x = np.random.uniform(0.0, 1.0, size=n)

# Step 2: inverse-transform sampling of the exponential distribution.
# For a uniform x on (0, 1), the exponential (rate 1) sample is y = -ln(x).
# This is the explicit inverse of the exponential CDF F(y) = 1 - exp(-y),
# and since 1 - x is also uniform on (0, 1) we may use x directly.
y = -np.log(x)

# Basic numerical checks: sample mean/variance/min should match Exp(1) (mean=1, var=1, min>=0)
print("Number of samples n:", n)
print("Sample mean of y (expected 1.0):", np.mean(y))
print("Sample variance of y (expected 1.0):", np.var(y))
print("Minimum of y (expected >= 0):", np.min(y))
print("Maximum of y:", np.max(y))

# Confirm the uniform input is indeed flat (mean ~ 0.5, constant density ~ 1)
print("Sample mean of uniform input x (expected 0.5):", np.mean(x))

# Step 3: build density histogram of the transformed samples
plt.figure(figsize=(9, 5))
counts, bins, _ = plt.hist(
    y, bins=50, density=True, color="steelblue",
    alpha=0.6, edgecolor="white", label="transformed samples (density)"
)

# Report the peak of the histogram density (should be near 1.0, the value of the exponential density at y=0)
print("Peak histogram density (expected near 1.0):", np.max(counts))

# Step 4: overlay the true exponential density f(y) = exp(-y) for y >= 0
grid = np.linspace(0.0, np.max(y), 400)
pdf = np.exp(-grid)
plt.plot(grid, pdf, color="crimson", linewidth=2.0, label=r"exponential density $e^{-y}$")

plt.title("Inverse-transform sampling: uniform -> exponential")
plt.xlabel("y")
plt.ylabel("density")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6A.3.1_s4.png")

# Explanation of the check:
print(
    "Why the check confirms the result: a flat (constant-density) uniform input "
    "turning into the decaying e^{-y} curve shows the transform y = -ln(x) maps the "
    "uniform measure exactly onto the exponential density, which is the defining property "
    "of correct inverse-transform sampling."
)
