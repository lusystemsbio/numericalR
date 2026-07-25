import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Fix seed for reproducibility of the numerical checks
rng = np.random.default_rng(0)

# Number of samples
n = 10000

# Step 1: draw n uniform random numbers on (0, 1)
x = rng.uniform(0.0, 1.0, size=n)

# Step 2: inverse-transform sampling of the exponential distribution.
# For a uniform x on (0,1), y = -ln(x) is Exponential(rate=1) on y >= 0.
# (This works because 1 - x is also uniform, so -ln(x) has the same law as -ln(1-x),
#  the standard inverse CDF of the unit-rate exponential.)
y = -np.log(x)

# Report basic statistics: for Exponential(1) the mean and std are both 1.
print(f"Sample size n: {n}")
print(f"Sample mean of y: {np.mean(y)}")
print(f"Sample std of y:  {np.std(y)}")
print(f"Theoretical mean (Exp rate=1): {1.0}")
print(f"Theoretical std  (Exp rate=1): {1.0}")
print(f"Sample minimum y: {np.min(y)}")
print(f"Sample maximum y: {np.max(y)}")

# Fraction of samples below 1: should approach 1 - e^{-1} for a decaying exponential.
frac_below_1 = np.mean(y < 1.0)
print(f"Fraction of y < 1 (observed): {frac_below_1}")
print(f"Fraction of y < 1 (expected 1 - e^-1): {1.0 - np.exp(-1.0)}")

# Build the figure: density histogram of the transformed samples vs the exponential density.
fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=(11, 4.5))

# Left panel: the flat uniform INPUT as a density histogram (the "before" check).
ax_in.hist(x, bins=40, range=(0, 1), density=True, color="tab:gray",
           edgecolor="black", alpha=0.7, label="uniform input")
ax_in.axhline(1.0, color="red", lw=2, label="uniform density = 1")
ax_in.set_title("Input: flat uniform on (0,1)")
ax_in.set_xlabel("x")
ax_in.set_ylabel("density")
ax_in.legend()

# Right panel: the transformed OUTPUT as a density histogram vs the exponential density.
ax_out.hist(y, bins=60, density=True, color="tab:blue",
            edgecolor="black", alpha=0.7, label="transformed samples y = -ln(x)")
grid = np.linspace(0.0, np.max(y), 400)
pdf = np.exp(-grid)  # Exponential(rate=1) density: f(y) = e^{-y}
ax_out.plot(grid, pdf, color="red", lw=2, label="exponential density e^{-y}")
ax_out.set_title("Output: decaying exponential density")
ax_out.set_xlabel("y")
ax_out.set_ylabel("density")
ax_out.legend()

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6A.3.1_s3.png")

# Explanation of the check:
print("Explanation: Seeing the flat uniform input turn into a histogram that hugs the "
      "e^{-y} curve confirms the transform works, because inverse-transform sampling is "
      "designed to reshape the constant uniform density into exactly the exponential's "
      "decaying density.")
