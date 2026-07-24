import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(0)

# --- Inverse-transform sampling of the exponential distribution ---
# Target: exponential density f(y) = exp(-y) for y >= 0 (rate lambda = 1).
# Its CDF is F(y) = 1 - exp(-y). The inverse CDF gives the transform.
# For a uniform x on (0,1), y = -ln(x) is exponentially distributed.
# (Using x directly rather than 1-x is valid since both are uniform on (0,1).)

n = 10000                                  # number of draws
x = rng.uniform(0.0, 1.0, size=n)          # 1) flat uniform input on (0,1)
y = -np.log(x)                             # 2) inverse-transform: y = -ln(x)

# --- Simple numerical checks ---
sample_mean = y.mean()                     # exponential(1) has mean 1
sample_var = y.var()                       # exponential(1) has variance 1
print(f"Number of samples n: {n}")
print(f"Uniform input mean (expected 0.5): {x.mean():.6f}")
print(f"Uniform input variance (expected {1/12:.6f}): {x.var():.6f}")
print(f"Transformed sample mean (expected 1.0): {sample_mean:.6f}")
print(f"Transformed sample variance (expected 1.0): {sample_var:.6f}")
print(f"Transformed sample min (should be >= 0): {y.min():.6f}")
print(f"Transformed sample max: {y.max():.6f}")

# --- Plot: density histogram of samples vs. exponential density ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: the flat uniform input, shown as a density histogram
ax1.hist(x, bins=40, density=True, color="tab:green", alpha=0.6,
         edgecolor="black", label="uniform input")
ax1.axhline(1.0, color="red", lw=2, label="uniform density = 1")
ax1.set_title("Input: flat uniform on (0,1)")
ax1.set_xlabel("x")
ax1.set_ylabel("density")
ax1.legend()

# Right: transformed samples vs. the target exponential density
counts, bins, _ = ax2.hist(y, bins=50, density=True, color="tab:blue",
                           alpha=0.6, edgecolor="black",
                           label="transformed samples")
yy = np.linspace(0, y.max(), 400)
ax2.plot(yy, np.exp(-yy), color="red", lw=2,
         label=r"exponential density $e^{-y}$")
ax2.set_title("Output: decaying exponential density")
ax2.set_xlabel("y")
ax2.set_ylabel("density")
ax2.legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.3.1_s1.png")

# Explanation: the check confirms the result because a flat (constant-density)
# uniform input becoming the characteristic monotonically decaying e^{-y} curve,
# matching the overlaid exponential density with sample mean and variance near 1,
# is exactly the behavior the inverse-transform must produce if it is correct.
print("Check: flat uniform input -> decaying exponential output confirms the transform,")
print("because only a correct inverse-CDF map turns constant density into e^{-y} with mean and variance 1.")
