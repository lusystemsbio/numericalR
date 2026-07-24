import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

rng = np.random.default_rng(0)

n = 10000  # number of accepted draws we want per stream

# We collect two independent normal streams from the polar Box-Muller method.
x_normals = []  # first output stream
y_normals = []  # second output stream

# Keep drawing until we have n accepted points.
while len(x_normals) < n:
    # Draw a candidate point uniform in the square [-1, 1]^2.
    x = rng.uniform(-1.0, 1.0)
    y = rng.uniform(-1.0, 1.0)
    # Squared radius of the candidate point.
    R2 = x * x + y * y
    # Accept only if it lies strictly inside the unit disk (excluding origin).
    if 0.0 < R2 < 1.0:
        # Common scaling factor from the polar transform.
        factor = np.sqrt(-2.0 * np.log(R2) / R2)
        # Two independent standard normals.
        x_normals.append(x * factor)
        y_normals.append(y * factor)

x_normals = np.array(x_normals)
y_normals = np.array(y_normals)

# Report basic sample statistics for each stream (should be ~0 mean, ~1 std).
print(f"Number of accepted draws per stream: {n}")
print(f"Stream 1 (x) sample mean: {x_normals.mean():.6f}")
print(f"Stream 1 (x) sample std : {x_normals.std(ddof=1):.6f}")
print(f"Stream 2 (y) sample mean: {y_normals.mean():.6f}")
print(f"Stream 2 (y) sample std : {y_normals.std(ddof=1):.6f}")

# Standard normal density for overlay.
def normal_pdf(t):
    return np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

grid = np.linspace(-4.0, 4.0, 400)
pdf = normal_pdf(grid)

# Plot: main histogram (stream 1) and a separate check of stream 2.
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(x_normals, bins=50, density=True, alpha=0.6,
             color="steelblue", label="samples (stream 1)")
axes[0].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[0].set_title("Density histogram vs. standard normal")
axes[0].set_xlabel("value")
axes[0].set_ylabel("density")
axes[0].legend()

axes[1].hist(y_normals, bins=50, density=True, alpha=0.6,
             color="seagreen", label="samples (stream 2)")
axes[1].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[1].set_title("Check: second output stream")
axes[1].set_xlabel("value")
axes[1].set_ylabel("density")
axes[1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.4.1_s2.png")

# Explanation of why the check confirms the result:
print("Explanation: Both output streams' density histograms tracking the same "
      "bell-shaped N(0,1) curve (with mean ~0 and std ~1) confirms the method, "
      "because the polar Box-Muller transform is only valid if each of its two "
      "returned values is independently standard normal.")
