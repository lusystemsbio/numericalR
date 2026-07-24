import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(12345)

n = 10000              # number of accepted draws we want (produces 2n normals)
out_x = []             # first output stream of normals
out_y = []             # second output stream of normals
attempts = 0           # count total uniform-square draws (for acceptance rate)

# Polar (Marsaglia) Box-Muller, implemented explicitly:
while len(out_x) < n:
    attempts += 1
    # draw a point uniform in the square [-1, 1]^2 from uniforms on [0,1]
    x = 2.0 * rng.random() - 1.0
    y = 2.0 * rng.random() - 1.0
    R2 = x * x + y * y
    # keep the point only if it lies strictly inside the unit disk, R2 in (0, 1)
    if 0.0 < R2 < 1.0:
        factor = np.sqrt(-2.0 * np.log(R2) / R2)  # common transform factor
        out_x.append(x * factor)                  # first independent normal
        out_y.append(y * factor)                  # second independent normal

out_x = np.array(out_x)
out_y = np.array(out_y)

# Acceptance rate should be about pi/4 ~ 0.7854 (area of unit disk / area of square)
acceptance_rate = n / attempts

# Sample statistics as a numerical sanity check against N(0, 1)
print("Number of accepted draws (n):", n)
print("Total uniform-square attempts:", attempts)
print("Empirical acceptance rate:", acceptance_rate)
print("Theoretical acceptance rate (pi/4):", np.pi / 4.0)
print("Stream X mean:", out_x.mean())
print("Stream X std:", out_x.std(ddof=1))
print("Stream Y mean:", out_y.mean())
print("Stream Y std:", out_y.std(ddof=1))

# For the main density-histogram check, use the first accepted stream (n samples)
samples = out_x

# Standard normal density for overlay
grid = np.linspace(-4.0, 4.0, 400)
pdf = np.exp(-0.5 * grid ** 2) / np.sqrt(2.0 * np.pi)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: main check -- density histogram of samples vs standard normal density
axes[0].hist(samples, bins=50, density=True, alpha=0.6,
             color="steelblue", edgecolor="white", label="samples (stream X)")
axes[0].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[0].set_title("Density histogram of generated normals vs N(0,1)")
axes[0].set_xlabel("value")
axes[0].set_ylabel("density")
axes[0].legend()

# Right: separate check -- both output streams follow the bell-shaped density
axes[1].hist(out_x, bins=50, density=True, alpha=0.5,
             color="steelblue", edgecolor="white", label="stream X")
axes[1].hist(out_y, bins=50, density=True, alpha=0.5,
             color="seagreen", edgecolor="white", label="stream Y")
axes[1].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[1].set_title("Both output streams vs N(0,1)")
axes[1].set_xlabel("value")
axes[1].set_ylabel("density")
axes[1].legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.4.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the density histograms of both independent output "
      "streams closely track the standard normal curve (with sample means ~0 and "
      "standard deviations ~1), the transformation is correctly producing "
      "independent N(0,1) variates.")
