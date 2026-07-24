import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(12345)

n = 10000  # number of accepted normal draws we want in the FIRST output stream

# Storage for the two independent output streams of normals
stream1 = []  # gets x * factor
stream2 = []  # gets y * factor

# Polar (Marsaglia) Box-Muller, implemented explicitly ---------------------
while len(stream1) < n:
    # 1) draw a point uniform in the square [-1, 1]^2 from uniforms on [0,1)
    x = 2.0 * rng.random() - 1.0
    y = 2.0 * rng.random() - 1.0

    # 2) compute squared radius
    r2 = x * x + y * y

    # 3) rejection step: keep only points strictly inside the unit disk,
    #    excluding the origin (need R2 in (0,1))
    if r2 <= 0.0 or r2 >= 1.0:
        continue  # reject and try again

    # 4) common factor sqrt(-2 * ln(R2) / R2)
    factor = np.sqrt(-2.0 * np.log(r2) / r2)

    # 5) two independent standard normals
    stream1.append(x * factor)
    stream2.append(y * factor)

stream1 = np.array(stream1)
stream2 = np.array(stream2)

# Combined set of samples (all accepted normals from both streams)
combined = np.concatenate([stream1, stream2])

# Standard normal density for overlay
def normal_pdf(t):
    return np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

# ---- Numerical summaries -------------------------------------------------
print("Requested accepted draws per stream (n):", n)
print("Stream 1 size:", stream1.size)
print("Stream 2 size:", stream2.size)
print("Stream 1 mean:", np.mean(stream1))
print("Stream 1 std (ddof=1):", np.std(stream1, ddof=1))
print("Stream 2 mean:", np.mean(stream2))
print("Stream 2 std (ddof=1):", np.std(stream2, ddof=1))
print("Combined mean:", np.mean(combined))
print("Combined std (ddof=1):", np.std(combined, ddof=1))
print("Target mean:", 0.0)
print("Target std:", 1.0)

# ---- Plot: density histograms vs standard normal density -----------------
grid = np.linspace(-4.5, 4.5, 400)
pdf = normal_pdf(grid)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: combined samples against the standard normal density
axes[0].hist(combined, bins=60, density=True, alpha=0.6,
             color="steelblue", edgecolor="white", label="samples (both streams)")
axes[0].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[0].set_title("Density histogram of samples vs N(0,1)")
axes[0].set_xlabel("value")
axes[0].set_ylabel("density")
axes[0].legend()

# Right panel: separate check that BOTH output streams are bell-shaped normal
axes[1].hist(stream1, bins=60, density=True, alpha=0.5,
             color="seagreen", edgecolor="white", label="stream 1 (x)")
axes[1].hist(stream2, bins=60, density=True, alpha=0.5,
             color="orange", edgecolor="white", label="stream 2 (y)")
axes[1].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[1].set_title("Both output streams follow the normal density")
axes[1].set_xlabel("value")
axes[1].set_ylabel("density")
axes[1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.4.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the density histograms of both independent output "
      "streams closely track the theoretical N(0,1) bell curve (with sample means "
      "near 0 and standard deviations near 1), the transform is correctly producing "
      "standard-normal variates from uniform inputs.")
