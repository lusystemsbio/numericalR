import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Target model: standard normal N(0, 1)
mu, sigma = 0.0, 1.0

n = 10000  # number of accepted normal draws we want per output stream
rng = np.random.default_rng(0)

# Two output streams of independent normals produced by the polar method
stream_x = []  # first normal from each accepted point
stream_y = []  # second normal from each accepted point

n_trials = 0  # total uniform points drawn (accepted + rejected)
while len(stream_x) < n:
    n_trials += 1
    # 1) draw a point uniform in the square [-1, 1]^2
    x = rng.uniform(-1.0, 1.0)
    y = rng.uniform(-1.0, 1.0)
    # 2) compute squared radius
    r2 = x * x + y * y
    # 3) reject unless R2 lies in the open unit disk (0, 1)
    if r2 <= 0.0 or r2 >= 1.0:
        continue
    # 4) common scaling factor for the accepted point
    factor = np.sqrt(-2.0 * np.log(r2) / r2)
    # 5) two independent standard normals
    stream_x.append(x * factor)
    stream_y.append(y * factor)

stream_x = np.asarray(stream_x)
stream_y = np.asarray(stream_y)

# Combined sample used for the main density histogram
samples = stream_x  # the primary test set of n accepted draws

# Reported numerical results
print(f"Requested accepted draws (n): {n}")
print(f"Total uniform points drawn (trials): {n_trials}")
print(f"Empirical acceptance rate: {n / n_trials:.4f}")
print(f"Theoretical acceptance rate (pi/4): {np.pi / 4:.4f}")
print(f"Stream X mean: {stream_x.mean():.4f}  (target {mu:.4f})")
print(f"Stream X std : {stream_x.std(ddof=1):.4f}  (target {sigma:.4f})")
print(f"Stream Y mean: {stream_y.mean():.4f}  (target {mu:.4f})")
print(f"Stream Y std : {stream_y.std(ddof=1):.4f}  (target {sigma:.4f})")

# Standard normal density for overlay
grid = np.linspace(-4.0, 4.0, 400)
pdf = (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((grid - mu) / sigma) ** 2)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: main test — density histogram of n accepted draws vs N(0,1)
axes[0].hist(samples, bins=50, density=True, alpha=0.6,
             color="steelblue", edgecolor="white", label="samples (density)")
axes[0].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[0].set_title(f"Polar Box-Muller: {n} accepted draws")
axes[0].set_xlabel("value")
axes[0].set_ylabel("density")
axes[0].legend()

# Right: separate check — both output streams against the same bell curve
axes[1].hist(stream_x, bins=50, density=True, alpha=0.5,
             color="steelblue", edgecolor="white", label="stream X")
axes[1].hist(stream_y, bins=50, density=True, alpha=0.5,
             color="seagreen", edgecolor="white", label="stream Y")
axes[1].plot(grid, pdf, "r-", lw=2, label="N(0,1) density")
axes[1].set_title("Both output streams vs N(0,1)")
axes[1].set_xlabel("value")
axes[1].set_ylabel("density")
axes[1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/6A.4.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because each output stream's density histogram closely tracks "
      "the standard normal curve (and the sample means/stds match 0 and 1), the "
      "transformation is confirmed to convert uniform inputs into correct N(0,1) samples.")
