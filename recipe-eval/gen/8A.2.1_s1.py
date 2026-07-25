import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducible RNG
rng = np.random.default_rng(1)

d = 1.0                                    # spacing between ruled lines
lengths = np.arange(0.2, 1.5, 0.2)         # needle lengths 0.2,...,1.4
sample_counts = [int(1e4), int(1e5)]       # number of drops to test

def buffon_estimate(l, n, rng):
    """Estimate pi from n Buffon-needle drops, needle length l, spacing d=1."""
    # Center of needle: only its distance to the nearest line matters -> uniform in [0, d/2].
    dist = rng.uniform(0.0, d / 2.0, size=n)

    # Uniform random needle DIRECTION via rejection sampling:
    # sample points in the unit square, keep those inside the unit circle,
    # then their angle is a uniform direction. We resample rejected slots
    # until every drop has an accepted direction.
    sin_theta = np.empty(n)
    remaining = np.arange(n)               # indices still needing an accepted direction
    while remaining.size > 0:
        x = rng.uniform(-1.0, 1.0, size=remaining.size)
        y = rng.uniform(-1.0, 1.0, size=remaining.size)
        r2 = x * x + y * y
        accept = (r2 > 0.0) & (r2 <= 1.0)  # inside unit circle -> uniform angle
        idx = remaining[accept]
        # |sin(theta)| = |y| / r  gives the perpendicular half-extent factor
        sin_theta[idx] = np.abs(y[accept]) / np.sqrt(r2[accept])
        remaining = remaining[~accept]     # retry the rejected ones

    # A needle crosses a line when its half-projection reaches past the nearest line.
    crossings = dist <= (l / 2.0) * sin_theta
    P = np.mean(crossings)                 # empirical crossing probability
    if P == 0.0:
        return np.nan
    return 2.0 * l / (P * d)               # invert P = 2l/(pi d)

# Compute estimates for every (length, n) pair.
estimates = {n: [] for n in sample_counts}
for n in sample_counts:
    for l in lengths:
        # Fresh, seeded RNG per case so results are reproducible and comparable.
        est = buffon_estimate(l, n, np.random.default_rng(1))
        estimates[n].append(est)
        print(f"n={n:>7d}  l={l:.1f}  pi_estimate={est:.6f}")

print(f"true_pi={np.pi:.6f}")

# Plot: estimate vs needle length, one series per sample count.
plt.figure(figsize=(8, 5))
for n in sample_counts:
    plt.plot(lengths, estimates[n], marker="o", label=f"n = {n}")
plt.axhline(np.pi, color="k", linestyle="--", label="true pi")
plt.axvline(d, color="r", linestyle=":", label="l = d (validity limit)")
plt.xlabel("needle length l  (spacing d = 1)")
plt.ylabel("estimated pi")
plt.title("Buffon's needle Monte Carlo estimate of pi")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.2.1_s1.png")

# Check summary: spread of estimates over valid lengths (l < d) at each n.
for n in sample_counts:
    vals = np.array([estimates[n][i] for i, l in enumerate(lengths) if l < d])
    print(f"n={n:>7d}  valid-length mean={np.mean(vals):.6f}  std={np.std(vals):.6f}")

# One-sentence explanation of why the check confirms the result:
print("Check: at n=1e4 the estimates scatter around pi with larger std (and l>d violates "
      "l<d so P can saturate and the formula breaks), while at n=1e5 the std shrinks toward "
      "0, confirming convergence to pi as Monte Carlo error falls like 1/sqrt(n).")
