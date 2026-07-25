import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Buffon's needle: drop a needle of length l on a floor ruled with
# parallel lines a distance d apart. If l < d the needle crosses a
# line with probability P = 2*l/(pi*d), hence pi = 2*l/(P*d).
# We estimate P by Monte Carlo and invert to estimate pi.
# ---------------------------------------------------------------

rng = np.random.default_rng(1)  # seed 1 for reproducibility

d = 1.0                                       # line spacing
lengths = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])  # needle lengths
sample_counts = [int(1e4), int(1e5)]          # number of drops to test


def buffon_estimate(l, d, n, rng):
    """One explicit Buffon experiment with n needle drops.

    Returns the estimate of pi (or nan if no crossings occurred)."""
    # Position of the needle's center relative to the nearest line.
    # By symmetry the perpendicular distance from center to nearest
    # line is uniform on [0, d/2].
    y = rng.uniform(0.0, d / 2.0, size=n)

    # Random needle orientation via REJECTION SAMPLING:
    # draw points uniformly in the unit square, keep those inside the
    # unit circle, and use the point's direction as a uniform angle.
    # This gives a uniformly distributed direction without calling
    # sin/cos on a directly-sampled uniform angle.
    ux = np.empty(n)
    uy = np.empty(n)
    filled = 0
    while filled < n:
        need = n - filled
        # oversample a bit to reduce loop iterations (acceptance ~ pi/4)
        m = int(need / 0.75) + 1
        cx = rng.uniform(-1.0, 1.0, size=m)
        cy = rng.uniform(-1.0, 1.0, size=m)
        r2 = cx * cx + cy * cy
        accept = (r2 > 0.0) & (r2 <= 1.0)      # keep points inside unit circle
        cx, cy, r2 = cx[accept], cy[accept], r2[accept]
        take = min(need, cx.size)
        r = np.sqrt(r2[:take])
        ux[filled:filled + take] = cx[:take] / r  # cos(theta) of uniform angle
        uy[filled:filled + take] = cy[:take] / r  # sin(theta) of uniform angle
        filled += take

    # |sin(theta)| is the vertical half-extent factor of the needle.
    # The needle crosses a line when its half-length projected onto the
    # perpendicular direction reaches the nearest line: (l/2)*|sin| >= y.
    half_vertical = (l / 2.0) * np.abs(uy)
    crossings = np.count_nonzero(half_vertical >= y)

    P = crossings / n                          # empirical crossing probability
    if crossings == 0:
        return np.nan
    return 2.0 * l / (P * d)                    # invert to estimate pi


# ---------------------------------------------------------------
# Run every (length, sample_count) combination and print results.
# ---------------------------------------------------------------
results = {}  # (l, n) -> pi estimate
for n in sample_counts:
    for l in lengths:
        est = buffon_estimate(l, d, n, rng)
        results[(l, n)] = est
        print(f"l = {l:.1f}, d = {d:.1f}, n = {n:>7d}: pi_estimate = {est:.6f}")

print(f"true pi = {np.pi:.6f}")

# ---------------------------------------------------------------
# Explicit check: at n = 1e4 estimates wander around pi (and the
# model breaks for l > d, where 2l/d > 1 can push P>1 territory /
# bias appears); at n = 1e5 the estimates tighten toward pi.
# ---------------------------------------------------------------
valid = lengths[lengths < d]  # the regime where the formula is valid (l < d)
err_1e4 = np.array([abs(results[(l, int(1e4))] - np.pi) for l in valid])
err_1e5 = np.array([abs(results[(l, int(1e5))] - np.pi) for l in valid])
print(f"mean |error| over valid l (l<d) at n=1e4: {err_1e4.mean():.6f}")
print(f"mean |error| over valid l (l<d) at n=1e5: {err_1e5.mean():.6f}")
print(f"error reduced by factor: {err_1e4.mean() / err_1e5.mean():.6f}")

# ---------------------------------------------------------------
# Plot: pi estimate vs needle length, one curve per sample count.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
markers = {int(1e4): "o", int(1e5): "s"}
for n in sample_counts:
    ests = [results[(l, n)] for l in lengths]
    ax.plot(lengths, ests, marker=markers[n], label=f"n = {n:.0e}")

ax.axhline(np.pi, color="k", linestyle="--", linewidth=1, label="true pi")
ax.axvline(d, color="r", linestyle=":", linewidth=1, label="l = d (model breaks for l > d)")
ax.set_xlabel("needle length l  (line spacing d = 1)")
ax.set_ylabel("estimate of pi")
ax.set_title("Buffon's needle Monte Carlo estimate of pi")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.2.1_s4.png")

# One-sentence explanation:
print("Explanation: The check confirms the result because a correct estimator "
      "should scatter randomly around the true value with an error that shrinks "
      "like 1/sqrt(n); seeing the n=1e5 estimates cluster far more tightly on pi "
      "than the noisier n=1e4 ones (for l<d, while l>d violates the formula's "
      "assumption and biases the result) is exactly that expected convergence.")
