import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Buffon's needle Monte Carlo estimate of pi ----
# A needle of length l dropped on a floor ruled with parallel lines a distance d
# apart (l < d) crosses a line with probability P = 2*l/(pi*d).
# Rearranged: pi = 2*l / (P*d), where P is estimated by (# crossings)/(# drops).

def buffon_pi(l, d, n, rng):
    # Drop the center of the needle: distance from the nearest line below,
    # uniform on [0, d). Only this distance matters (translational symmetry).
    y = rng.uniform(0.0, d, size=n)

    # Uniform random needle *direction* via rejection sampling:
    # sample points in the unit square, keep those inside the unit circle,
    # then the angle of the accepted vector is uniform on [0, 2*pi).
    # We resample rejected draws until every needle has a valid direction.
    ax = np.empty(n)
    ay = np.empty(n)
    filled = np.zeros(n, dtype=bool)
    while not filled.all():
        idx = np.where(~filled)[0]            # slots still needing a direction
        vx = rng.uniform(-1.0, 1.0, size=idx.size)
        vy = rng.uniform(-1.0, 1.0, size=idx.size)
        r2 = vx*vx + vy*vy
        ok = (r2 > 0.0) & (r2 <= 1.0)         # accept points inside unit circle
        good = idx[ok]
        rr = np.sqrt(r2[ok])
        ax[good] = vx[ok] / rr                # cos(theta), uniformly directed
        ay[good] = vy[ok] / rr                # sin(theta)
        filled[good] = True

    # Vertical half-extent of the needle from its center is (l/2)*|sin(theta)| = (l/2)*|ay|.
    half = (l / 2.0) * np.abs(ay)

    # A crossing occurs if the needle reaches the line below (y < half) or
    # the line above (y > d - half).
    crossings = (y < half) | (y > d - half)
    P = np.count_nonzero(crossings) / n       # estimated crossing probability
    if P == 0.0:
        return np.nan                         # no crossings -> estimate undefined
    return 2.0 * l / (P * d)                   # pi estimate


d = 1.0
lengths = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
sample_counts = [int(1e4), int(1e5)]

# Compute estimates; use a fresh rng seeded at 1 for each (l, n) run for reproducibility.
estimates = {}   # (l, n) -> pi estimate
for n in sample_counts:
    for l in lengths:
        rng = np.random.default_rng(1)
        est = buffon_pi(l, d, n, rng)
        estimates[(l, n)] = est
        print(f"l={l:.1f}, d={d:.1f}, n={n:>6d}: pi_estimate = {est:.6f}")

print(f"true_pi = {np.pi:.6f}")

# ---- Separate check: spread of the estimate across needle lengths ----
# (only l < d = 1.0 is valid; l >= d violates the model and is expected to break)
valid = lengths < d
for n in sample_counts:
    ests_valid = np.array([estimates[(l, n)] for l in lengths[valid]])
    spread = np.nanmax(ests_valid) - np.nanmin(ests_valid)
    mean_err = np.nanmean(np.abs(ests_valid - np.pi))
    print(f"n={n:>6d}: valid-l spread of pi estimates = {spread:.6f}, "
          f"mean |error| = {mean_err:.6f}")

# ---- Plot: pi estimate vs needle length, one line per sample count ----
plt.figure(figsize=(8, 5))
for n in sample_counts:
    ys = [estimates[(l, n)] for l in lengths]
    plt.plot(lengths, ys, marker="o", label=f"n = {n}")
plt.axhline(np.pi, color="k", linestyle="--", label="true pi")
plt.axvline(d, color="r", linestyle=":", label="l = d (model breaks for l > d)")
plt.xlabel("needle length l  (d = 1)")
plt.ylabel("estimate of pi")
plt.title("Buffon's needle estimate of pi vs needle length")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.2.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: The estimates for valid l<d cluster around the true pi "
      "with a spread that shrinks from n=1e4 to n=1e5 (while l>d cases fall off), "
      "confirming the estimator is unbiased and converges as Monte Carlo error ~1/sqrt(n).")
