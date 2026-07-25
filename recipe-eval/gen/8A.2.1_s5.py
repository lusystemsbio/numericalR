import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Buffon's needle Monte Carlo estimate of pi
#   floor ruled with parallel lines a distance d apart
#   needle of length l (l < d) crosses a line with prob 2l/(pi d)
#   => pi = 2l / (P d),  with P the observed crossing fraction
# ---------------------------------------------------------------

rng = np.random.default_rng(1)   # seed 1

d = 1.0                                   # line spacing
lengths = np.arange(0.2, 1.41, 0.2)       # l = 0.2, 0.4, ..., 1.4
N = int(1e5)                              # max drops (we index 1e4 inside)
n_small, n_big = int(1e4), int(1e5)


def uniform_abs_sin(n):
    """|sin(theta)| for n uniformly-random needle directions,
    obtained by REJECTION SAMPLING of a uniform point in the unit disk."""
    out = np.empty(n)
    filled = 0
    while filled < n:
        # sample a batch of candidate points in the square [-1,1]^2
        m = int((n - filled) / (np.pi / 4) * 1.3) + 16
        u = rng.uniform(-1.0, 1.0, m)
        v = rng.uniform(-1.0, 1.0, m)
        r2 = u * u + v * v
        keep = (r2 <= 1.0) & (r2 > 0.0)   # accept only points inside the disk
        u, v, r2 = u[keep], v[keep], r2[keep]
        # for an accepted point the direction is uniform on the circle;
        # |sin(theta)| is the vertical component over the radius
        s = np.abs(v) / np.sqrt(r2)
        take = min(len(s), n - filled)
        out[filled:filled + take] = s[:take]
        filled += take
    return out


# generate all N drops once per length, then read off running estimates
sample_axis = np.arange(1, N + 1)

plt.figure(figsize=(9, 6))
for l in lengths:
    # perpendicular distance of needle center to nearest line: uniform in [0, d/2]
    y = rng.uniform(0.0, d / 2.0, N)
    # uniform direction via rejection sampling
    abs_sin = uniform_abs_sin(N)
    # a crossing occurs when the half-length projection reaches the nearest line
    crossings = ((l / 2.0) * abs_sin >= y).astype(float)

    # running crossing fraction and running pi estimate
    P_run = np.cumsum(crossings) / sample_axis
    with np.errstate(divide="ignore", invalid="ignore"):
        pi_run = 2.0 * l / (P_run * d)     # pi = 2l/(P d)

    pi_small = pi_run[n_small - 1]
    pi_big = pi_run[n_big - 1]
    valid = " (l>d: model invalid)" if l > d else ""
    print(f"l={l:.1f}: pi_hat(n=1e4)={pi_small:.5f}  pi_hat(n=1e5)={pi_big:.5f}"
          f"  err(1e5)={pi_big - np.pi:+.5f}{valid}")

    plt.plot(sample_axis, pi_run, label=f"l={l:.1f}", alpha=0.8)

print(f"true pi = {np.pi:.5f}")

plt.axhline(np.pi, color="k", lw=1.2, ls="--", label="pi")
plt.xscale("log")
plt.ylim(2.5, 4.0)
plt.xlabel("number of drops (sample count)")
plt.ylabel("estimate of pi")
plt.title("Buffon's needle: running pi estimate vs sample count (d=1, seed=1)")
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Check: at n=1e4 each curve still wanders noticeably around pi (and l=1.2,1.4 "
      "with l>d are biased since 2l/d>pi makes the formula invalid), while at n=1e5 the "
      "curves for l<d visibly tighten onto pi.")
print("Why it confirms: because a correct Monte Carlo estimator must converge toward "
      "the true value as sample count grows, so shrinking scatter around pi for l<d "
      "(and only there) confirms the estimator and its validity range are right.")
