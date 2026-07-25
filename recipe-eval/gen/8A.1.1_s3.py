import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Monte Carlo estimator of pi by area sampling ---------------------------
# Model: fraction of uniform points in [-1,1]^2 that land inside the unit
# circle equals (area of circle)/(area of square) = pi/4, so pi = 4 * P.

rng = np.random.default_rng(1)  # seed 1 for reproducibility


def estimate_pi(n, rng):
    # Draw n uniform points in the square [-1,1]^2 (2 coords each).
    pts = rng.uniform(-1.0, 1.0, size=(n, 2))
    # A point is inside the unit circle if x^2 + y^2 <= 1.
    inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) <= 1.0
    # Fraction inside times 4 estimates pi.
    p = inside.mean()
    return 4.0 * p, pts, inside


# --- Single estimate at n = 1e4 ---------------------------------------------
n_single = int(1e4)
pi_single, pts_single, inside_single = estimate_pi(n_single, np.random.default_rng(1))
print(f"Single estimate at n=1e4:            pi ~= {pi_single:.6f}")
print(f"Absolute error at n=1e4:             {abs(pi_single - np.pi):.6f}")

# --- Convergence run to n = 1e5 ---------------------------------------------
# Reuse one stream of points and accumulate the running estimate so we can
# watch it converge as the sample count grows.
n_max = int(1e5)
pts = np.random.default_rng(1).uniform(-1.0, 1.0, size=(n_max, 2))
inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) <= 1.0  # boolean hit/miss per point
counts = np.arange(1, n_max + 1)                   # 1,2,...,n_max
running_p = np.cumsum(inside) / counts             # running fraction inside
running_pi = 4.0 * running_p                       # running estimate of pi

pi_final = running_pi[-1]
print(f"Running estimate at n=1e4:           pi ~= {running_pi[n_single - 1]:.6f}")
print(f"Final estimate at n=1e5:             pi ~= {pi_final:.6f}")
print(f"Absolute error at n=1e5:             {abs(pi_final - np.pi):.6f}")
print(f"True value of pi:                     {np.pi:.6f}")

# --- Separate check: off at 1e4, settles near pi by 1e5 ---------------------
err_1e4 = abs(running_pi[n_single - 1] - np.pi)
err_1e5 = abs(pi_final - np.pi)
print(f"Error shrank from n=1e4 to n=1e5:    {err_1e4:.6f} -> {err_1e5:.6f}")
# The 1/sqrt(n) rate predicts error should drop by ~sqrt(10) ~= 3.16x.
print(f"Error ratio (1e4 err / 1e5 err):     {err_1e4 / err_1e5:.3f}")
print(f"Expected ~sqrt(10) drop:             {np.sqrt(10):.3f}")

# --- Plots ------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# Left: running estimate of pi vs sample count.
ax1.plot(counts, running_pi, lw=0.8, color="steelblue", label="running estimate")
ax1.axhline(np.pi, color="crimson", ls="--", lw=1.2, label="true pi")
ax1.axvline(n_single, color="gray", ls=":", lw=1.0, label="n=1e4")
ax1.set_xscale("log")
ax1.set_xlabel("sample count n")
ax1.set_ylabel("estimate of pi")
ax1.set_ylim(np.pi - 0.4, np.pi + 0.4)
ax1.set_title("Running Monte Carlo estimate of pi")
ax1.legend()

# Right: scatter of the single-estimate points, colored in/out of circle.
ax2.scatter(pts_single[inside_single, 0], pts_single[inside_single, 1],
            s=2, color="seagreen", label="inside")
ax2.scatter(pts_single[~inside_single, 0], pts_single[~inside_single, 1],
            s=2, color="salmon", label="outside")
theta = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2)
ax2.set_aspect("equal")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title(f"Sampled points (n={n_single})")
ax2.legend(markerscale=4)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.1.1_s3.png")

# One-sentence explanation of why the check confirms the result.
print("Explanation: seeing the error shrink by roughly sqrt(10) as n grows "
      "from 1e4 to 1e5 confirms the estimator converges to pi at the "
      "characteristic 1/sqrt(n) Monte Carlo rate rather than by luck.")
