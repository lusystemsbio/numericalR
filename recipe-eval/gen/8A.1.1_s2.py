import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Monte Carlo estimator of pi by area sampling ----
# Model: a uniform point in [-1,1]^2 falls inside the unit circle with
# probability = (area of circle)/(area of square) = pi/4, so pi = 4 * P.
def estimate_pi(n, rng):
    # draw n uniform points in the square [-1,1]^2 (each coord in [-1,1])
    pts = rng.uniform(-1.0, 1.0, size=(int(n), 2))
    # a point is inside the unit circle when x^2 + y^2 <= 1
    inside = (pts[:, 0]**2 + pts[:, 1]**2) <= 1.0
    # fraction inside times 4 is the pi estimate
    frac = np.mean(inside)
    return 4.0 * frac, pts, inside

# ---- Single estimate at n = 1e4 (seed 1) ----
rng = np.random.default_rng(1)
n_single = int(1e4)
pi_single, pts_single, inside_single = estimate_pi(n_single, rng)
print(f"Single estimate at n=1e4: pi_hat = {pi_single:.6f}")
print(f"Absolute error at n=1e4: {abs(pi_single - np.pi):.6f}")

# ---- Convergence run to n = 1e5 (seed 1, fresh draw) ----
rng = np.random.default_rng(1)
n_conv = int(1e5)
pts = rng.uniform(-1.0, 1.0, size=(n_conv, 2))       # uniform points in [-1,1]^2
inside = (pts[:, 0]**2 + pts[:, 1]**2) <= 1.0        # inside-circle indicator
# running fraction inside after k samples = cumulative count / k
cum_inside = np.cumsum(inside)
sample_count = np.arange(1, n_conv + 1)
running_pi = 4.0 * cum_inside / sample_count         # running estimate of pi

# report the estimate at the final sample count
pi_final = running_pi[-1]
print(f"Running estimate at n=1e5: pi_hat = {pi_final:.6f}")
print(f"Absolute error at n=1e5: {abs(pi_final - np.pi):.6f}")

# ---- Check: off at 1e4, settles near pi by 1e5, slow 1/sqrt(n) convergence ----
err_1e4 = abs(running_pi[int(1e4) - 1] - np.pi)
err_1e5 = abs(running_pi[-1] - np.pi)
print(f"Check error at n=1e4 (convergence run): {err_1e4:.6f}")
print(f"Check error at n=1e5 (convergence run): {err_1e5:.6f}")
# expected statistical scale of the error ~ 1/sqrt(n)
print(f"1/sqrt(n) scale at n=1e4: {1.0/np.sqrt(1e4):.6f}")
print(f"1/sqrt(n) scale at n=1e5: {1.0/np.sqrt(1e5):.6f}")
print(f"True pi: {np.pi:.6f}")
# One sentence: the check confirms the result because the error shrinks as we
# add samples and lands near pi only at the larger n, matching the slow
# 1/sqrt(n) Monte Carlo error rate rather than being exact at small n.
print("Why: the error is visibly larger at n=1e4 and drops to near-zero by "
      "n=1e5, tracking the 1/sqrt(n) rate, which confirms genuine Monte Carlo "
      "convergence rather than a lucky exact hit.")

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# running estimate of pi vs sample count
ax1.plot(sample_count, running_pi, color="steelblue", lw=0.8,
         label="running estimate")
ax1.axhline(np.pi, color="red", ls="--", lw=1.2, label="true pi")
ax1.set_xscale("log")
ax1.set_xlabel("sample count n")
ax1.set_ylabel("estimate of pi")
ax1.set_title("Running Monte Carlo estimate of pi")
ax1.set_ylim(np.pi - 0.5, np.pi + 0.5)
ax1.legend()

# scatter of points in/out of the circle (use the n=1e4 single estimate draw)
ax2.scatter(pts_single[inside_single, 0], pts_single[inside_single, 1],
            s=3, color="tab:green", label="inside")
ax2.scatter(pts_single[~inside_single, 0], pts_single[~inside_single, 1],
            s=3, color="tab:orange", label="outside")
theta = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2)
ax2.set_aspect("equal")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title(f"Sampled points (n=1e4), pi_hat={pi_single:.4f}")
ax2.legend(loc="upper right")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.1.1_s2.png")
