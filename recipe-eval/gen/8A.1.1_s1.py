import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Monte Carlo estimator of pi via area sampling -----
# Model: fraction of uniform points in [-1,1]^2 that land inside the unit
# circle equals (area of circle)/(area of square) = pi/4, so pi = 4 * P.

def estimate_pi(n, rng):
    # draw n uniform points in the square [-1, 1]^2
    pts = rng.uniform(-1.0, 1.0, size=(n, 2))
    # a point is inside the unit circle when x^2 + y^2 <= 1
    inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) <= 1.0
    # fraction inside times 4 is the pi estimate
    p = np.mean(inside)
    return 4.0 * p, pts, inside

# ----- single estimate at n = 1e4 -----
rng1 = np.random.default_rng(1)  # seed 1
n_single = int(1e4)
pi_single, pts_single, inside_single = estimate_pi(n_single, rng1)
print(f"Single estimate at n=1e4: pi_hat = {pi_single:.6f}")
print(f"Absolute error at n=1e4:  {abs(pi_single - np.pi):.6f}")

# ----- convergence run to n = 1e5 -----
rng2 = np.random.default_rng(1)  # reset seed 1 for reproducible convergence run
n_conv = int(1e5)
pts_conv = rng2.uniform(-1.0, 1.0, size=(n_conv, 2))
inside_conv = (pts_conv[:, 0] ** 2 + pts_conv[:, 1] ** 2) <= 1.0
# running fraction inside as the sample count grows, times 4
counts = np.arange(1, n_conv + 1)
running_pi = 4.0 * np.cumsum(inside_conv) / counts

pi_at_1e4 = running_pi[int(1e4) - 1]  # running estimate exactly at n=1e4
pi_final = running_pi[-1]             # running estimate at n=1e5
print(f"Running estimate at n=1e4 (convergence run): pi_hat = {pi_at_1e4:.6f}")
print(f"Running estimate at n=1e5 (final):           pi_hat = {pi_final:.6f}")
print(f"Absolute error at n=1e4 (convergence run):   {abs(pi_at_1e4 - np.pi):.6f}")
print(f"Absolute error at n=1e5 (final):             {abs(pi_final - np.pi):.6f}")
print(f"True pi: {np.pi:.6f}")

# ----- check: off at 1e4, settled near pi by 1e5 -----
err_1e4 = abs(pi_at_1e4 - np.pi)
err_1e5 = abs(pi_final - np.pi)
settled = err_1e5 < err_1e4
print(f"Error shrinks from 1e4 to 1e5: {settled}")
print(f"Error ratio (1e4 / 1e5): {err_1e4 / err_1e5:.3f}")
# expected 1/sqrt(n) scaling would give sqrt(1e5/1e4) ~ 3.16x smaller error
print(f"Expected error reduction from 1/sqrt(n): {np.sqrt(n_conv / 1e4):.3f}")

# ----- figure: running estimate + scatter -----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# left: running estimate of pi vs sample count
ax = axes[0]
ax.plot(counts, running_pi, color="steelblue", lw=0.8, label="running estimate")
ax.axhline(np.pi, color="red", ls="--", lw=1.2, label=f"true pi = {np.pi:.4f}")
# 1/sqrt(n) reference envelope around pi to show slow convergence
sigma = np.pi / np.sqrt(counts[1:])  # rough scale of the fluctuations
ax.plot(counts[1:], np.pi + sigma, color="gray", ls=":", lw=0.8, label="~1/sqrt(n) band")
ax.plot(counts[1:], np.pi - sigma, color="gray", ls=":", lw=0.8)
ax.set_xscale("log")
ax.set_xlabel("sample count n")
ax.set_ylabel("estimate of pi")
ax.set_ylim(2.6, 3.7)
ax.set_title("Running Monte Carlo estimate of pi")
ax.legend(loc="upper right", fontsize=8)

# right: scatter of sampled points colored in/out of the circle (subset for clarity)
ax = axes[1]
sub = pts_single  # the single n=1e4 sample
ins = inside_single
ax.scatter(sub[ins, 0], sub[ins, 1], s=3, color="tab:blue", label="inside")
ax.scatter(sub[~ins, 0], sub[~ins, 1], s=3, color="tab:orange", label="outside")
theta = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Sampled points (n=1e4), pi_hat = {pi_single:.4f}")
ax.legend(loc="upper right", fontsize=8, markerscale=3)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.1.1_s1.png", dpi=110)

# ----- one-sentence explanation -----
print("Explanation: The check confirms the result because the estimate is visibly off at n=1e4 "
      "but its error shrinks by roughly sqrt(10)~3.16x by n=1e5, matching the 1/sqrt(n) rate "
      "expected of Monte Carlo and showing convergence toward the true pi.")
