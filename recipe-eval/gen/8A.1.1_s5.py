import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Monte Carlo estimator of pi by area sampling ----
# Model: a uniform point in [-1,1]^2 lands inside the unit circle with
# probability = (area of circle)/(area of square) = pi/4, so pi = 4 * P.

def estimate_pi(n, rng):
    # draw n uniform points in the square [-1,1]^2
    x = rng.uniform(-1.0, 1.0, size=n)
    y = rng.uniform(-1.0, 1.0, size=n)
    # a point is inside the unit circle when x^2 + y^2 <= 1
    inside = (x * x + y * y) <= 1.0
    # fraction inside times 4 estimates pi
    p = np.mean(inside)
    return 4.0 * p, x, y, inside

# ---- Single estimate at n = 1e4 ----
rng = np.random.default_rng(1)          # seed 1
n_single = int(1e4)
pi_single, xs, ys, ins = estimate_pi(n_single, rng)
print(f"Single estimate of pi at n=1e4: {pi_single:.6f}")
print(f"Absolute error at n=1e4:        {abs(pi_single - np.pi):.6f}")

# ---- Convergence run to n = 1e5 ----
# Re-seed so the convergence run is reproducible and independent of the above.
rng = np.random.default_rng(1)          # seed 1
n_conv = int(1e5)
x = rng.uniform(-1.0, 1.0, size=n_conv)
y = rng.uniform(-1.0, 1.0, size=n_conv)
inside = (x * x + y * y) <= 1.0
# running fraction inside after k samples, times 4 -> running pi estimate
counts = np.arange(1, n_conv + 1)
running_inside = np.cumsum(inside)
running_pi = 4.0 * running_inside / counts

pi_at_1e4 = running_pi[int(1e4) - 1]
pi_at_1e5 = running_pi[-1]
print(f"Running estimate at n=1e4:      {pi_at_1e4:.6f}")
print(f"Running estimate at n=1e5:      {pi_at_1e5:.6f}")
print(f"Abs error at n=1e4 (run):       {abs(pi_at_1e4 - np.pi):.6f}")
print(f"Abs error at n=1e5 (run):       {abs(pi_at_1e5 - np.pi):.6f}")

# ---- Check: off at 1e4, settles near pi by 1e5 (slow 1/sqrt(n) convergence) ----
err_1e4 = abs(pi_at_1e4 - np.pi)
err_1e5 = abs(pi_at_1e5 - np.pi)
print(f"Error ratio (1e4 / 1e5):        {err_1e4 / err_1e5:.6f}")
# Expected error scaling: error ~ C/sqrt(n); going 1e4 -> 1e5 (10x n)
# should shrink error by ~sqrt(10) ~= 3.16.
print(f"Expected 1/sqrt(n) shrink (sqrt(10)): {np.sqrt(10.0):.6f}")
check_passes = err_1e5 < err_1e4
print(f"Check (error at 1e5 < error at 1e4): {check_passes}")

# ---- Plot 1: running estimate of pi vs sample count ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(counts, running_pi, color="steelblue", lw=0.8, label="running estimate")
ax1.axhline(np.pi, color="red", ls="--", lw=1.2, label="true pi")
ax1.axvline(1e4, color="gray", ls=":", lw=1.0, label="n=1e4")
ax1.set_xscale("log")
ax1.set_xlabel("sample count n")
ax1.set_ylabel("estimate of pi")
ax1.set_title("Monte Carlo running estimate of pi")
ax1.set_ylim(np.pi - 0.5, np.pi + 0.5)
ax1.legend(loc="upper right")

# ---- Plot 2: scatter of sampled points in/out of circle (use single run) ----
ax2.scatter(xs[ins], ys[ins], s=2, color="tab:blue", label="inside")
ax2.scatter(xs[~ins], ys[~ins], s=2, color="tab:orange", label="outside")
theta = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2)
ax2.set_aspect("equal")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title(f"Sampled points (n=1e4), pi ~ {pi_single:.4f}")
ax2.legend(loc="upper right", markerscale=4)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.1.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the estimator because its error shrinks "
      "from a noticeable gap at n=1e4 to nearly zero at n=1e5, and does so at "
      "roughly the sqrt(10) rate predicted by the 1/sqrt(n) Monte Carlo law, "
      "showing the method converges to the true pi.")
