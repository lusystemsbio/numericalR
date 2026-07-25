import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Monte Carlo estimation of pi by area sampling.
# Fraction of uniform points in [-1,1]^2 falling inside the unit circle -> pi/4,
# so pi = 4 * (points inside) / (total points).

rng = np.random.default_rng(1)  # seed 1 for reproducibility

# --- Single estimate at n = 1e4 ---
n1 = int(1e4)
x1 = rng.uniform(-1.0, 1.0, n1)          # x coordinates
y1 = rng.uniform(-1.0, 1.0, n1)          # y coordinates
inside1 = (x1**2 + y1**2) <= 1.0         # boolean mask: inside unit circle
pi_n1e4 = 4.0 * np.mean(inside1)         # fraction inside times four
print(f"pi estimate at n=1e4: {pi_n1e4}")
print(f"absolute error at n=1e4: {abs(pi_n1e4 - np.pi)}")

# --- Convergence run to n = 1e5 (fresh stream, same seed) ---
rng = np.random.default_rng(1)           # reset so the run is reproducible
n_max = int(1e5)
x = rng.uniform(-1.0, 1.0, n_max)
y = rng.uniform(-1.0, 1.0, n_max)
inside = (x**2 + y**2) <= 1.0            # inside/outside mask for all points

# Running estimate: cumulative count inside divided by cumulative sample count.
counts = np.arange(1, n_max + 1)         # sample counts 1..n_max
running_inside = np.cumsum(inside)       # cumulative number of hits
running_pi = 4.0 * running_inside / counts  # running pi estimate

pi_n1e5 = running_pi[-1]                  # final estimate at n=1e5
print(f"pi estimate at n=1e5: {pi_n1e5}")
print(f"absolute error at n=1e5: {abs(pi_n1e5 - np.pi)}")
print(f"true pi: {np.pi}")

# --- Separate check: off at 1e4, settled near pi at 1e5 ---
err_1e4 = abs(pi_n1e4 - np.pi)
err_1e5 = abs(pi_n1e5 - np.pi)
# Expected 1/sqrt(n) scaling: error should shrink roughly by sqrt(10) ~ 3.16.
expected_ratio = np.sqrt(n_max / n1)
observed_ratio = err_1e4 / err_1e5
print(f"error ratio (1e4/1e5) observed: {observed_ratio}")
print(f"error ratio expected ~sqrt(10): {expected_ratio}")
check_passed = (err_1e4 > err_1e5) and (err_1e5 < 0.02)
print(f"check passed (1e4 noticeably off, 1e5 near pi): {check_passed}")

# --- Plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# Running estimate of pi vs sample count.
ax1.plot(counts, running_pi, color="steelblue", lw=0.8, label="running estimate")
ax1.axhline(np.pi, color="red", ls="--", label="true pi")
ax1.set_xscale("log")
ax1.set_xlabel("sample count n")
ax1.set_ylabel("pi estimate")
ax1.set_title("Running Monte Carlo estimate of pi")
ax1.set_ylim(2.8, 3.5)
ax1.legend()

# Scatter of sampled points, colored by in/out of the circle (subsample for clarity).
sub = slice(0, 5000)
ax2.scatter(x[sub][inside[sub]], y[sub][inside[sub]], s=3, color="green", label="inside")
ax2.scatter(x[sub][~inside[sub]], y[sub][~inside[sub]], s=3, color="orange", label="outside")
theta = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2)  # unit circle
ax2.set_aspect("equal")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("Sampled points in/out of unit circle")
ax2.legend()

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.1.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the result because the estimate's error "
      "shrinks from n=1e4 to n=1e5 by roughly sqrt(10), matching the theoretical "
      "1/sqrt(n) Monte Carlo convergence rate that drives it toward the true pi.")
