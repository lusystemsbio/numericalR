import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Function to integrate
def f(x):
    return np.exp(-x**2)

# Integration parameters
x1, x2 = 0.0, 3.0          # integration interval [0, 3]
exact = np.sqrt(np.pi) / 2 # exact value of integral over [0, inf)

# --- Midpoint rule (deterministic) ---
n_mid = 1000
h = (x2 - x1) / n_mid                       # width of each subinterval
midpoints = x1 + (np.arange(n_mid) + 0.5) * h  # centers of each subinterval
midpoint_estimate = h * np.sum(f(midpoints))   # sum of area of each rectangle

# --- Monte Carlo integration (uniform sampling) ---
n_mc = int(1e4)
rng = np.random.default_rng(1)               # seed 1 for reproducibility
x_samples = rng.uniform(x1, x2, size=n_mc)   # uniform samples in [x1, x2]
mc_estimate = (x2 - x1) / n_mc * np.sum(f(x_samples))  # average value * width

# --- Errors ---
midpoint_error = abs(midpoint_estimate - exact)
mc_error = abs(mc_estimate - exact)

# --- Separate check: both land near sqrt(pi)/2 ---
tol = 0.05
check_midpoint = midpoint_error < tol
check_mc = mc_error < tol
both_near = check_midpoint and check_mc

# --- Print results ---
print(f"Exact value sqrt(pi)/2:        {exact:.6f}")
print(f"Midpoint estimate (n={n_mid}):    {midpoint_estimate:.6f}")
print(f"Monte Carlo estimate (n={n_mc}): {mc_estimate:.6f}")
print(f"Midpoint absolute error:       {midpoint_error:.6e}")
print(f"Monte Carlo absolute error:    {mc_error:.6e}")
print(f"Midpoint within {tol} of exact: {check_midpoint}")
print(f"Monte Carlo within {tol} of exact: {check_mc}")
print(f"Both estimates near sqrt(pi)/2: {both_near}")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: Because both independent methods (a deterministic "
      "rule and a random sampler) converge to sqrt(pi)/2 = 0.886, agreement "
      "confirms the estimate is correct rather than an artifact of one method.")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
methods = ["Midpoint\n(n=1000)", "Monte Carlo\n(n=10000)"]
estimates = [midpoint_estimate, mc_estimate]
ax.bar(methods, estimates, color=["steelblue", "indianred"], alpha=0.8, width=0.5)
ax.axhline(exact, color="black", linestyle="--", label=f"Exact = {exact:.4f}")
for i, est in enumerate(estimates):
    ax.text(i, est + 0.005, f"{est:.4f}", ha="center")
ax.set_ylabel("Integral estimate")
ax.set_title("Integrating exp(-x^2) on [0,3]: Midpoint vs Monte Carlo")
ax.set_ylim(0.80, 0.92)
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.3.1_s3.png")
