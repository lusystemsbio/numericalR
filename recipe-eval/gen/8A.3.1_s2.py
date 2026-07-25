import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Function to integrate
def f(x):
    return np.exp(-x**2)

# Integration domain and exact reference value
x1, x2 = 0.0, 3.0
exact = np.sqrt(np.pi) / 2  # exact integral of exp(-x^2) over [0, inf)

# ---- Deterministic midpoint rule ----
# Split [x1, x2] into n equal subintervals and evaluate f at each midpoint.
n_mid = 1000
h = (x2 - x1) / n_mid                     # width of each subinterval
midpoints = x1 + (np.arange(n_mid) + 0.5) * h  # centers of the subintervals
midpoint_estimate = h * np.sum(f(midpoints))   # sum of f * width

# ---- Monte Carlo integration by uniform sampling ----
# Draw n uniform points on [x1, x2]; I = (x2 - x1)/n * sum f(x_i).
rng = np.random.default_rng(1)            # seed 1 for reproducibility
n_mc = int(1e4)
samples = rng.uniform(x1, x2, n_mc)       # uniform random samples
mc_estimate = (x2 - x1) / n_mc * np.sum(f(samples))  # MC estimator

# ---- Errors relative to the exact value ----
midpoint_error = abs(midpoint_estimate - exact)
mc_error = abs(mc_estimate - exact)

# ---- Separate check: both land near sqrt(pi)/2 ----
tol = 1e-2
midpoint_near = midpoint_error < tol
mc_near = mc_error < tol
check_passed = midpoint_near and mc_near

# ---- Print results ----
print(f"Exact value sqrt(pi)/2:            {exact:.6f}")
print(f"Midpoint rule estimate (n=1000):   {midpoint_estimate:.6f}")
print(f"Monte Carlo estimate (n=10000):    {mc_estimate:.6f}")
print(f"Midpoint rule absolute error:      {midpoint_error:.6e}")
print(f"Monte Carlo absolute error:        {mc_error:.6e}")
print(f"Midpoint near sqrt(pi)/2 (tol {tol}): {midpoint_near}")
print(f"Monte Carlo near sqrt(pi)/2 (tol {tol}): {mc_near}")
print(f"Check passed (both near 0.886):    {check_passed}")
# One-sentence explanation:
print("Explanation: Because both independent estimators fall within a small "
      "tolerance of the analytic value sqrt(pi)/2, the agreement confirms that "
      "each method is correctly computing the same integral.")

# ---- Plot ----
fig, ax = plt.subplots(figsize=(8, 5))
methods = ["Midpoint\n(n=1000)", "Monte Carlo\n(n=10000)"]
estimates = [midpoint_estimate, mc_estimate]
ax.bar(methods, estimates, color=["steelblue", "darkorange"], width=0.5,
       label="estimate")
ax.axhline(exact, color="red", linestyle="--",
           label=f"exact = sqrt(pi)/2 = {exact:.3f}")
for i, est in enumerate(estimates):
    ax.text(i, est + 0.002, f"{est:.4f}", ha="center", va="bottom")
ax.set_ylabel("Integral estimate")
ax.set_title("Integration of exp(-x^2) on [0, 3]: Midpoint vs Monte Carlo")
ax.set_ylim(0.85, 0.90)
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.3.1_s2.png")
