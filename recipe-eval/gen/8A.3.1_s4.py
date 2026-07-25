import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# The integrand: f(x) = exp(-x^2)
def f(x):
    return np.exp(-x**2)

# Integration interval and parameters
x1, x2 = 0.0, 3.0        # finite interval used for both methods
n_mid = 1000             # number of midpoint subintervals
n_mc = int(1e4)          # number of Monte Carlo samples
np.random.seed(1)        # fixed seed for reproducibility

# Exact value of the integral of exp(-x^2) over [0, inf)
exact = np.sqrt(np.pi) / 2.0

# --- Deterministic midpoint rule ---
# Split [x1, x2] into n_mid equal panels of width h, evaluate f at each
# panel center, and sum: integral ~ h * sum(f(midpoints)).
h = (x2 - x1) / n_mid                       # width of each subinterval
midpoints = x1 + (np.arange(n_mid) + 0.5) * h  # centers of the subintervals
midpoint_estimate = h * np.sum(f(midpoints))

# --- Monte Carlo integration by uniform sampling ---
# Draw n_mc uniform points on [x1, x2]; the estimator is the interval
# length times the mean of f: I = (x2 - x1)/n * sum f(x_i).
xs = np.random.uniform(x1, x2, size=n_mc)   # uniform samples on [x1, x2]
mc_estimate = (x2 - x1) / n_mc * np.sum(f(xs))

# Errors relative to the exact value
mid_error = abs(midpoint_estimate - exact)
mc_error = abs(mc_estimate - exact)

# --- Print numerical results ---
print(f"Exact value sqrt(pi)/2:            {exact:.6f}")
print(f"Midpoint-rule estimate (n=1000):   {midpoint_estimate:.6f}")
print(f"Monte Carlo estimate (n=10000):    {mc_estimate:.6f}")
print(f"Midpoint-rule absolute error:      {mid_error:.6e}")
print(f"Monte Carlo absolute error:        {mc_error:.6e}")

# --- Separate check: both land near sqrt(pi)/2 = 0.886 ---
tol = 0.05  # loose tolerance; MC has larger, dimension-independent error
mid_ok = mid_error < tol
mc_ok = mc_error < tol
print(f"Midpoint estimate near 0.886 (tol={tol}): {mid_ok}")
print(f"Monte Carlo estimate near 0.886 (tol={tol}): {mc_ok}")
print(f"Both estimates near sqrt(pi)/2:    {mid_ok and mc_ok}")
# This check confirms the result because both independent methods converging
# to the same known analytic value sqrt(pi)/2 shows the numerics are correct
# rather than accidentally agreeing on a wrong number.

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
labels = ["Midpoint\n(n=1000)", "Monte Carlo\n(n=10000)"]
estimates = [midpoint_estimate, mc_estimate]
ax.bar(labels, estimates, color=["steelblue", "indianred"], alpha=0.8)
ax.axhline(exact, color="black", linestyle="--", label=f"Exact = {exact:.4f}")
ax.set_ylabel("Integral estimate")
ax.set_title("Integral of exp(-x^2): estimates vs sqrt(pi)/2")
for i, v in enumerate(estimates):
    ax.text(i, v + 0.005, f"{v:.4f}", ha="center")
ax.set_ylim(0.7, 1.0)
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.3.1_s4.png")
