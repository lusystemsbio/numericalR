import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Integrand
def f(x):
    return np.exp(-x**2)

# Integration settings
x1, x2 = 0.0, 3.0          # finite interval approximating [0, inf)
exact = np.sqrt(np.pi) / 2  # exact value of integral over [0, inf)

# --- Deterministic midpoint rule (implemented explicitly) ---
n_mid = 1000
h = (x2 - x1) / n_mid                       # width of each subinterval
midpoints = x1 + (np.arange(n_mid) + 0.5) * h  # centers of subintervals
midpoint_est = h * np.sum(f(midpoints))     # sum of f at midpoints times width

# --- Monte Carlo integration by uniform sampling (implemented explicitly) ---
rng = np.random.default_rng(1)              # seed 1 for reproducibility
n_mc = int(1e4)
x_samples = rng.uniform(x1, x2, size=n_mc)  # uniform samples in [x1, x2]
# I = (x2 - x1)/n * sum f(x_i)  --> average of f times interval length
mc_est = (x2 - x1) / n_mc * np.sum(f(x_samples))

# Errors relative to the exact value
mid_err = abs(midpoint_est - exact)
mc_err = abs(mc_est - exact)

# --- Report ---
print(f"Exact value sqrt(pi)/2:              {exact:.6f}")
print(f"Midpoint rule estimate (n={n_mid}):    {midpoint_est:.6f}")
print(f"Monte Carlo estimate (n={n_mc}):      {mc_est:.6f}")
print(f"Midpoint rule absolute error:        {mid_err:.6e}")
print(f"Monte Carlo absolute error:          {mc_err:.6e}")

# --- Separate check: both should land near sqrt(pi)/2 = 0.886 ---
tol = 0.05
mid_ok = mid_err < tol
mc_ok = mc_err < tol
print(f"Midpoint within {tol} of exact:        {mid_ok}")
print(f"Monte Carlo within {tol} of exact:     {mc_ok}")
print(f"Both near sqrt(pi)/2 = 0.886:        {mid_ok and mc_ok}")
# Note: midpoint rule is more accurate for this smooth 1D integrand,
# while Monte Carlo's error scales as 1/sqrt(n) independent of dimension.

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
labels = ["Midpoint\n(n=1000)", "Monte Carlo\n(n=1e4)"]
estimates = [midpoint_est, mc_est]
ax.bar(labels, estimates, color=["steelblue", "darkorange"], alpha=0.8)
ax.axhline(exact, color="red", linestyle="--", label=f"exact = {exact:.4f}")
ax.set_ylabel("Integral estimate")
ax.set_title("Integration of exp(-x^2) on [0, 3]")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.3.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: because both independent methods converge to the "
      "same known value sqrt(pi)/2 = 0.886, agreement confirms each method "
      "is correctly implemented rather than sharing a coincidental error.")
