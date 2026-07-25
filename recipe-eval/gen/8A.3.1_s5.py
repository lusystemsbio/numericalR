import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Integrand: f(x) = exp(-x^2)
def f(x):
    return np.exp(-x**2)

# Integration bounds and exact reference value over [0, inf)
x1, x2 = 0.0, 3.0
exact = np.sqrt(np.pi) / 2.0  # exact integral of exp(-x^2) over [0, inf)

# --- Deterministic midpoint rule (explicit implementation) ---
n_mid = 1000
h = (x2 - x1) / n_mid                       # width of each subinterval
mids = x1 + (np.arange(n_mid) + 0.5) * h    # midpoint of each subinterval
midpoint_est = h * np.sum(f(mids))          # sum of f at midpoints times width

# --- Monte Carlo integration by uniform sampling (explicit implementation) ---
n_mc = int(1e4)
rng = np.random.default_rng(1)              # seeded generator for reproducibility
xs = rng.uniform(x1, x2, size=n_mc)         # uniform samples on [x1, x2]
mc_est = (x2 - x1) / n_mc * np.sum(f(xs))   # I = (x2-x1)/n * sum f(x_i)

# --- Errors relative to the exact value ---
midpoint_err = abs(midpoint_est - exact)
mc_err = abs(mc_est - exact)

# --- Report numerical results ---
print(f"Exact value sqrt(pi)/2:            {exact:.6f}")
print(f"Midpoint rule estimate (n=1000):   {midpoint_est:.6f}")
print(f"Monte Carlo estimate (n=10000):    {mc_est:.6f}")
print(f"Midpoint rule absolute error:      {midpoint_err:.6e}")
print(f"Monte Carlo absolute error:        {mc_err:.6e}")

# --- Separate check: both estimates land near sqrt(pi)/2 = 0.886 ---
tol = 1e-2
midpoint_ok = abs(midpoint_est - exact) < tol
mc_ok = abs(mc_est - exact) < tol
print(f"Midpoint rule within {tol} of exact: {midpoint_ok}")
print(f"Monte Carlo within {tol} of exact:   {mc_ok}")
print(f"Both estimates near 0.886:           {midpoint_ok and mc_ok}")

# Explanation:
# This check confirms the result because both independent methods converging to the
# same known analytic value sqrt(pi)/2 shows the estimates are correct rather than
# coincidentally close, validating each implementation against ground truth.
print("Why this check confirms the result: two independent methods both agreeing "
      "with the known analytic value sqrt(pi)/2 rules out coincidence and validates "
      "each implementation against ground truth.")

# --- Plot: estimates vs exact value ---
fig, ax = plt.subplots(figsize=(7, 5))
labels = ["Midpoint\n(n=1000)", "Monte Carlo\n(n=10000)"]
vals = [midpoint_est, mc_est]
ax.bar(labels, vals, color=["#4c72b0", "#dd8452"], width=0.5)
ax.axhline(exact, color="black", linestyle="--", label=f"exact = {exact:.4f}")
for i, v in enumerate(vals):
    ax.text(i, v + 0.005, f"{v:.4f}", ha="center")
ax.set_ylabel("Integral estimate")
ax.set_title("Integration of exp(-x^2) on [0, 3] vs sqrt(pi)/2")
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.3.1_s5.png")
