import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model: I = integral over [0, 100] of f(x) = exp(-x^2) dx ----
# The true value is very close to the Gaussian integral sqrt(pi)/2,
# since exp(-x^2) is negligible beyond a few units.
from math import sqrt, pi, erf
a, b = 0.0, 100.0
n = int(1e4)          # samples per estimate
reps = 100            # replicates per method
true_val = (sqrt(pi) / 2.0) * erf(b)   # exact integral of exp(-x^2) on [0, b]

def f(x):
    return np.exp(-x**2)

# Seed once for reproducibility of the whole experiment.
rng = np.random.default_rng(1)

uniform_estimates = np.empty(reps)
importance_estimates = np.empty(reps)

for r in range(reps):
    # ---- Uniform Monte Carlo ----
    # Draw x ~ Uniform(a, b); here p(x) = 1/(b-a) is constant.
    xu = rng.uniform(a, b, size=n)
    # Estimator is the mean of f(x)/p(x) = f(x)*(b-a).
    uniform_estimates[r] = np.mean(f(xu) * (b - a))

    # ---- Importance sampling with p(x) = exp(-x) ----
    # Inverse-transform sampling: if u ~ Uniform(0,1), then
    # x = -ln(1-u) is distributed as Exp(1) with density p(x)=exp(-x).
    u = rng.uniform(0.0, 1.0, size=n)
    xi = -np.log(1.0 - u)
    # Estimator is the mean of the importance weight f(x)/p(x).
    # exp(-x) has support on [0, inf), matching the tail we care about;
    # points beyond b=100 contribute ~0 to f so truncation is irrelevant.
    pi_dens = np.exp(-xi)
    importance_estimates[r] = np.mean(f(xi) / pi_dens)

# ---- Report numerical results ----
print(f"True integral value:                {true_val:.10f}")
print(f"Uniform MC   mean of estimates:     {np.mean(uniform_estimates):.10f}")
print(f"Uniform MC   std  of estimates:     {np.std(uniform_estimates, ddof=1):.10f}")
print(f"Uniform MC   variance of estimates: {np.var(uniform_estimates, ddof=1):.10e}")
print(f"Importance   mean of estimates:     {np.mean(importance_estimates):.10f}")
print(f"Importance   std  of estimates:     {np.std(importance_estimates, ddof=1):.10f}")
print(f"Importance   variance of estimates: {np.var(importance_estimates, ddof=1):.10e}")
print(f"Variance reduction factor (uniform/importance): "
      f"{np.var(uniform_estimates, ddof=1) / np.var(importance_estimates, ddof=1):.4f}")

# ---- Separate check: where do the points land relative to f ? ----
# f(x) = exp(-x^2) is essentially zero for x > ~4. Count how many of the
# n sampled points fall in the "useful" region [0, 4] for one replicate.
xu_check = rng.uniform(a, b, size=n)
u_check = rng.uniform(0.0, 1.0, size=n)
xi_check = -np.log(1.0 - u_check)
useful = 4.0
frac_uniform_useful = np.mean(xu_check <= useful)
frac_importance_useful = np.mean(xi_check <= useful)
print(f"Fraction of uniform points in [0,{useful}] (where f matters):    "
      f"{frac_uniform_useful:.4f}")
print(f"Fraction of importance points in [0,{useful}] (where f matters): "
      f"{frac_importance_useful:.4f}")

# ---- Box plots of the two sets of estimates ----
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([importance_estimates, uniform_estimates],
           labels=["Importance\n(p=exp(-x))", "Uniform\nU(0,100)"],
           showmeans=True)
ax.axhline(true_val, color="red", linestyle="--", label=f"true = {true_val:.5f}")
ax.set_ylabel("Estimate of the integral")
ax.set_title("Monte Carlo estimates over 100 replicates (n=1e4, range [0,100])")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.4.1_s5.png")

# ---- One-sentence explanation ----
print("Explanation: The check confirms the result because almost all "
      "importance-sampling points land in [0,4] where f is large while "
      "almost no uniform points do, so the exponential density concentrates "
      "effort where the integrand contributes and yields far tighter, "
      "less-variable estimates than uniform sampling.")
