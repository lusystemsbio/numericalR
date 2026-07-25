import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Model: I = integral_0^100 exp(-x^2) dx  (effectively the Gaussian tail
# integral; true value ~ sqrt(pi)/2 since the range is very wide)
# ----------------------------------------------------------------------
def f(x):
    return np.exp(-x**2)

a, b = 0.0, 100.0          # integration range
n = int(1e4)               # samples per replicate
reps = 100                 # replicates for each method
true_val = np.sqrt(np.pi) / 2.0   # reference value of int_0^inf exp(-x^2) dx

rng = np.random.default_rng(1)     # seed 1

uniform_estimates = np.empty(reps)
importance_estimates = np.empty(reps)

for r in range(reps):
    # --- Uniform Monte Carlo -----------------------------------------
    # Draw x uniformly on [a,b]; density is 1/(b-a), so the estimator is
    # (b-a) * mean(f(x)).  Most points land where f is negligible.
    xu = rng.uniform(a, b, n)
    uniform_estimates[r] = (b - a) * np.mean(f(xu))

    # --- Importance sampling with p(x) = exp(-x) ---------------------
    # Inverse-transform sampling: if u ~ Uniform(0,1) then x = -ln(1-u)
    # is Exponential(1), concentrated near 0 where f is large.
    u = rng.uniform(0.0, 1.0, n)
    xi = -np.log(1.0 - u)                 # draw from p by inverse transform
    p = np.exp(-xi)                        # proposal density p(x)=exp(-x)
    # Only keep points inside [a,b] (exponential support is [0,inf))
    inside = (xi >= a) & (xi <= b)
    weights = np.where(inside, f(xi) / p, 0.0)   # f(x)/p(x)
    importance_estimates[r] = np.mean(weights)    # mean of f/p over draws

# ----------------------------------------------------------------------
# Report numerical results
# ----------------------------------------------------------------------
print(f"True value (sqrt(pi)/2): {true_val:.10f}")
print(f"Uniform     mean estimate: {np.mean(uniform_estimates):.10f}")
print(f"Uniform     std  estimate: {np.std(uniform_estimates, ddof=1):.10e}")
print(f"Uniform     variance      : {np.var(uniform_estimates, ddof=1):.10e}")
print(f"Importance  mean estimate: {np.mean(importance_estimates):.10f}")
print(f"Importance  std  estimate: {np.std(importance_estimates, ddof=1):.10e}")
print(f"Importance  variance      : {np.var(importance_estimates, ddof=1):.10e}")
print(f"Variance reduction factor (uniform/importance): "
      f"{np.var(uniform_estimates, ddof=1)/np.var(importance_estimates, ddof=1):.4f}")

# --- Separate check: how many sampled points fall where f is non-negligible
thresh = 1e-6
u_useful = np.mean(f(rng.uniform(a, b, n)) > thresh)
xi_check = -np.log(1.0 - rng.uniform(0.0, 1.0, n))
i_useful = np.mean(f(xi_check[(xi_check >= a) & (xi_check <= b)]) > thresh)
print(f"Fraction of uniform    points with f(x)>{thresh:g}: {u_useful:.4f}")
print(f"Fraction of importance points with f(x)>{thresh:g}: {i_useful:.4f}")

# ----------------------------------------------------------------------
# Box plots
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([importance_estimates, uniform_estimates],
           labels=["Importance sampling", "Uniform MC"])
ax.axhline(true_val, color="red", linestyle="--", label=f"true = {true_val:.4f}")
ax.set_ylabel("Estimate of the integral")
ax.set_title("Importance sampling vs uniform Monte Carlo (100 replicates)")
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.4.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the result because uniform sampling "
      "wastes almost all its points on the vast region where f is essentially "
      "zero (a tiny useful fraction, hence noisy estimates), whereas the "
      "exponential proposal places nearly all points where f is large, giving "
      "the far smaller variance seen in the box plots.")
