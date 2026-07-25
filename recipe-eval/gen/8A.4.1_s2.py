import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Target: I = integral_0^100 f(x) dx, with f(x) = exp(-x^2)
# Reference (essentially exact, since f is negligible past ~6):
#   integral_0^inf exp(-x^2) dx = sqrt(pi)/2
# ---------------------------------------------------------------
def f(x):
    return np.exp(-x**2)

a, b = 0.0, 100.0          # integration range
n = int(1e4)               # samples per estimate
reps = 100                 # replicates of each method
true_val = np.sqrt(np.pi) / 2.0
print(f"True value (sqrt(pi)/2): {true_val:.10f}")

rng = np.random.default_rng(1)   # seed 1

# storage for the estimates from each replicate
unif_ests = np.empty(reps)
imp_ests = np.empty(reps)

for r in range(reps):
    # --- Uniform Monte Carlo ---
    # draw x uniformly on [a,b]; density is 1/(b-a),
    # so I_hat = mean( f(x) / (1/(b-a)) ) = (b-a) * mean(f(x))
    xu = rng.uniform(a, b, n)
    unif_ests[r] = (b - a) * np.mean(f(xu))

    # --- Importance Sampling with p(x) = exp(-x) on [0, inf) ---
    # inverse transform: if U ~ Uniform(0,1), then x = -log(U) ~ Exp(1)
    u = rng.uniform(0.0, 1.0, n)
    xi = -np.log(u)                       # samples from p(x) = exp(-x)
    p = np.exp(-xi)                       # importance density values
    inside = (xi >= a) & (xi <= b)        # keep points in the target range
    w = np.where(inside, f(xi) / p, 0.0)  # weights f(x)/p(x), zero outside
    imp_ests[r] = np.mean(w)              # mean of f(x_i)/p(x_i)

# ---------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------
print(f"Uniform  MC: mean estimate = {np.mean(unif_ests):.10f}")
print(f"Uniform  MC: std  of estimates = {np.std(unif_ests, ddof=1):.10e}")
print(f"Uniform  MC: variance of estimates = {np.var(unif_ests, ddof=1):.10e}")
print(f"Importance : mean estimate = {np.mean(imp_ests):.10f}")
print(f"Importance : std  of estimates = {np.std(imp_ests, ddof=1):.10e}")
print(f"Importance : variance of estimates = {np.var(imp_ests, ddof=1):.10e}")

var_ratio = np.var(unif_ests, ddof=1) / np.var(imp_ests, ddof=1)
print(f"Variance reduction factor (uniform var / importance var): {var_ratio:.4f}")

# --- Separate check: how many drawn points land where f is non-negligible ---
# For one replicate, count points with f(x) > 1e-6 for each scheme.
xu_check = rng.uniform(a, b, n)
u_check = rng.uniform(0.0, 1.0, n)
xi_check = -np.log(u_check)
thresh = 1e-6
unif_useful = np.mean(f(xu_check) > thresh)
imp_useful = np.mean((xi_check <= b) & (f(xi_check) > thresh))
print(f"Fraction of uniform points where f>1e-6:    {unif_useful:.4f}")
print(f"Fraction of importance points where f>1e-6: {imp_useful:.4f}")

# ---------------------------------------------------------------
# Box plots of the two sets of estimates
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([imp_ests, unif_ests], labels=["Importance", "Uniform"],
           showmeans=True)
ax.axhline(true_val, color="red", linestyle="--", label=f"True = {true_val:.5f}")
ax.set_ylabel("Estimate of I")
ax.set_title("Importance sampling vs uniform Monte Carlo (100 replicates)")
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.4.1_s2.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# The near-zero fraction of useful uniform points versus the large useful
# fraction of importance points shows uniform sampling wastes almost all
# draws where f is negligible while importance sampling concentrates draws
# where f is large, which is exactly why its estimate variance is far smaller.
# ---------------------------------------------------------------
print("Check: uniform wastes points where f~0 (tiny useful fraction) while "
      "importance concentrates them where f is large, explaining its much "
      "lower variance.")
