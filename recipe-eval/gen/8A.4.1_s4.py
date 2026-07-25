import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Problem setup -----------------------------------------------------
# Integrand: f(x) = exp(-x^2), integrated over the wide range [0, 100].
# True value over [0, inf) is sqrt(pi)/2; on [0,100] it is effectively the same
# because f is negligible for x beyond a few units.
f = lambda x: np.exp(-x**2)
a, b = 0.0, 100.0
n = int(1e4)          # samples per estimate
reps = 100            # replicates of each method
true_val = np.sqrt(np.pi) / 2.0
np.random.seed(1)     # reproducibility

# Importance density p(x) = exp(-x) on [0, inf).
p = lambda x: np.exp(-x)

uniform_estimates = np.empty(reps)
importance_estimates = np.empty(reps)

for r in range(reps):
    # ---- Uniform Monte Carlo -----------------------------------------
    # Draw x uniformly on [a,b]; estimator is (b-a) * mean(f(x)).
    xu = np.random.uniform(a, b, n)          # uniform draws
    uniform_estimates[r] = (b - a) * np.mean(f(xu))

    # ---- Importance sampling -----------------------------------------
    # Inverse transform for p(x)=exp(-x): if u~U(0,1) then x=-ln(u) ~ Exp(1).
    u = np.random.uniform(0.0, 1.0, n)
    xi = -np.log(u)                          # draws from p
    within = (xi >= a) & (xi <= b)           # keep points inside [0,100]
    # Weighted estimator: mean of f(x)/p(x), zeroing contributions outside range.
    w = np.where(within, f(xi) / p(xi), 0.0)
    importance_estimates[r] = np.mean(w)

# ---- Summary statistics -----------------------------------------------
print(f"True value (sqrt(pi)/2)            : {true_val:.10f}")
print(f"Uniform    mean over replicates    : {np.mean(uniform_estimates):.10f}")
print(f"Uniform    std  over replicates    : {np.std(uniform_estimates, ddof=1):.10f}")
print(f"Uniform    variance over replicates: {np.var(uniform_estimates, ddof=1):.10e}")
print(f"Importance mean over replicates    : {np.mean(importance_estimates):.10f}")
print(f"Importance std  over replicates    : {np.std(importance_estimates, ddof=1):.10f}")
print(f"Importance variance over replicates: {np.var(importance_estimates, ddof=1):.10e}")
print(f"Variance reduction factor (unif/imp): "
      f"{np.var(uniform_estimates, ddof=1)/np.var(importance_estimates, ddof=1):.4f}")

# ---- Concentration check ----------------------------------------------
# Fraction of points landing where f is non-negligible (say x < 5), one rep.
np.random.seed(1)
xu = np.random.uniform(a, b, n)
u = np.random.uniform(0.0, 1.0, n)
xi = -np.log(u)
frac_unif_useful = np.mean(xu < 5.0)
frac_imp_useful = np.mean(xi < 5.0)
print(f"Fraction of uniform points with x<5    : {frac_unif_useful:.4f}")
print(f"Fraction of importance points with x<5 : {frac_imp_useful:.4f}")

# ---- Box plot ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([importance_estimates, uniform_estimates],
           labels=["Importance", "Uniform"])
ax.axhline(true_val, color="red", linestyle="--", label="true = sqrt(pi)/2")
ax.set_ylabel("Estimate of integral")
ax.set_title("Importance sampling vs uniform Monte Carlo (n=1e4, 100 reps)")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.4.1_s4.png")

# ---- One-sentence explanation -----------------------------------------
print("Explanation: The check confirms the result because nearly all importance "
      "points fall in x<5 where f is large (giving tightly clustered estimates), "
      "while almost all uniform points fall where f is negligible and are wasted, "
      "producing the wide, high-variance box plot.")
