import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model: I = integral_0^100 f(x) dx, with f(x) = exp(-x^2) ---
def f(x):
    return np.exp(-x**2)

a, b = 0.0, 100.0          # wide integration range
n = int(1e4)               # samples per estimate
reps = 100                 # replicates per method
rng = np.random.default_rng(1)   # seed 1

# Reference value: integral_0^100 exp(-x^2) dx ~= sqrt(pi)/2 (tail beyond 100 is 0)
true_val = np.sqrt(np.pi) / 2.0
print(f"True integral (sqrt(pi)/2): {true_val:.10f}")

uniform_estimates = np.empty(reps)
importance_estimates = np.empty(reps)

for r in range(reps):
    # --- Uniform Monte Carlo ---
    # Draw x uniformly on [a, b]; estimator = (b-a) * mean(f(x))
    xu = rng.uniform(a, b, size=n)
    uniform_estimates[r] = (b - a) * np.mean(f(xu))

    # --- Importance sampling with p(x) = exp(-x) on [0, inf) ---
    # Inverse transform: if U~Uniform(0,1), then x = -ln(U) has density exp(-x)
    u = rng.uniform(0.0, 1.0, size=n)
    xi = -np.log(u)
    # Weighted estimator: mean of f(x)/p(x), where p(x) = exp(-x)
    # f(x)/p(x) = exp(-x^2) / exp(-x) = exp(-x^2 + x)
    weights = f(xi) / np.exp(-xi)
    # keep only draws inside [a, b] (exponential rarely exceeds 100; contribution beyond is ~0)
    weights = np.where((xi >= a) & (xi <= b), weights, 0.0)
    importance_estimates[r] = np.mean(weights)

# --- Report summary statistics ---
print(f"Uniform    mean estimate: {np.mean(uniform_estimates):.10f}")
print(f"Uniform    variance:      {np.var(uniform_estimates):.6e}")
print(f"Uniform    std dev:       {np.std(uniform_estimates):.6e}")
print(f"Importance mean estimate: {np.mean(importance_estimates):.10f}")
print(f"Importance variance:      {np.var(importance_estimates):.6e}")
print(f"Importance std dev:       {np.std(importance_estimates):.6e}")
print(f"Variance reduction factor (uniform/importance): "
      f"{np.var(uniform_estimates)/np.var(importance_estimates):.6e}")

# --- Check: where do the sampled points land relative to where f is large? ---
# f(x) is negligible (< 1e-6) beyond x ~ 3.7; count fraction of points that are "useful".
thresh = 3.7
last_u = rng.uniform(a, b, size=n)
last_i = -np.log(rng.uniform(0.0, 1.0, size=n))
print(f"Uniform    fraction of points with x < {thresh}: "
      f"{np.mean(last_u < thresh):.6f}")
print(f"Importance fraction of points with x < {thresh}: "
      f"{np.mean(last_i < thresh):.6f}")

# --- Box plots comparing the two methods ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([uniform_estimates, importance_estimates],
           labels=["Uniform", "Importance"])
ax.axhline(true_val, color="red", linestyle="--", label="true value")
ax.set_ylabel("Estimate of integral")
ax.set_title("Uniform vs Importance Sampling (100 replicates)")
ax.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8A.4.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the result because uniform sampling places "
      "almost all points where f(x) is essentially zero (wasting them), while the "
      "exponential importance density concentrates points near x=0 where f is large, "
      "so its estimates cluster tightly around the true value with far smaller variance.")
