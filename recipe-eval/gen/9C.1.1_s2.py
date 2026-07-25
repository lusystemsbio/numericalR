import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Rastrigin function: global minimum f=0 at the origin, many local minima
# ---------------------------------------------------------------
def rastrigin(x):
    return 20 + np.sum(x**2 - 10 * np.cos(2 * np.pi * x))

# ---------------------------------------------------------------
# GA parameters
# ---------------------------------------------------------------
np.random.seed(1)               # reproducibility
n_vars = 2                      # two variables
lo, hi = -5.12, 5.12            # search bounds per variable
pop_size = 50                   # population size
cx_rate = 0.7                   # crossover rate
n_gen = 100                     # number of generations
mut_step = 0.10 * (hi - lo)     # mutation step = 10% of the range

# ---------------------------------------------------------------
# Initialize population uniformly in the box, evaluate fitness
# ---------------------------------------------------------------
pop = np.random.uniform(lo, hi, size=(pop_size, n_vars))
fit = np.array([rastrigin(ind) for ind in pop])

best_history = []               # best score per generation

# ---------------------------------------------------------------
# Evolution loop
# ---------------------------------------------------------------
for gen in range(n_gen):
    # --- Elitist selection: keep the fittest members, rank by fitness (lower is better)
    order = np.argsort(fit)
    pop = pop[order]
    fit = fit[order]

    # --- Build the next generation, keeping the elite (best) individual unchanged
    new_pop = [pop[0].copy()]
    while len(new_pop) < pop_size:
        # tournament-style parent picks biased toward the fitter (front of sorted pop)
        i = np.random.randint(0, pop_size // 2)
        j = np.random.randint(0, pop_size // 2)
        parent1 = pop[i].copy()
        parent2 = pop[j].copy()

        # --- Crossover: swap whole variables between the two parents
        child = parent1.copy()
        if np.random.rand() < cx_rate:
            for k in range(n_vars):
                if np.random.rand() < 0.5:
                    child[k] = parent2[k]   # take this whole variable from parent2

        # --- Mutation: nudge one randomly chosen variable by a Gaussian step
        k = np.random.randint(0, n_vars)
        child[k] += np.random.normal(0.0, mut_step)
        child = np.clip(child, lo, hi)      # keep inside the bounds

        new_pop.append(child)

    # --- Replace population and re-evaluate fitness
    pop = np.array(new_pop)
    fit = np.array([rastrigin(ind) for ind in pop])

    best_history.append(fit.min())

# ---------------------------------------------------------------
# Report results
# ---------------------------------------------------------------
best_idx = np.argmin(fit)
best_x = pop[best_idx]
best_score = fit[best_idx]

print(f"Best score (final generation): {best_score:.6f}")
print(f"Best solution x1: {best_x[0]:.6f}")
print(f"Best solution x2: {best_x[1]:.6f}")
print(f"Best score at generation 1:   {best_history[0]:.6f}")
print(f"Best score at generation 50:  {best_history[49]:.6f}")
print(f"Best score at generation 100: {best_history[-1]:.6f}")

# --- Separate check: has the best score dropped toward zero (the global minimum)?
tol = 1e-2
converged = best_score < tol
print(f"Convergence tolerance: {tol:.6f}")
print(f"Converged toward global minimum (best_score < tol): {converged}")
# One-sentence explanation:
# Because f=0 is the true global minimum of the Rastrigin function and occurs
# only at the origin, a final best score driven below tolerance confirms the GA
# escaped the many local minima and located that global optimum.
print("Explanation: since f=0 is the unique global minimum at the origin, a best "
      "score falling near zero confirms the GA found the global optimum rather than a local one.")

# ---------------------------------------------------------------
# Plot best score versus generation
# ---------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(range(1, n_gen + 1), best_history, color="darkred", lw=1.5)
plt.axhline(0.0, color="gray", ls="--", lw=1, label="global minimum f=0")
plt.xlabel("Generation")
plt.ylabel("Best score")
plt.title("Genetic Algorithm on the Rastrigin Function")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.1.1_s2.png")
