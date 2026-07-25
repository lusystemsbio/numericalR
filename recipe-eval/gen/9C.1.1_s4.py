import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Problem setup ----
rng = np.random.default_rng(1)          # seed 1 for reproducibility
lo, hi = -5.12, 5.12                    # bounds for each variable
ndim = 2                                # two variables
rng_span = hi - lo                      # range width per variable

def rastrigin(x):
    # f = 20 + sum(xi^2 - 10 cos(2 pi xi)); global min f=0 at origin
    return 10 * ndim + np.sum(x**2 - 10 * np.cos(2 * np.pi * x), axis=-1)

# ---- GA parameters ----
pop_size = 50
crossover_rate = 0.7
generations = 100
mut_step = 0.10 * rng_span               # mutation nudges by 10% of the range
n_elite = 2                              # keep the fittest members (elitism)

# ---- Initialize population uniformly in the box ----
pop = rng.uniform(lo, hi, size=(pop_size, ndim))
best_history = np.empty(generations)     # best score per generation

def tournament(fitness):
    # pick the better of two random individuals (selection pressure)
    a, b = rng.integers(0, pop_size, size=2)
    return a if fitness[a] < fitness[b] else b

for gen in range(generations):
    fitness = rastrigin(pop)             # lower is better (minimization)
    order = np.argsort(fitness)          # rank individuals best-first
    best_history[gen] = fitness[order[0]]

    # Elitist selection: carry the fittest members over unchanged
    new_pop = pop[order[:n_elite]].copy()

    # Fill the rest of the next generation
    while len(new_pop) < pop_size:
        p1 = pop[tournament(fitness)].copy()
        p2 = pop[tournament(fitness)].copy()

        # Crossover: swap whole variables between parents
        if rng.random() < crossover_rate:
            for i in range(ndim):
                if rng.random() < 0.5:
                    p1[i], p2[i] = p2[i], p1[i]

        # Mutation: nudge one variable by a Gaussian step, then clip to bounds
        for child in (p1, p2):
            j = rng.integers(0, ndim)
            child[j] += rng.normal(0, mut_step)
            np.clip(child, lo, hi, out=child)

        new_pop = np.vstack([new_pop, p1, p2])

    pop = new_pop[:pop_size]             # keep population size fixed

# ---- Final best solution ----
final_fitness = rastrigin(pop)
best_idx = int(np.argmin(final_fitness))
best_x = pop[best_idx]
best_score = final_fitness[best_idx]

# ---- Report results ----
print(f"Best score (final generation): {best_score:.6e}")
print(f"Best x1: {best_x[0]:.6f}")
print(f"Best x2: {best_x[1]:.6f}")
print(f"Distance of best point from origin: {np.linalg.norm(best_x):.6e}")
print(f"Best score at generation 1:   {best_history[0]:.6e}")
print(f"Best score at generation 50:  {best_history[49]:.6e}")
print(f"Best score at generation 100: {best_history[-1]:.6e}")

# Convergence check: the best score should trend toward zero
converged = best_score < 1e-2
print(f"Converged toward zero (best score < 1e-2): {converged}")
# One-sentence explanation:
print("This check confirms the result because the Rastrigin function equals zero only at "
      "its unique global minimum (the origin), so a best score driven toward zero means the "
      "GA escaped the local minima and located that global optimum.")

# ---- Plot best score versus generation ----
plt.figure(figsize=(8, 5))
plt.plot(np.arange(1, generations + 1), best_history, color="tab:blue")
plt.yscale("log")
plt.xlabel("Generation")
plt.ylabel("Best score (log scale)")
plt.title("GA on Rastrigin: best score vs generation")
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.1.1_s4.png")
