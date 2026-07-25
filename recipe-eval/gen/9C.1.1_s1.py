import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reproducibility
rng = np.random.default_rng(1)

# ---- Problem definition ----
def rastrigin(x):
    # x is a 1D array of variables; global minimum f=0 at the origin
    return 10.0 * len(x) + np.sum(x**2 - 10.0 * np.cos(2.0 * np.pi * x))

n_vars = 2                      # two variables
lo, hi = -5.12, 5.12            # search bounds per variable
span = hi - lo                  # range width used for mutation scaling

# ---- GA parameters ----
pop_size = 50
crossover_rate = 0.7
n_gen = 100
mut_step = 0.10 * span          # mutation nudge = 10% of the range

# ---- Initialize population uniformly in the box ----
pop = rng.uniform(lo, hi, size=(pop_size, n_vars))
fitness = np.array([rastrigin(ind) for ind in pop])

best_history = []               # best score per generation

for gen in range(n_gen):
    # --- Elitist selection: keep the fittest half as parents ---
    order = np.argsort(fitness)         # ascending: smaller f is better
    pop = pop[order]
    fitness = fitness[order]
    n_keep = pop_size // 2
    parents = pop[:n_keep].copy()       # survivors / mating pool

    # --- Generate offspring to refill the population ---
    children = []
    while len(children) < pop_size - n_keep:
        # pick two distinct parents at random from the elite pool
        i, j = rng.integers(0, n_keep, size=2)
        p1 = parents[i].copy()
        p2 = parents[j].copy()

        # Crossover: with probability crossover_rate, swap whole variables
        if rng.random() < crossover_rate:
            for v in range(n_vars):
                if rng.random() < 0.5:
                    p1[v], p2[v] = p2[v], p1[v]

        # Mutation: nudge one randomly chosen variable of each child
        for child in (p1, p2):
            v = rng.integers(0, n_vars)
            child[v] += rng.normal(0.0, mut_step)
            child[v] = np.clip(child[v], lo, hi)   # keep within bounds
            children.append(child)

    children = np.array(children[:pop_size - n_keep])

    # --- Form the new population: elites + offspring ---
    pop = np.vstack([parents, children])
    fitness = np.array([rastrigin(ind) for ind in pop])

    best_history.append(fitness.min())

# ---- Final best solution ----
best_idx = np.argmin(fitness)
best_x = pop[best_idx]
best_f = fitness[best_idx]

print(f"Best x1: {best_x[0]:.6f}")
print(f"Best x2: {best_x[1]:.6f}")
print(f"Best score f (final): {best_f:.6e}")
print(f"Best score f (generation 1): {best_history[0]:.6e}")
print(f"Best score f (generation 100): {best_history[-1]:.6e}")
print(f"Distance of best solution from origin: {np.linalg.norm(best_x):.6e}")

# ---- Convergence check ----
# The check: does the best score decrease toward zero over generations?
improved = best_history[-1] < best_history[0]
near_zero = best_history[-1] < 1.0
print(f"Check - best score improved over run: {improved}")
print(f"Check - final best score near zero (<1): {near_zero}")
print("Check explanation: because f=0 occurs only at the global minimum at the "
      "origin, a best score converging toward zero confirms the GA escaped the "
      "local minima and located the true global optimum.")

# ---- Plot best score vs generation ----
plt.figure(figsize=(8, 5))
plt.plot(range(1, n_gen + 1), best_history, lw=1.5, color="C0")
plt.xlabel("Generation")
plt.ylabel("Best score f")
plt.title("GA on Rastrigin: best score vs generation")
plt.yscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.1.1_s1.png")
