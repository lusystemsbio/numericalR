import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Reproducibility ---
rng = np.random.default_rng(1)

# --- Problem definition: Rastrigin function in 2D ---
# f(x) = 20 + sum(x_i^2 - 10*cos(2*pi*x_i)); global min f=0 at the origin.
def rastrigin(x):
    x = np.asarray(x)
    return 10.0 * len(x) + np.sum(x**2 - 10.0 * np.cos(2.0 * np.pi * x))

# --- GA parameters ---
n_vars      = 2                       # number of continuous variables
lower, upper = -5.12, 5.12            # search box [-5.12, 5.12]^2
pop_size    = 50                      # population size
p_cross     = 0.7                     # crossover rate
n_gen       = 100                     # number of generations
mut_step    = 0.10 * (upper - lower)  # mutation step = 10% of the range
n_elite     = 2                       # keep the fittest members (elitism)

# --- Initialize population uniformly in the box ---
pop = rng.uniform(lower, upper, size=(pop_size, n_vars))

# Evaluate initial fitness (lower is better since we minimize)
fitness = np.array([rastrigin(ind) for ind in pop])

best_history = []  # best score per generation

# --- Main GA loop ---
for gen in range(n_gen):
    # Sort population by fitness (ascending: best first)
    order = np.argsort(fitness)
    pop, fitness = pop[order], fitness[order]

    # Record the best score of this generation
    best_history.append(fitness[0])

    # Start next generation with elites (fittest members kept unchanged)
    new_pop = [pop[i].copy() for i in range(n_elite)]

    # Fill the rest of the population with offspring
    while len(new_pop) < pop_size:
        # Tournament selection: pick 2 parents, each the better of a random pair
        def tournament():
            a, b = rng.integers(0, pop_size, size=2)
            return pop[a] if fitness[a] < fitness[b] else pop[b]
        parent1 = tournament().copy()
        parent2 = tournament().copy()

        # Crossover: with probability p_cross, swap whole variables between parents
        child = parent1.copy()
        if rng.random() < p_cross:
            for j in range(n_vars):
                if rng.random() < 0.5:      # swap this whole variable
                    child[j] = parent2[j]

        # Mutation: nudge one randomly chosen variable by a Gaussian step
        j = rng.integers(0, n_vars)
        child[j] += rng.normal(0.0, mut_step)

        # Keep child inside the box
        child = np.clip(child, lower, upper)

        new_pop.append(child)

    # Replace population and re-evaluate fitness
    pop = np.array(new_pop)
    fitness = np.array([rastrigin(ind) for ind in pop])

# --- Final best solution ---
order = np.argsort(fitness)
pop, fitness = pop[order], fitness[order]
best_history.append(fitness[0])
best_solution, best_score = pop[0], fitness[0]

# --- Report results ---
print(f"Best solution x1: {best_solution[0]:.6f}")
print(f"Best solution x2: {best_solution[1]:.6f}")
print(f"Best score (final): {best_score:.6f}")
print(f"Best score (initial gen 0): {best_history[0]:.6f}")
print(f"Distance of best solution from origin: {np.linalg.norm(best_solution):.6f}")

# --- Convergence check ---
# The check: the best score must fall from its initial large value toward ~0.
improvement = best_history[0] - best_history[-1]
converged = best_score < 1.0  # within a small tolerance of the global minimum f=0
print(f"Improvement from first to last generation: {improvement:.6f}")
print(f"Converged to global minimum (best score < 1.0): {converged}")
# One-sentence explanation:
print("Explanation: because the Rastrigin minimum is exactly 0 at the origin, "
      "a best score driven toward 0 confirms the GA escaped the many local minima "
      "and located the true global optimum.")

# --- Plot best score vs generation ---
plt.figure(figsize=(8, 5))
plt.plot(range(len(best_history)), best_history, marker=".", lw=1)
plt.xlabel("Generation")
plt.ylabel("Best score f(x)")
plt.title("GA on 2D Rastrigin: best score vs generation")
plt.yscale("log")
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.1.1_s5.png")
