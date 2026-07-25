import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Problem setup ---
np.random.seed(1)                       # reproducibility
DIM = 2                                  # number of variables
LOW, HIGH = -5.12, 5.12                  # search bounds per variable
RANGE = HIGH - LOW                       # width of each variable's domain

def rastrigin(x):
    # Rastrigin: global min f=0 at the origin, surrounded by many local minima
    return 10.0 * DIM + np.sum(x**2 - 10.0 * np.cos(2.0 * np.pi * x), axis=-1)

# --- GA parameters ---
POP = 50                                 # population size
CROSS_RATE = 0.7                         # probability two parents exchange a variable
GENS = 100                               # number of generations
MUT_STEP = 0.10 * RANGE                  # mutation nudge = 10% of the range
ELITE = 1                                # number of fittest kept unchanged each generation

# --- Initialize population uniformly in the box ---
pop = np.random.uniform(LOW, HIGH, size=(POP, DIM))

def evaluate(P):
    return np.array([rastrigin(ind) for ind in P])

best_history = []                        # best score per generation

for gen in range(GENS):
    fitness = evaluate(pop)              # lower is better (minimization)
    order = np.argsort(fitness)          # sort individuals best -> worst
    pop = pop[order]
    fitness = fitness[order]
    best_history.append(fitness[0])      # record current best

    # Elitist selection: carry the fittest members over unchanged
    new_pop = [pop[i].copy() for i in range(ELITE)]

    # Fill the rest of the next generation
    while len(new_pop) < POP:
        # Tournament selection: pick 2 random contenders, keep the better one (twice)
        def select():
            i, j = np.random.randint(0, POP, size=2)
            return pop[i].copy() if fitness[i] < fitness[j] else pop[j].copy()
        p1, p2 = select(), select()

        # Crossover: with prob CROSS_RATE swap a whole variable between the two parents
        c1, c2 = p1.copy(), p2.copy()
        if np.random.rand() < CROSS_RATE:
            k = np.random.randint(0, DIM)     # variable index to swap
            c1[k], c2[k] = p2[k], p1[k]

        # Mutation: nudge one variable of each child by a Gaussian step, clip to bounds
        for child in (c1, c2):
            m = np.random.randint(0, DIM)     # which variable to nudge
            child[m] += np.random.randn() * MUT_STEP
            child[m] = np.clip(child[m], LOW, HIGH)

        new_pop.append(c1)
        if len(new_pop) < POP:
            new_pop.append(c2)

    pop = np.array(new_pop)

# --- Final evaluation ---
final_fitness = evaluate(pop)
best_idx = np.argmin(final_fitness)
best_solution = pop[best_idx]
best_score = final_fitness[best_idx]
best_history.append(best_score)

# --- Report results ---
print(f"Best score (final):            {best_score:.6e}")
print(f"Best solution x1:              {best_solution[0]:.6e}")
print(f"Best solution x2:              {best_solution[1]:.6e}")
print(f"Best score at generation 0:    {best_history[0]:.6e}")
print(f"Best score at generation 50:   {best_history[50]:.6e}")
print(f"Best score at final generation:{best_history[-1]:.6e}")
print(f"Distance of best sol from origin: {np.linalg.norm(best_solution):.6e}")

# --- Convergence check ---
# The best score dropping toward zero confirms the result because f=0 is attained
# only at the global minimum (the origin), so approaching zero means the GA escaped
# the local minima and located that global optimum.
converged = best_score < 1e-1
print(f"Converged toward global minimum (best < 0.1): {converged}")
print("Why: f=0 occurs only at the origin, so a best score approaching zero "
      "confirms the GA found the global minimum rather than a local one.")

# --- Plot best score versus generation ---
plt.figure(figsize=(8, 5))
plt.plot(range(len(best_history)), best_history, marker=".", color="crimson")
plt.xlabel("Generation")
plt.ylabel("Best score f(x)")
plt.title("GA on Rastrigin function: best score vs generation")
plt.yscale("log")
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.1.1_s3.png")
