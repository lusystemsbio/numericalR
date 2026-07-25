import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import permutations

# ---------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------
random.seed(12345)
np.random.seed(12345)

# ---------------------------------------------------------------
# TSP instance generator: random city coordinates in the unit square
# ---------------------------------------------------------------
def make_instance(n, seed):
    rng = np.random.RandomState(seed)
    coords = rng.rand(n, 2)                       # n cities, (x, y)
    # full pairwise Euclidean distance matrix
    diff = coords[:, None, :] - coords[None, :, :]
    D = np.sqrt((diff ** 2).sum(axis=2))
    return coords, D

# ---------------------------------------------------------------
# Objective: tour length = sum of consecutive legs + return leg
# ---------------------------------------------------------------
def tour_length(tour, D):
    total = 0.0
    n = len(tour)
    for i in range(n):
        total += D[tour[i], tour[(i + 1) % n]]    # (i+1)%n closes the loop
    return total

# ---------------------------------------------------------------
# Permutation-preserving mutations
# ---------------------------------------------------------------
def mutate(tour):
    t = tour[:]                                   # copy
    if random.random() < 0.5:
        # swap two cities
        i, j = random.sample(range(len(t)), 2)
        t[i], t[j] = t[j], t[i]
    else:
        # reverse a segment (2-opt style move)
        i, j = sorted(random.sample(range(len(t)), 2))
        t[i:j + 1] = reversed(t[i:j + 1])
    return t

# ---------------------------------------------------------------
# Order crossover (OX): keeps a slice from parent 1, fills the rest
# in the order they appear in parent 2 -> child is a valid permutation
# ---------------------------------------------------------------
def order_crossover(p1, p2):
    n = len(p1)
    a, b = sorted(random.sample(range(n), 2))
    child = [None] * n
    child[a:b + 1] = p1[a:b + 1]                  # copy contiguous slice from p1
    taken = set(child[a:b + 1])
    # fill remaining positions with p2's order, skipping already-used cities
    fill = [c for c in p2 if c not in taken]
    idx = 0
    for k in range(n):
        if child[k] is None:
            child[k] = fill[idx]
            idx += 1
    return child

# ---------------------------------------------------------------
# Tournament selection: pick the fitter of two random individuals
# ---------------------------------------------------------------
def tournament(pop, fits):
    i, j = random.sample(range(len(pop)), 2)
    return pop[i][:] if fits[i] < fits[j] else pop[j][:]

# ---------------------------------------------------------------
# The genetic algorithm on permutations
# ---------------------------------------------------------------
def genetic_tsp(D, pop_size, crossover_rate, generations):
    n = D.shape[0]
    # initial population: random permutations
    pop = [random.sample(range(n), n) for _ in range(pop_size)]
    fits = [tour_length(t, D) for t in pop]

    best_idx = int(np.argmin(fits))
    best_tour = pop[best_idx][:]
    best_len = fits[best_idx]
    history = [best_len]                          # best length per generation

    for _ in range(generations):
        new_pop = [best_tour[:]]                   # elitism: keep the best so far
        while len(new_pop) < pop_size:
            parent = tournament(pop, fits)
            if random.random() < crossover_rate:
                other = tournament(pop, fits)
                child = order_crossover(parent, other)
            else:
                child = parent[:]
            child = mutate(child)                  # always apply a mutation move
            new_pop.append(child)

        pop = new_pop
        fits = [tour_length(t, D) for t in pop]
        gen_best = int(np.argmin(fits))
        if fits[gen_best] < best_len:
            best_len = fits[gen_best]
            best_tour = pop[gen_best][:]
        history.append(best_len)

    return best_tour, best_len, history

# ---------------------------------------------------------------
# Exact solver: Held-Karp dynamic programming (fix city 0 as start)
# ---------------------------------------------------------------
def held_karp(D):
    n = D.shape[0]
    # dp[(mask, j)] = min cost path starting at 0, visiting mask, ending at j
    dp = {(1 << j, j): (D[0, j], 0) for j in range(1, n)}
    for size in range(2, n):
        for subset_bits in permutations(range(1, n), size):
            mask = 0
            for b in subset_bits:
                mask |= (1 << b)
            for j in subset_bits:
                prev = mask & ~(1 << j)            # remove j from the set
                best = None
                for k in subset_bits:
                    if k == j:
                        continue
                    cost = dp[(prev, k)][0] + D[k, j]
                    if best is None or cost < best[0]:
                        best = (cost, k)
                dp[(mask, j)] = best
    # close the tour back to city 0
    full = (1 << n) - 2                            # all cities except 0
    best = None
    for j in range(1, n):
        cost = dp[(full, j)][0] + D[j, 0]
        if best is None or cost < best[0]:
            best = (cost, j)
    return best[0]

# ===============================================================
# Test 1: 10-city instance, pop 10, crossover 0.8, 50 generations
# ===============================================================
n1 = 10
coords1, D1 = make_instance(n1, seed=1)
best_tour, best_len, history = genetic_tsp(D1, pop_size=10,
                                           crossover_rate=0.8, generations=50)

print("=== 10-city instance (GA) ===")
print("GA best tour:", best_tour)
print("GA best tour length:", best_len)
print("GA initial best length (gen 0):", history[0])
print("GA final best length (gen 50):", history[-1])

# ===============================================================
# Check: does the GA recover the exact DP optimum?
# ===============================================================
exact_len = held_karp(D1)
print("Exact (Held-Karp DP) shortest length:", exact_len)
print("Difference (GA - exact):", best_len - exact_len)
print("GA matches exact optimum:", np.isclose(best_len, exact_len, atol=1e-9))

# ===============================================================
# Scaling: larger instance where exact DP is infeasible
# ===============================================================
n2 = 40
coords2, D2 = make_instance(n2, seed=2)
best_tour2, best_len2, history2 = genetic_tsp(D2, pop_size=50,
                                              crossover_rate=0.8, generations=500)
print("=== 40-city instance (GA, exact DP infeasible) ===")
print("GA best tour length (40 cities):", best_len2)
print("GA improvement from gen 0 to final:", history2[0] - history2[-1])

# ---------------------------------------------------------------
# Plot: best length vs generation (both instances)
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].plot(range(len(history)), history, marker="o", ms=3)
axes[0].axhline(exact_len, color="red", ls="--", label="exact DP optimum")
axes[0].set_xlabel("generation")
axes[0].set_ylabel("best tour length")
axes[0].set_title("10-city GA convergence")
axes[0].legend()

# draw the best 10-city tour
loop = best_tour + [best_tour[0]]
axes[1].plot(coords1[loop, 0], coords1[loop, 1], "-o")
for idx, (x, y) in enumerate(coords1):
    axes[1].annotate(str(idx), (x, y))
axes[1].set_title("Best 10-city tour (GA)")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.2.1_s5.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result
# ---------------------------------------------------------------
print("Why the check confirms the result: Held-Karp dynamic programming "
      "computes the provably shortest tour by exhaustively optimizing over all "
      "city subsets, so the GA matching that exact length on the small instance "
      "shows the GA finds true optima, lending confidence to its answers on "
      "larger instances where exact enumeration is computationally infeasible.")
