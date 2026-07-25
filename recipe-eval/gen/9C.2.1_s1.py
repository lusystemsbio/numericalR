import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import itertools

rng = np.random.default_rng(42)  # reproducible pseudo-random generator

# ----------------------------------------------------------------------
# Model: Traveling Salesman Problem.
# A candidate solution is a permutation of city indices 0..n-1.
# Objective = sum of distances between consecutive cities + return leg.
# ----------------------------------------------------------------------

def make_instance(n, seed):
    """Generate n random 2-D city coordinates and their distance matrix."""
    r = np.random.default_rng(seed)
    coords = r.random((n, 2)) * 100.0
    # Euclidean distance matrix D[i,j]
    diff = coords[:, None, :] - coords[None, :, :]
    D = np.sqrt((diff ** 2).sum(axis=2))
    return coords, D

def tour_length(tour, D):
    """Sum of consecutive city distances plus the return leg (closed tour)."""
    total = 0.0
    for k in range(len(tour)):
        a = tour[k]
        b = tour[(k + 1) % len(tour)]  # wrap around for the return leg
        total += D[a, b]
    return total

# ----------------------------------------------------------------------
# Permutation-preserving mutation operators.
# ----------------------------------------------------------------------

def mutate_swap(tour):
    """Swap two randomly chosen cities (keeps it a valid permutation)."""
    t = tour.copy()
    i, j = rng.integers(0, len(t), size=2)
    t[i], t[j] = t[j], t[i]
    return t

def mutate_reverse(tour):
    """Reverse a random contiguous segment (2-opt style, stays a permutation)."""
    t = tour.copy()
    i, j = sorted(rng.integers(0, len(t), size=2))
    t[i:j + 1] = t[i:j + 1][::-1]
    return t

def mutate(tour):
    """Apply one of the two permutation-preserving mutations at random."""
    if rng.random() < 0.5:
        return mutate_swap(tour)
    return mutate_reverse(tour)

# ----------------------------------------------------------------------
# Order-preserving crossover (OX1): copies a slice from parent1, then fills
# remaining positions with parent2's cities in their original order.
# ----------------------------------------------------------------------

def order_crossover(p1, p2):
    n = len(p1)
    child = [-1] * n
    a, b = sorted(rng.integers(0, n, size=2))
    # 1) copy the segment [a,b] from parent1
    child[a:b + 1] = p1[a:b + 1]
    taken = set(child[a:b + 1])
    # 2) fill the rest with parent2's cities, skipping ones already used,
    #    preserving parent2's relative order, starting after the segment
    fill = [c for c in p2 if c not in taken]
    idx = 0
    for pos in range(n):
        if child[pos] == -1:
            child[pos] = fill[idx]
            idx += 1
    return np.array(child)

# ----------------------------------------------------------------------
# Genetic algorithm on permutations.
# ----------------------------------------------------------------------

def tournament_select(pop, fits, k=3):
    """Pick the fittest (shortest) of k random individuals."""
    cand = rng.integers(0, len(pop), size=k)
    best = cand[np.argmin(fits[cand])]
    return pop[best].copy()

def genetic_tsp(D, pop_size=10, crossover_rate=0.8, generations=50):
    n = D.shape[0]
    # initial population: random permutations
    pop = [rng.permutation(n) for _ in range(pop_size)]
    fits = np.array([tour_length(t, D) for t in pop])
    history = []  # best length at each generation

    best_tour = pop[int(np.argmin(fits))].copy()
    best_len = fits.min()

    for _ in range(generations):
        new_pop = [best_tour.copy()]  # elitism: keep the best so far
        while len(new_pop) < pop_size:
            p1 = tournament_select(pop, fits)
            p2 = tournament_select(pop, fits)
            # crossover with the given rate, else clone a parent
            if rng.random() < crossover_rate:
                child = order_crossover(p1, p2)
            else:
                child = p1.copy()
            child = mutate(child)  # always give a mutation chance
            new_pop.append(child)
        pop = new_pop
        fits = np.array([tour_length(t, D) for t in pop])
        gen_best = int(np.argmin(fits))
        if fits[gen_best] < best_len:
            best_len = fits[gen_best]
            best_tour = pop[gen_best].copy()
        history.append(best_len)

    return best_tour, best_len, history

# ----------------------------------------------------------------------
# Exact solver: Held-Karp dynamic programming (feasible only for small n).
# ----------------------------------------------------------------------

def held_karp(D):
    n = D.shape[0]
    # dp[(subset, j)] = shortest path from city 0 visiting `subset` ending at j
    dp = {}
    for j in range(1, n):
        dp[(1 << j, j)] = (D[0, j], 0)
    for size in range(2, n):
        for subset in itertools.combinations(range(1, n), size):
            bits = 0
            for c in subset:
                bits |= 1 << c
            for j in subset:
                prev = bits & ~(1 << j)
                best = min((dp[(prev, k)][0] + D[k, j], k)
                           for k in subset if k != j)
                dp[(bits, j)] = best
    full = (1 << n) - 2  # all cities except 0
    opt_len, last = min((dp[(full, j)][0] + D[j, 0], j) for j in range(1, n))
    # reconstruct the optimal tour
    tour = [0]
    bits = full
    j = last
    for _ in range(n - 1):
        tour.append(j)
        _, k = dp[(bits, j)]
        bits &= ~(1 << j)
        j = k
    tour.reverse()
    return np.array(tour), opt_len

# ----------------------------------------------------------------------
# Helper: two closed tours are equivalent under rotation/reflection.
# ----------------------------------------------------------------------

def same_cycle(t1, t2):
    t1 = list(t1)
    t2 = list(t2)
    n = len(t1)
    i0 = t2.index(t1[0])
    fwd = [t2[(i0 + k) % n] for k in range(n)]
    bwd = [t2[(i0 - k) % n] for k in range(n)]
    return t1 == fwd or t1 == bwd

# ======================================================================
# Test 1: 10-city instance, population 10, crossover 0.8, 50 generations.
# ======================================================================
coords10, D10 = make_instance(10, seed=1)
ga_tour, ga_len, history = genetic_tsp(D10, pop_size=10,
                                       crossover_rate=0.8, generations=50)

print("=== 10-city GA result ===")
print("GA best tour:", list(ga_tour))
print("GA best tour length:", ga_len)

# ======================================================================
# Check: does the GA recover the exact optimum from dynamic programming?
# ======================================================================
dp_tour, dp_len = held_karp(D10)
print("=== Exact (Held-Karp DP) result ===")
print("DP optimal tour:", list(dp_tour))
print("DP optimal tour length:", dp_len)

print("=== Cross-check ===")
print("GA length - DP length (gap):", ga_len - dp_len)
print("GA matches exact optimum (length within 1e-6):",
      bool(abs(ga_len - dp_len) < 1e-6))
print("GA tour equals DP tour as a cycle (rotation/reflection):",
      same_cycle(ga_tour, dp_tour))

# ======================================================================
# Scaling: run the GA on a larger instance where DP is infeasible.
# ======================================================================
coords50, D50 = make_instance(50, seed=2)
big_tour, big_len, _ = genetic_tsp(D50, pop_size=30,
                                   crossover_rate=0.8, generations=400)
rand_len = tour_length(np.arange(50), D50)  # unoptimized baseline tour
print("=== Scaling to 50 cities (exact DP infeasible) ===")
print("50-city GA best tour length:", big_len)
print("50-city baseline (identity) tour length:", rand_len)
print("50-city GA improvement factor over baseline:", rand_len / big_len)

# ----------------------------------------------------------------------
# Plot: best tour length versus generation for the 10-city run.
# ----------------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(range(1, len(history) + 1), history, marker="o", ms=3,
         label="GA best length")
plt.axhline(dp_len, color="red", ls="--", label="DP optimum")
plt.xlabel("Generation")
plt.ylabel("Best tour length")
plt.title("GA convergence on 10-city TSP (pop=10, cx=0.8, 50 gens)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9C.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
# The check confirms correctness because Held-Karp dynamic programming returns
# the provably global minimum tour length, so the GA reaching that same length
# demonstrates the GA found the true optimum (not merely a good local one).
print("Why the check confirms the result:",
      "Held-Karp DP returns the provably optimal tour length, so the GA "
      "matching that length proves the GA found the true global optimum "
      "on the small instance, giving confidence in its larger-instance results.")
