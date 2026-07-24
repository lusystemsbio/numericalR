import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations, permutations

# --- City coordinates (n = 10) ---
x = np.array([0, -28.87, -14.66, -29.06, -36.04, -50.48, -50.59, -0.14, -21.50, -43.07])
y = np.array([0, 0, 43.39, 43.22, 21.61, -7.37, 21.59, 28.73, -7.32, -14.55])
n = len(x)

# --- Build the Euclidean distance matrix d[i][j] ---
coords = np.column_stack((x, y))
d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))

# --- Held-Karp dynamic programming ---
# We fix city 0 as start/end. States are (last_city, set_of_visited_cities_excluding_start).
# g[(x_last, S)] = length of shortest path that starts at 0, visits exactly the cities in S,
#                  and ends at city x_last (x_last is in S). S is a frozenset of city indices.
# Recurrence: g(x, S) = min over i in S\{x} of ( g(i, S\{x}) + d(i, x) ),
#             base case S = {x}: g(x, {x}) = d(0, x).

g = {}          # maps (last, frozenset) -> cost
parent = {}     # maps (last, frozenset) -> previous last city, for reconstruction

others = list(range(1, n))  # cities other than the start city 0

# Base cases: paths from 0 directly to a single city x
for xc in others:
    g[(xc, frozenset([xc]))] = d[0][xc]
    parent[(xc, frozenset([xc]))] = 0

# Grow subsets by increasing size
for size in range(2, n):  # subset sizes from 2 up to n-1 (all non-start cities)
    for subset in combinations(others, size):
        S = frozenset(subset)
        for xc in subset:                     # xc is the endpoint of the path
            Sprev = S - {xc}                  # subset before arriving at xc
            best_cost = float("inf")
            best_prev = None
            for i in subset:                  # i is the city visited just before xc
                if i == xc:
                    continue
                cost = g[(i, Sprev)] + d[i][xc]
                if cost < best_cost:
                    best_cost = cost
                    best_prev = i
            g[(xc, S)] = best_cost
            parent[(xc, S)] = best_prev

# Close the tour: return to city 0 from the last city, over the full set of other cities
full = frozenset(others)
best_cost = float("inf")
best_last = None
for xc in others:
    cost = g[(xc, full)] + d[xc][0]
    if cost < best_cost:
        best_cost = cost
        best_last = xc

# --- Reconstruct the optimal tour by walking parents backwards ---
tour = [0]
S = full
last = best_last
path = []
while last != 0:
    path.append(last)
    prev = parent[(last, S)]
    S = S - {last}
    last = prev
path.reverse()
tour = [0] + path + [0]

print("Optimal tour (0-indexed):", tour)
print("Optimal tour (1-indexed):", [c + 1 for c in tour])
print("Held-Karp shortest tour length:", best_cost)

# --- Independent check: brute-force over all permutations of the other cities ---
# For n = 10 this is 9! = 362880 tours, fully enumerable.
brute_best = float("inf")
brute_tour = None
for perm in permutations(others):
    length = d[0][perm[0]] + sum(d[perm[k]][perm[k + 1]] for k in range(len(perm) - 1)) + d[perm[-1]][0]
    if length < brute_best:
        brute_best = length
        brute_tour = (0,) + perm + (0,)

print("Brute-force shortest tour length:", brute_best)
print("Brute-force tour (1-indexed):", [c + 1 for c in brute_tour])
print("Held-Karp matches brute force:", np.isclose(best_cost, brute_best))
# One-sentence explanation: the brute-force check enumerates every possible tour and takes the
# minimum, so agreement with it proves Held-Karp found the true global optimum, not just a good tour.

# --- Draw the shortest tour ---
fig, ax = plt.subplots(figsize=(7, 7))
tx = x[tour]
ty = y[tour]
ax.plot(tx, ty, "-o", color="steelblue", markersize=8, zorder=1)
for idx in range(n):
    ax.annotate(str(idx + 1), (x[idx], y[idx]),
                textcoords="offset points", xytext=(6, 6), fontsize=11)
ax.plot(x[0], y[0], "s", color="red", markersize=12, label="start/end (city 1)", zorder=2)
ax.set_title(f"Held-Karp shortest TSP tour (length = {best_cost:.4f})")
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_aspect("equal", adjustable="datalim")
ax.legend(); ax.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/9B.3.1_s5.png")
