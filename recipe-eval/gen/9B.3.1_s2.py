import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
from itertools import combinations

# ---- City coordinates (city 1 is index 0) ----
x = (0, -28.87, -14.66, -29.06, -36.04, -50.48, -50.59, -0.14, -21.50, -43.07)
y = (0, 0, 43.39, 43.22, 21.61, -7.37, 21.59, 28.73, -7.32, -14.55)
n = len(x)

# ---- Build the Euclidean distance matrix ----
d = [[math.hypot(x[i] - x[j], y[i] - y[j]) for j in range(n)] for i in range(n)]

# ---- Held-Karp DP ----
# We fix city 0 as start/end. States are (last_city, set_of_visited_besides_start).
# g[(x_city, S)] = min cost of a path that starts at 0, visits exactly the cities in S,
# and ends at x_city (with x_city in S). Sets are encoded as bitmasks.
g = {}
parent = {}  # to reconstruct the optimal tour

# Base case: paths 0 -> i using only city i (S = {i}).
for i in range(1, n):
    S = 1 << i
    g[(i, S)] = d[0][i]
    parent[(i, S)] = 0

# Grow subsets by increasing size, applying the recurrence
# g(x,S) = min over i in S\{x} of ( g(i, S\{x}) + d(i,x) ).
for size in range(2, n):
    for subset in combinations(range(1, n), size):
        S = 0
        for c in subset:
            S |= (1 << c)
        for last in subset:                 # 'last' is x in the recurrence
            prev_S = S & ~(1 << last)        # S without x
            best_cost = math.inf
            best_prev = None
            for i in subset:                 # i ranges over S\{x}
                if i == last:
                    continue
                cost = g[(i, prev_S)] + d[i][last]
                if cost < best_cost:
                    best_cost = cost
                    best_prev = i
            g[(last, S)] = best_cost
            parent[(last, S)] = best_prev

# ---- Close the tour: return to city 0 from the full set ----
full = (1 << n) - 2  # all cities 1..n-1 visited (bit 0 excluded)
best_total = math.inf
best_last = None
for last in range(1, n):
    cost = g[(last, full)] + d[last][0]
    if cost < best_total:
        best_total = cost
        best_last = last

# ---- Reconstruct the optimal tour by walking the parent pointers ----
tour = [0]
S = full
last = best_last
path = []
while last != 0:
    path.append(last)
    plast = parent[(last, S)]
    S &= ~(1 << last)
    last = plast
path.reverse()
tour = [0] + path + [0]  # 1-based cities: add 1 when printing

# ---- Brute-force check over all permutations (exact, so it must match) ----
from itertools import permutations
brute_best = math.inf
brute_tour = None
for perm in permutations(range(1, n)):
    seq = (0,) + perm + (0,)
    length = sum(d[seq[k]][seq[k+1]] for k in range(len(seq) - 1))
    if length < brute_best:
        brute_best = length
        brute_tour = seq

# ---- Print numerical results ----
print("Number of cities n:", n)
print("Held-Karp shortest tour length:", best_total)
print("Held-Karp tour (1-based cities):", [c + 1 for c in tour])
print("Brute-force shortest tour length:", brute_best)
print("Brute-force tour (1-based cities):", [c + 1 for c in brute_tour])
print("Lengths match (Held-Karp == brute force):", math.isclose(best_total, brute_best))

# ---- Plot the shortest tour ----
tx = [x[c] for c in tour]
ty = [y[c] for c in tour]
plt.figure(figsize=(7, 7))
plt.plot(tx, ty, "-o", color="tab:blue")
for c in range(n):
    plt.annotate(str(c + 1), (x[c], y[c]), textcoords="offset points", xytext=(6, 6))
plt.title("Held-Karp shortest TSP tour (length = %.4f)" % best_total)
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect("equal", adjustable="box")
plt.grid(True)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/9B.3.1_s2.png")

# The brute-force check enumerates every possible tour and picks the minimum, so its
# agreement with Held-Karp confirms the DP found the true global optimum, not just a
# local one.
print("Explanation: brute force examines every permutation of cities and takes the")
print("global minimum, so matching it proves the Held-Karp tour is exactly optimal.")
