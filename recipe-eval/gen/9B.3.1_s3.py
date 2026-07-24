import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations, permutations
import math

# --- City coordinates (city 1 is index 0) ---
x = (0, -28.87, -14.66, -29.06, -36.04, -50.48, -50.59, -0.14, -21.50, -43.07)
y = (0, 0, 43.39, 43.22, 21.61, -7.37, 21.59, 28.73, -7.32, -14.55)
n = len(x)

# --- Build the Euclidean distance matrix d[i][j] ---
d = [[math.hypot(x[i] - x[j], y[i] - y[j]) for j in range(n)] for i in range(n)]

# --- Held-Karp DP ---
# We fix city 0 as start/end. States are (subset S of visited cities, endpoint x).
# g[(S, x)] = length of shortest path that starts at city 0, visits exactly the
# cities in S (S includes x but not 0), and ends at city x.
# Recurrence: g(x, S) = min over i in S\{x} of ( g(i, S\{x}) + d[i][x] ),
# with base case g(x, {x}) = d[0][x].
g = {}
parent = {}  # remember the predecessor i that achieved each minimum, to rebuild the tour

# Base cases: go directly from city 0 to each other city x
for k in range(1, n):
    g[(frozenset([k]), k)] = d[0][k]
    parent[(frozenset([k]), k)] = 0

# Fill DP over subsets of increasing size (all subsets of {1..n-1})
others = range(1, n)
for size in range(2, n):
    for subset in combinations(others, size):
        S = frozenset(subset)
        for x_end in subset:
            S_prev = S - {x_end}          # subset before arriving at x_end
            best, best_i = math.inf, None
            for i in S_prev:              # try each possible predecessor i in S
                cost = g[(S_prev, i)] + d[i][x_end]
                if cost < best:
                    best, best_i = cost, i
            g[(S, x_end)] = best
            parent[(S, x_end)] = best_i

# Close the tour: return to city 0 from the last city, over the full set {1..n-1}
full = frozenset(others)
best_len, last = math.inf, None
for x_end in others:
    cost = g[(full, x_end)] + d[x_end][0]
    if cost < best_len:
        best_len, last = cost, x_end

# --- Reconstruct the optimal tour by walking parent pointers backward ---
tour = [0]
S = full
cur = last
rev = []
while cur is not None and cur != 0:
    rev.append(cur)
    prev = parent[(S, cur)]
    S = S - {cur}
    cur = prev
tour += rev[::-1][::-1]        # rev is from last back toward start; reverse to forward order
tour = [0] + rev[::-1] + [0]   # start at 0, forward path, back to 0

print("Distance matrix (rows/cols = cities 1..10):")
for row in d:
    print("  " + "  ".join(f"{v:6.2f}" for v in row))

print(f"Held-Karp shortest tour length: {best_len:.6f}")
tour_1based = [c + 1 for c in tour]
print("Held-Karp shortest tour (city numbers, 1-based): " + " -> ".join(map(str, tour_1based)))

# --- Independent check: exhaustive brute-force over all permutations ---
# Fix city 0 as start; permute the remaining n-1 cities; take the closed-tour minimum.
brute_best, brute_tour = math.inf, None
for perm in permutations(others):
    length = d[0][perm[0]] + sum(d[perm[k]][perm[k + 1]] for k in range(len(perm) - 1)) + d[perm[-1]][0]
    if length < brute_best:
        brute_best, brute_tour = length, (0,) + perm + (0,)

print(f"Brute-force shortest tour length: {brute_best:.6f}")
print("Brute-force shortest tour (city numbers, 1-based): " +
      " -> ".join(str(c + 1) for c in brute_tour))
print(f"Held-Karp vs brute-force length difference: {abs(best_len - brute_best):.3e}")
print(f"Check passed (lengths match to 1e-6): {abs(best_len - brute_best) < 1e-6}")

# --- Plot the shortest tour through the ten cities ---
plt.figure(figsize=(7, 7))
tx = [x[c] for c in tour]
ty = [y[c] for c in tour]
plt.plot(tx, ty, "-o", color="tab:blue", zorder=1)
for c in range(n):
    plt.scatter(x[c], y[c], color="tab:red", zorder=2)
    plt.annotate(str(c + 1), (x[c], y[c]), textcoords="offset points", xytext=(6, 6))
plt.title(f"Shortest TSP tour (Held-Karp), length = {best_len:.2f}")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/9B.3.1_s3.png")

# The brute-force search evaluates every possible closed route, so its minimum is the
# true optimum by definition; matching it confirms the Held-Karp result is exact.
print("Explanation: brute force evaluates every possible closed route, so its minimum "
      "is the true optimum by definition, and Held-Karp matching it confirms exactness.")
