import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
from itertools import permutations

# City coordinates (city 1 is index 0)
x = (0, -28.87, -14.66, -29.06, -36.04, -50.48, -50.59, -0.14, -21.50, -43.07)
y = (0, 0, 43.39, 43.22, 21.61, -7.37, 21.59, 28.73, -7.32, -14.55)
n = len(x)

# Build the Euclidean distance matrix d[i][j]
d = [[math.hypot(x[i] - x[j], y[i] - y[j]) for j in range(n)] for i in range(n)]

# --- Held-Karp dynamic programming ---
# We fix city 0 as start/end. States are (subset S of cities to still-visit-context, endpoint).
# g[(S, i)] = length of shortest path that starts at city 0, visits exactly the set S,
#             and ends at city i (with i in S). S is a bitmask over cities 1..n-1.

g = {}      # cost table
parent = {} # to reconstruct the tour

# Base case: paths from 0 directly to i (S contains only i)
for i in range(1, n):
    S = 1 << (i - 1)          # bit (i-1) marks that city i is included
    g[(S, i)] = d[0][i]       # go straight from start city 0 to i
    parent[(S, i)] = 0

# Recurrence: grow the subset one city at a time.
# g(x, S) = min over i in S\{x} of ( g(i, S\{x}) + d(i, x) )
for size in range(2, n):                          # subsets of increasing size
    for subset in permutations(range(1, n), size):
        subset = tuple(sorted(subset))
        if subset[0] != subset[0]:  # dummy, keep sorted-unique handling below
            pass
        # skip duplicates: only process each combination once
        if list(subset) != sorted(set(subset)):
            continue
        S = 0
        for c in subset:
            S |= 1 << (c - 1)
        for last in subset:                       # 'x': the endpoint city
            prev_S = S & ~(1 << (last - 1))       # S without x
            best = math.inf
            best_prev = None
            for i in subset:                      # candidate previous city i in S
                if i == last:
                    continue
                cost = g[(prev_S, i)] + d[i][last]
                if cost < best:
                    best = cost
                    best_prev = i
            g[(S, last)] = best
            parent[(S, last)] = best_prev

# Close the tour: return to city 0 from the last city, over the full set.
full = (1 << (n - 1)) - 1
best_len = math.inf
best_last = None
for i in range(1, n):
    cost = g[(full, i)] + d[i][0]
    if cost < best_len:
        best_len = cost
        best_last = i

# Reconstruct the optimal tour by walking parent pointers backward.
tour = [0]
S = full
last = best_last
path = []
while last != 0:
    path.append(last)
    plast = parent[(S, last)]
    S = S & ~(1 << (last - 1))
    last = plast
tour = [0] + path[::-1] + [0]

# --- Brute-force check: try every permutation of the other 9 cities ---
brute_len = math.inf
brute_tour = None
for perm in permutations(range(1, n)):
    route = (0,) + perm + (0,)
    length = sum(d[route[k]][route[k + 1]] for k in range(n))
    if length < brute_len:
        brute_len = length
        brute_tour = route

# Report results
print("Held-Karp optimal tour (city numbers, 1-indexed):",
      [c + 1 for c in tour])
print("Held-Karp optimal tour length:", best_len)
print("Brute-force optimal tour length:", brute_len)
print("Brute-force optimal tour (1-indexed):", [c + 1 for c in brute_tour])
print("Check (Held-Karp length matches brute-force length):",
      math.isclose(best_len, brute_len))
print("Absolute difference between the two lengths:", abs(best_len - brute_len))

# Explanation of the check:
print("Why the check confirms the result: the brute-force search evaluates every "
      "possible permutation of the cities, so its minimum is the exact optimum, and "
      "Held-Karp matching it proves Held-Karp found that same exact shortest tour.")

# --- Draw the shortest tour ---
fig, ax = plt.subplots(figsize=(7, 7))
tx = [x[c] for c in tour]
ty = [y[c] for c in tour]
ax.plot(tx, ty, '-o', color='steelblue', zorder=1)
for c in range(n):
    ax.scatter(x[c], y[c], color='crimson', zorder=2)
    ax.annotate(str(c + 1), (x[c], y[c]),
                textcoords="offset points", xytext=(6, 6), fontsize=11)
ax.set_title("Held-Karp shortest tour (length = %.4f)" % best_len)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
ax.grid(True, linestyle=":", alpha=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/9B.3.1_s1.png")
