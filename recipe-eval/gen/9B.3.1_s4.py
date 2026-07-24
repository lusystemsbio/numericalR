import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
from itertools import permutations

# --- City coordinates (city 1 is index 0) ---
x = (0, -28.87, -14.66, -29.06, -36.04, -50.48, -50.59, -0.14, -21.50, -43.07)
y = (0, 0, 43.39, 43.22, 21.61, -7.37, 21.59, 28.73, -7.32, -14.55)
n = len(x)

# --- Build the Euclidean distance matrix d[i][j] ---
d = [[math.hypot(x[i] - x[j], y[i] - y[j]) for j in range(n)] for i in range(n)]

# =====================================================================
# Held-Karp dynamic program (exact TSP)
# We fix city 0 as start/end. States are (endpoint x, set S of intermediate
# cities visited, not counting the start city 0). g(x, S) = shortest path
# from start 0, through exactly the cities in S, ending at x.
# Recurrence: g(x, S) = min over i in S of ( g(i, S\{i}) + d[i][x] ).
# We represent S as a bitmask over cities 1..n-1.
# =====================================================================

# g maps (endpoint, mask) -> cost of the best path 0 -> ... -> endpoint
g = {}
# parent stores the predecessor endpoint, for reconstructing the tour
parent = {}

# Base case: paths from 0 directly to a single city k (S = {k}).
for k in range(1, n):
    mask = 1 << (k - 1)          # bit for city k
    g[(k, mask)] = d[0][k]       # cost 0 -> k
    parent[(k, mask)] = 0

# Build up over subsets of increasing size (subsets of cities 1..n-1).
full = (1 << (n - 1)) - 1        # mask with all intermediate cities set
for size in range(2, n):
    for subset in range(1, full + 1):
        # only consider subsets with exactly `size` cities
        if bin(subset).count("1") != size:
            continue
        for k in range(1, n):
            bit_k = 1 << (k - 1)
            if not (subset & bit_k):     # k must be the current endpoint, in subset
                continue
            prev = subset ^ bit_k        # subset without k
            best_cost = math.inf
            best_prev = -1
            # try every possible predecessor i in prev
            for i in range(1, n):
                if not (prev & (1 << (i - 1))):
                    continue
                cand = g[(i, prev)] + d[i][k]
                if cand < best_cost:
                    best_cost = cand
                    best_prev = i
            g[(k, subset)] = best_cost
            parent[(k, subset)] = best_prev

# Close the tour: return from the last city k back to start 0.
best_len = math.inf
last = -1
for k in range(1, n):
    cost = g[(k, full)] + d[k][0]
    if cost < best_len:
        best_len = cost
        last = k

# Reconstruct the ordered tour by walking parents backwards.
tour = []
mask = full
k = last
while k != 0:
    tour.append(k)
    pk = parent[(k, mask)]
    mask ^= (1 << (k - 1))
    k = pk
tour.append(0)
tour.reverse()          # now starts at 0
tour_full = tour + [0]  # close the loop back to city 1

# Report using 1-based city labels.
tour_labels = [c + 1 for c in tour_full]
print("Held-Karp shortest tour (city labels):", tour_labels)
print("Held-Karp shortest tour length:", best_len)

# =====================================================================
# Independent check: brute-force over all permutations of cities 1..n-1.
# Because n = 10 is small, we can enumerate every possible tour and take
# the minimum, giving a ground-truth exact answer.
# =====================================================================
bf_len = math.inf
bf_tour = None
for perm in permutations(range(1, n)):
    route = (0,) + perm + (0,)
    length = sum(d[route[i]][route[i + 1]] for i in range(len(route) - 1))
    if length < bf_len:
        bf_len = length
        bf_tour = route

print("Brute-force shortest tour (city labels):", [c + 1 for c in bf_tour])
print("Brute-force shortest tour length:", bf_len)
print("Held-Karp matches brute force (length):", math.isclose(best_len, bf_len))

# One sentence explaining why the check confirms the result:
# Because the brute-force search enumerates every possible tour and keeps
# the minimum, it is by definition the exact optimum, so Held-Karp agreeing
# with it confirms Held-Karp also found the exact shortest tour.
print("Why the check confirms it: brute force enumerates all tours and keeps the minimum,")
print("so it is the exact optimum by definition; Held-Karp matching it proves Held-Karp is exact too.")

# --- Draw the shortest tour ---
fig, ax = plt.subplots(figsize=(8, 7))
tx = [x[c] for c in tour_full]
ty = [y[c] for c in tour_full]
ax.plot(tx, ty, "-o", color="tab:blue", zorder=1)
for c in range(n):
    ax.scatter(x[c], y[c], color="tab:red", zorder=2)
    ax.annotate(str(c + 1), (x[c], y[c]), textcoords="offset points",
                xytext=(6, 6), fontsize=11, fontweight="bold")
ax.set_title("Held-Karp shortest TSP tour (length = {:.4f})".format(best_len))
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal", adjustable="datalim")
ax.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/9B.3.1_s4.png")
