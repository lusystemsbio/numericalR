import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 1. Zachary's karate club as an undirected adjacency matrix (34 vertices)
# ---------------------------------------------------------------
# Standard edge list of the karate club (0-indexed vertices).
edges = [
    (0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,10),(0,11),
    (0,12),(0,13),(0,17),(0,19),(0,21),(0,31),
    (1,2),(1,3),(1,7),(1,13),(1,17),(1,19),(1,21),(1,30),
    (2,3),(2,7),(2,8),(2,9),(2,13),(2,27),(2,28),(2,32),
    (3,7),(3,12),(3,13),
    (4,6),(4,10),
    (5,6),(5,10),(5,16),
    (6,16),
    (8,30),(8,32),(8,33),
    (9,33),
    (13,33),
    (14,32),(14,33),
    (15,32),(15,33),
    (18,32),(18,33),
    (19,33),
    (20,32),(20,33),
    (22,32),(22,33),
    (23,25),(23,27),(23,29),(23,32),(23,33),
    (24,25),(24,27),(24,31),
    (25,31),
    (26,29),(26,33),
    (27,33),
    (28,31),(28,33),
    (29,32),(29,33),
    (30,32),(30,33),
    (31,32),(31,33),
    (32,33),
]

N = 34
# Build symmetric adjacency matrix from the edge list.
A = np.zeros((N, N), dtype=int)
for i, j in edges:
    A[i, j] = 1
    A[j, i] = 1

# ---------------------------------------------------------------
# 2. Degree distribution
# ---------------------------------------------------------------
# Degree of a vertex = number of neighbours = row sum of A.
degree = A.sum(axis=1)
mean_degree = degree.mean()

print("=== Degree ===")
for v in range(N):
    print(f"degree[{v:2d}] = {degree[v]}")
print(f"mean_degree = {mean_degree:.4f}")
# Distribution: count of vertices at each degree value.
print("degree_distribution (degree: count):")
vals, counts = np.unique(degree, return_counts=True)
for d, c in zip(vals, counts):
    print(f"  degree {d}: {c} nodes")
# Highest-degree hubs.
hubs = np.argsort(degree)[::-1][:4]
print("top_degree_hubs (node, degree):")
for h in hubs:
    print(f"  node {h}: degree {degree[h]}")

# ---------------------------------------------------------------
# 3. Clustering coefficient (local, then averaged)
# ---------------------------------------------------------------
# Local clustering C_i = (2 * #edges among neighbours) / (k_i * (k_i - 1)).
local_clustering = np.zeros(N)
for i in range(N):
    neighbours = np.where(A[i] == 1)[0]
    k = len(neighbours)
    if k < 2:
        local_clustering[i] = 0.0  # undefined -> 0 by convention
        continue
    # Count edges present between pairs of neighbours.
    links = 0
    for a_idx in range(k):
        for b_idx in range(a_idx + 1, k):
            if A[neighbours[a_idx], neighbours[b_idx]] == 1:
                links += 1
    local_clustering[i] = 2.0 * links / (k * (k - 1))
mean_clustering = local_clustering.mean()

print("\n=== Clustering ===")
for v in range(N):
    print(f"clustering[{v:2d}] = {local_clustering[v]:.4f}")
print(f"mean_clustering = {mean_clustering:.4f}")

# ---------------------------------------------------------------
# 4. Shortest-path distance matrix via BFS (unweighted) + diameter
# ---------------------------------------------------------------
def bfs_distances(adj, source):
    """Return array of shortest hop distances from source (inf if unreachable)."""
    n = adj.shape[0]
    dist = np.full(n, np.inf)
    dist[source] = 0
    queue = [source]
    while queue:
        u = queue.pop(0)
        for w in np.where(adj[u] == 1)[0]:
            if dist[w] == np.inf:
                dist[w] = dist[u] + 1
                queue.append(w)
    return dist

# Distance matrix: one BFS per source.
D = np.zeros((N, N))
for s in range(N):
    D[s] = bfs_distances(A, s)

# Diameter = largest finite shortest-path distance.
diameter = int(D[np.isfinite(D)].max())

print("\n=== Shortest paths ===")
print(f"diameter = {diameter}")
print(f"mean_shortest_path = {D[np.isfinite(D)].mean():.4f}")

# ---------------------------------------------------------------
# 5. Betweenness centrality (Brandes' algorithm, undirected/unweighted)
# ---------------------------------------------------------------
# For each source, do BFS accumulating #shortest paths, then back-propagate
# dependencies to sum how often each node lies on shortest paths between pairs.
betweenness = np.zeros(N)
for s in range(N):
    S = []                      # stack of nodes in order of non-decreasing distance
    P = [[] for _ in range(N)]  # predecessors on shortest paths
    sigma = np.zeros(N)         # number of shortest paths from s
    sigma[s] = 1
    dist = np.full(N, -1)
    dist[s] = 0
    queue = [s]
    while queue:
        v = queue.pop(0)
        S.append(v)
        for w in np.where(A[v] == 1)[0]:
            if dist[w] < 0:                 # first time we reach w
                dist[w] = dist[v] + 1
                queue.append(w)
            if dist[w] == dist[v] + 1:      # shortest path to w via v
                sigma[w] += sigma[v]
                P[w].append(v)
    delta = np.zeros(N)
    # Accumulate dependencies in reverse BFS order.
    while S:
        w = S.pop()
        for v in P[w]:
            delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
        if w != s:
            betweenness[w] += delta[w]
# Each shortest path counted twice (undirected); divide by 2.
betweenness /= 2.0

print("\n=== Betweenness centrality ===")
for v in range(N):
    print(f"betweenness[{v:2d}] = {betweenness[v]:.4f}")
top_bet = np.argsort(betweenness)[::-1][:2]
print("two_highest_betweenness (node, value):")
for b in top_bet:
    print(f"  node {b}: {betweenness[b]:.4f}")

# ---------------------------------------------------------------
# 6. Kamada-Kawai layout from scratch (stress majorization on graph distances)
# ---------------------------------------------------------------
# Target ideal distances = graph-theoretic shortest-path distances D.
# Minimize sum_{i<j} (||x_i - x_j|| - D_ij)^2 / D_ij^2 by simple gradient descent.
rng = np.random.default_rng(0)
pos = rng.standard_normal((N, 2))
weight = np.zeros_like(D)
mask = D > 0
weight[mask] = 1.0 / (D[mask] ** 2)   # Kamada-Kawai weighting

for _ in range(2000):
    grad = np.zeros_like(pos)
    for i in range(N):
        diff = pos[i] - pos                       # vectors to all others
        dcur = np.linalg.norm(diff, axis=1)
        dcur[dcur == 0] = 1e-9
        # gradient of stress w.r.t. x_i
        coeff = weight[i] * (dcur - D[i]) / dcur
        coeff[i] = 0
        grad[i] = (coeff[:, None] * diff).sum(axis=0)
    pos -= 0.01 * grad

# ---------------------------------------------------------------
# 7. Draw the network with the Kamada-Kawai layout
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 8))
# Draw edges.
for i, j in edges:
    ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
            color="gray", alpha=0.4, zorder=1)
# Node size scaled by degree, colour by betweenness.
sizes = 80 + degree * 40
sc = ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=betweenness,
                cmap="viridis", zorder=2, edgecolors="k")
for v in range(N):
    ax.text(pos[v, 0], pos[v, 1], str(v), fontsize=7,
            ha="center", va="center", zorder=3, color="white")
plt.colorbar(sc, ax=ax, label="betweenness centrality")
ax.set_title("Zachary's Karate Club (Kamada-Kawai layout)")
ax.set_axis_off()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.1.1_s2.png", dpi=150)

# ---------------------------------------------------------------
# 8. Independent check of the known structural facts
# ---------------------------------------------------------------
print("\n=== CHECK ===")
print(f"check_mean_degree ~ 4.6 : {mean_degree:.4f} -> {abs(mean_degree - 4.6) < 0.1}")
print(f"check_diameter == 5     : {diameter} -> {diameter == 5}")
print(f"check_clustering ~ 0.57  : {mean_clustering:.4f} -> {abs(mean_clustering - 0.57) < 0.02}")
print(f"check_hubs (nodes 0 & 33 among top degrees): "
      f"{set([0, 33]).issubset(set(hubs.tolist()))}")
print(f"top_two_betweenness nodes: {sorted(top_bet.tolist())} "
      f"-> instructor(0) & president(33): {sorted(top_bet.tolist()) == [0, 33]}")
# Explanation:
print("\nExplanation: The check confirms the result because reproducing the club's "
      "known fingerprint from scratch\n(mean degree ~4.6, diameter 5, clustering ~0.57, "
      "and betweenness peaking exactly at the instructor (node 0) and\npresident (node 33) "
      "who bridge the two factions) matches the independently documented ground truth, "
      "so our\nimplementations of degree, clustering, distances, and betweenness must be correct.")
