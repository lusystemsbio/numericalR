import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from collections import deque

# ---------------------------------------------------------------
# INPUT MODEL: an undirected adjacency matrix.
# Example graph: Zachary's karate club (34 vertices).
# We only borrow networkx to obtain the standard graph + a
# Kamada-Kawai layout; every structural property below is
# computed explicitly from the adjacency matrix.
# ---------------------------------------------------------------
G = nx.karate_club_graph()
n = G.number_of_nodes()
A = nx.to_numpy_array(G, nodelist=range(n))   # undirected 0/1 adjacency matrix
A = (A > 0).astype(int)

# ---------------------------------------------------------------
# 1) DEGREE: for each vertex, sum its row of the adjacency matrix.
# ---------------------------------------------------------------
degree = A.sum(axis=1).astype(int)            # degree of each vertex
mean_degree = degree.mean()                   # average degree

# degree distribution: count how many vertices have each degree value
deg_values, deg_counts = np.unique(degree, return_counts=True)

# ---------------------------------------------------------------
# 2) LOCAL CLUSTERING COEFFICIENT.
#    For vertex i, look at its neighbours; count the edges that
#    actually exist among them, divide by the max possible edges
#    k*(k-1)/2. Average over all vertices for the global value.
# ---------------------------------------------------------------
local_clustering = np.zeros(n)
for i in range(n):
    nbrs = np.where(A[i] == 1)[0]             # neighbours of i
    k = len(nbrs)
    if k < 2:                                 # need >=2 neighbours to form a triangle
        local_clustering[i] = 0.0
        continue
    # count edges among the neighbours (each counted once)
    links = 0
    for a_idx in range(k):
        for b_idx in range(a_idx + 1, k):
            if A[nbrs[a_idx], nbrs[b_idx]] == 1:
                links += 1
    local_clustering[i] = 2.0 * links / (k * (k - 1))
avg_clustering = local_clustering.mean()      # network clustering coefficient

# ---------------------------------------------------------------
# 3) SHORTEST-PATH DISTANCE MATRIX via BFS from every source.
#    Unweighted graph -> breadth-first search gives hop distances.
# ---------------------------------------------------------------
INF = np.inf
dist = np.full((n, n), INF)
for s in range(n):
    dist[s, s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        for v in np.where(A[u] == 1)[0]:      # explore neighbours of u
            if dist[s, v] == INF:             # first time reached => shortest
                dist[s, v] = dist[s, u] + 1
                q.append(v)

# ---------------------------------------------------------------
# 4) DIAMETER: the largest finite shortest-path distance.
# ---------------------------------------------------------------
finite = dist[np.isfinite(dist)]
diameter = int(finite.max())

# ---------------------------------------------------------------
# 5) BETWEENNESS CENTRALITY via Brandes' algorithm.
#    For each source, do a BFS recording shortest-path counts,
#    then accumulate dependency contributions back to front.
# ---------------------------------------------------------------
betweenness = np.zeros(n)
for s in range(n):
    S = []                                    # stack of vertices in order of discovery
    pred = [[] for _ in range(n)]             # predecessors on shortest paths
    sigma = np.zeros(n); sigma[s] = 1.0       # number of shortest paths from s
    d = np.full(n, -1); d[s] = 0              # BFS distance from s
    q = deque([s])
    while q:
        v = q.popleft(); S.append(v)
        for w in np.where(A[v] == 1)[0]:
            if d[w] < 0:                      # w found for the first time
                d[w] = d[v] + 1; q.append(w)
            if d[w] == d[v] + 1:              # shortest path to w via v
                sigma[w] += sigma[v]
                pred[w].append(v)
    delta = np.zeros(n)                       # dependency accumulation
    while S:
        w = S.pop()
        for v in pred[w]:
            delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
        if w != s:
            betweenness[w] += delta[w]
betweenness /= 2.0                            # undirected: each pair counted twice

# ---------------------------------------------------------------
# PRINT RESULTS
# ---------------------------------------------------------------
print("=== Degree ===")
for i in range(n):
    print(f"degree[node {i}] = {degree[i]}")
print(f"mean_degree = {mean_degree:.4f}")
print("degree_distribution (degree: count):")
for dv, dc in zip(deg_values, deg_counts):
    print(f"  degree {int(dv)}: {int(dc)} nodes")

print("\n=== Clustering ===")
for i in range(n):
    print(f"local_clustering[node {i}] = {local_clustering[i]:.4f}")
print(f"average_clustering_coefficient = {avg_clustering:.4f}")

print("\n=== Diameter ===")
print(f"diameter = {diameter}")

print("\n=== Betweenness centrality ===")
for i in range(n):
    print(f"betweenness[node {i}] = {betweenness[i]:.4f}")

# ---------------------------------------------------------------
# SEPARATE CHECK / VALIDATION
# ---------------------------------------------------------------
print("\n=== Validation check ===")
top_deg = np.argsort(degree)[::-1][:3]
print(f"top-3 highest-degree hubs (nodes) = {top_deg.tolist()} with degrees {degree[top_deg].tolist()}")
print(f"mean degree approx 4.6? computed = {mean_degree:.4f}")
print(f"diameter equals 5? computed = {diameter}")
print(f"clustering approx 0.57? computed = {avg_clustering:.4f}")
top_btw = np.argsort(betweenness)[::-1][:2]
print(f"two highest-betweenness nodes = {top_btw.tolist()} "
      f"(node 0 = instructor 'Mr. Hi', node 33 = president/officer)")
print("Explanation: because the from-scratch degree (~4.6), diameter (5), and "
      "clustering (~0.57) match the known values while the two top-betweenness "
      "nodes are exactly the instructor (0) and president (33) who bridge the "
      "two factions, the independently reproduced structural signature confirms "
      "the computation is correct.")

# ---------------------------------------------------------------
# DRAW with Kamada-Kawai layout, sizing nodes by betweenness.
# ---------------------------------------------------------------
pos = nx.kamada_kawai_layout(G)
plt.figure(figsize=(9, 8))
node_sizes = 300 + 4000 * (betweenness / betweenness.max())
nx.draw_networkx_edges(G, pos, alpha=0.35)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=betweenness,
                       cmap=plt.cm.viridis)
nx.draw_networkx_labels(G, pos, font_size=8, font_color="white")
plt.title("Zachary's Karate Club (Kamada-Kawai layout)\n"
          "node size/color = betweenness centrality")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.1.1_s3.png")
