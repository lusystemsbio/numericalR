import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# ---------------------------------------------------------------
# Input model: an undirected adjacency matrix A (A[i,j]=1 if edge).
# We take Zachary's karate club (34 vertices) as the example graph,
# but everything structural below is computed FROM the matrix by hand.
# networkx is used ONLY to source the graph and for layout/drawing.
# ---------------------------------------------------------------
G_src = nx.karate_club_graph()
n = G_src.number_of_nodes()
A = np.zeros((n, n), dtype=int)
for i, j in G_src.edges():
    A[i, j] = 1
    A[j, i] = 1  # symmetric: undirected

# ---------------------------------------------------------------
# 1) Degree distribution: degree of a node = row sum of A.
# ---------------------------------------------------------------
degrees = A.sum(axis=1)                    # number of neighbors per node
mean_degree = degrees.mean()
# tally how many nodes have each degree value
deg_vals, deg_counts = np.unique(degrees, return_counts=True)

# ---------------------------------------------------------------
# 2) Clustering coefficient (local): for each node, fraction of
#    possible edges among its neighbors that actually exist.
#    C_i = 2*(#links among neighbors) / (k_i*(k_i-1)).
# ---------------------------------------------------------------
clustering = np.zeros(n)
for i in range(n):
    nbrs = np.where(A[i] == 1)[0]          # neighbor indices
    k = len(nbrs)
    if k < 2:
        clustering[i] = 0.0                # undefined -> 0 by convention
        continue
    links = 0
    for a_idx in range(k):                 # count edges between neighbor pairs
        for b_idx in range(a_idx + 1, k):
            if A[nbrs[a_idx], nbrs[b_idx]] == 1:
                links += 1
    clustering[i] = 2.0 * links / (k * (k - 1))
avg_clustering = clustering.mean()

# ---------------------------------------------------------------
# 3) Shortest-path distance matrix via BFS from every node
#    (unweighted graph -> BFS gives geodesic distances).
# ---------------------------------------------------------------
INF = np.iinfo(np.int32).max
dist = np.full((n, n), INF, dtype=int)
for s in range(n):
    dist[s, s] = 0
    queue = [s]
    while queue:
        u = queue.pop(0)                   # FIFO queue -> breadth-first
        for v in np.where(A[u] == 1)[0]:
            if dist[s, v] == INF:          # first time reached = shortest
                dist[s, v] = dist[s, u] + 1
                queue.append(v)

# ---------------------------------------------------------------
# 4) Diameter: the largest finite shortest-path distance.
# ---------------------------------------------------------------
finite = dist[dist < INF]
diameter = int(finite.max())

# ---------------------------------------------------------------
# 5) Betweenness centrality via Brandes' algorithm.
#    For each source s, one BFS accumulates the number of shortest
#    paths (sigma) and back-propagates dependencies (delta).
# ---------------------------------------------------------------
betweenness = np.zeros(n)
for s in range(n):
    S = []                                 # stack: nodes in order of removal
    P = [[] for _ in range(n)]             # predecessors on shortest paths
    sigma = np.zeros(n); sigma[s] = 1.0    # # of shortest paths s->v
    d = np.full(n, -1, dtype=int); d[s] = 0
    queue = [s]
    while queue:                           # BFS from s
        v = queue.pop(0)
        S.append(v)
        for w in np.where(A[v] == 1)[0]:
            if d[w] < 0:                   # w found for the first time
                queue.append(w)
                d[w] = d[v] + 1
            if d[w] == d[v] + 1:           # shortest edge on a geodesic
                sigma[w] += sigma[v]
                P[w].append(v)
    delta = np.zeros(n)                    # dependency accumulation
    while S:                               # process in reverse BFS order
        w = S.pop()
        for v in P[w]:
            delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
        if w != s:
            betweenness[w] += delta[w]
betweenness /= 2.0                         # each pair counted twice (undirected)

# ---------------------------------------------------------------
# Print all computed results.
# ---------------------------------------------------------------
print("=== Degree ===")
for i in range(n):
    print(f"node {i:2d} degree: {degrees[i]}")
for dv, dc in zip(deg_vals, deg_counts):
    print(f"degree {dv:2d}: {dc} nodes")
print(f"mean degree: {mean_degree:.4f}")
hubs = np.argsort(degrees)[::-1][:5]
print(f"top-5 highest-degree hubs (node:degree): "
      + ", ".join(f"{h}:{degrees[h]}" for h in hubs))

print("\n=== Clustering ===")
for i in range(n):
    print(f"node {i:2d} clustering: {clustering[i]:.4f}")
print(f"average clustering coefficient: {avg_clustering:.4f}")

print("\n=== Diameter ===")
print(f"diameter: {diameter}")

print("\n=== Betweenness centrality ===")
for i in range(n):
    print(f"node {i:2d} betweenness: {betweenness[i]:.4f}")
top2 = np.argsort(betweenness)[::-1][:2]
print(f"two highest-betweenness nodes: {int(top2[0])} and {int(top2[1])}")

# ---------------------------------------------------------------
# Separate independent CHECK against known ground truth values.
# ---------------------------------------------------------------
print("\n=== CHECK ===")
check_mean = abs(mean_degree - 4.6) < 0.1
check_diam = (diameter == 5)
check_clust = abs(avg_clustering - 0.57) < 0.02
# In nx.karate_club_graph, node 0 = instructor "Mr. Hi", node 33 = "Officer"/president
check_bridge = set(int(x) for x in top2) == {0, 33}
print(f"has high-degree hubs (node0={degrees[0]}, node33={degrees[33]}): "
      f"{degrees[0] > 10 and degrees[33] > 10}")
print(f"mean degree ~4.6: {check_mean} (got {mean_degree:.4f})")
print(f"diameter == 5: {check_diam} (got {diameter})")
print(f"clustering ~0.57: {check_clust} (got {avg_clustering:.4f})")
print(f"top-2 betweenness == instructor(0) & president(33): {check_bridge} "
      f"(got {int(top2[0])}, {int(top2[1])})")
print(f"ALL CHECKS PASS: {check_mean and check_diam and check_clust and check_bridge}")
# One-sentence explanation:
print("Explanation: matching the known mean degree, diameter, clustering, and the "
      "instructor/president as top bridges reproduces Zachary's documented statistics, "
      "which confirms the from-scratch structural computations are correct.")

# ---------------------------------------------------------------
# Draw the network with a Kamada-Kawai layout; size/color by degree,
# highlight the two highest-betweenness bridge nodes.
# ---------------------------------------------------------------
pos = nx.kamada_kawai_layout(G_src)
plt.figure(figsize=(10, 8))
node_colors = [betweenness[i] for i in range(n)]
node_sizes = [80 + 40 * degrees[i] for i in range(n)]
nx.draw_networkx_edges(G_src, pos, alpha=0.3)
nodes = nx.draw_networkx_nodes(G_src, pos, node_size=node_sizes,
                               node_color=node_colors, cmap="viridis")
nx.draw_networkx_labels(G_src, pos, font_size=8)
# outline the two faction-bridging leaders
nx.draw_networkx_nodes(G_src, pos, nodelist=list(int(x) for x in top2),
                       node_size=[node_sizes[i] for i in top2],
                       node_color="none", edgecolors="red", linewidths=2.5)
plt.colorbar(nodes, label="betweenness centrality")
plt.title("Zachary's Karate Club (Kamada-Kawai layout)\n"
          "size ~ degree, color ~ betweenness, red = instructor & president")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.1.1_s1.png")
