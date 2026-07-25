import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# 1. Build Zachary's karate club adjacency matrix (34 vertices, undirected)
# ---------------------------------------------------------------------------
# Edge list of the classic Zachary karate club network (1-indexed nodes).
edges = [
    (1,2),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(1,9),(1,11),(1,12),(1,13),
    (1,14),(1,18),(1,20),(1,22),(1,32),(2,3),(2,4),(2,8),(2,14),(2,18),
    (2,20),(2,22),(2,31),(3,4),(3,8),(3,9),(3,10),(3,14),(3,28),(3,29),
    (3,33),(4,8),(4,13),(4,14),(5,7),(5,11),(6,7),(6,11),(6,17),(7,17),
    (9,31),(9,33),(9,34),(10,34),(14,34),(15,33),(15,34),(16,33),(16,34),
    (19,33),(19,34),(20,34),(21,33),(21,34),(23,33),(23,34),(24,26),(24,28),
    (24,30),(24,33),(24,34),(25,26),(25,28),(25,32),(26,32),(27,30),(27,34),
    (28,34),(29,32),(29,34),(30,33),(30,34),(31,33),(31,34),(32,33),(32,34),
    (33,34),
]

n = 34
A = np.zeros((n, n), dtype=int)
for u, v in edges:
    A[u-1, v-1] = 1        # convert to 0-indexed and make symmetric
    A[v-1, u-1] = 1

# ---------------------------------------------------------------------------
# 2. Degree distribution  (degree = number of neighbors = row sum)
# ---------------------------------------------------------------------------
degrees = A.sum(axis=1)
mean_degree = degrees.mean()

print("=== Degree ===")
for i in range(n):
    print(f"Node {i+1:2d} degree: {degrees[i]}")
print(f"Mean degree: {mean_degree:.4f}")
# Identify hubs: the highest-degree vertices
hub_order = np.argsort(degrees)[::-1]
print(f"Top-5 highest-degree nodes (hubs): "
      f"{[(int(i+1), int(degrees[i])) for i in hub_order[:5]]}")

# ---------------------------------------------------------------------------
# 3. Local clustering coefficient (fraction of a node's neighbor pairs that
#    are themselves connected), then the network average.
# ---------------------------------------------------------------------------
clustering = np.zeros(n)
for i in range(n):
    nbrs = np.where(A[i] == 1)[0]      # neighbors of i
    k = len(nbrs)
    if k < 2:
        clustering[i] = 0.0            # undefined -> 0 for degree < 2
        continue
    # count edges present among the neighbors
    links = 0
    for a in range(k):
        for b in range(a+1, k):
            if A[nbrs[a], nbrs[b]] == 1:
                links += 1
    # possible pairs = k*(k-1)/2
    clustering[i] = 2.0 * links / (k * (k - 1))

mean_clustering = clustering.mean()
print("\n=== Clustering ===")
for i in range(n):
    print(f"Node {i+1:2d} clustering: {clustering[i]:.4f}")
print(f"Mean clustering coefficient: {mean_clustering:.4f}")

# ---------------------------------------------------------------------------
# 4. Shortest-path distance matrix via BFS from every source (unweighted).
# ---------------------------------------------------------------------------
INF = np.inf
dist = np.full((n, n), INF)

def bfs(src):
    """Return array of shortest-path hop counts from src to all nodes."""
    d = np.full(n, INF)
    d[src] = 0
    queue = [src]
    while queue:
        cur = queue.pop(0)
        for nb in np.where(A[cur] == 1)[0]:
            if d[nb] == INF:          # first time reached = shortest
                d[nb] = d[cur] + 1
                queue.append(nb)
    return d

for s in range(n):
    dist[s] = bfs(s)

# Diameter = largest finite shortest-path distance
diameter = int(dist[np.isfinite(dist)].max())
print("\n=== Diameter ===")
print(f"Diameter (longest shortest path): {diameter}")

# ---------------------------------------------------------------------------
# 5. Betweenness centrality via Brandes' algorithm (from scratch, unweighted).
#    Counts fraction of all shortest paths passing through each node.
# ---------------------------------------------------------------------------
betweenness = np.zeros(n)
for s in range(n):
    # single-source shortest-path counting
    S = []                              # stack of nodes in order of finishing
    P = [[] for _ in range(n)]          # predecessors on shortest paths
    sigma = np.zeros(n)                 # number of shortest paths s->v
    d = np.full(n, -1)                  # distance s->v (-1 = unvisited)
    sigma[s] = 1
    d[s] = 0
    queue = [s]
    while queue:
        v = queue.pop(0)
        S.append(v)
        for w in np.where(A[v] == 1)[0]:
            if d[w] < 0:                # found w for the first time
                d[w] = d[v] + 1
                queue.append(w)
            if d[w] == d[v] + 1:        # shortest path to w via v
                sigma[w] += sigma[v]
                P[w].append(v)
    # accumulation: back-propagate dependencies
    delta = np.zeros(n)
    while S:
        w = S.pop()
        for v in P[w]:
            delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
        if w != s:
            betweenness[w] += delta[w]

# For undirected graph each shortest path counted twice -> divide by 2
betweenness /= 2.0

print("\n=== Betweenness centrality ===")
for i in range(n):
    print(f"Node {i+1:2d} betweenness: {betweenness[i]:.4f}")
btw_order = np.argsort(betweenness)[::-1]
print(f"Two highest-betweenness nodes: "
      f"{[(int(i+1), round(float(betweenness[i]),4)) for i in btw_order[:2]]}")
print("(Node 1 = instructor 'Mr. Hi', Node 34 = president/officer)")

# ---------------------------------------------------------------------------
# 6. Draw the network with a Kamada-Kawai layout (implemented from scratch).
#    Kamada-Kawai places nodes so Euclidean distance ~ graph-theoretic
#    distance, by minimizing a spring energy; we do gradient descent.
# ---------------------------------------------------------------------------
rng = np.random.default_rng(0)
# target distances L = shortest-path distances; ideal spring lengths
L = dist.copy()
# desired lengths and stiffness (Kamada-Kawai standard: k = 1/L^2)
with np.errstate(divide='ignore'):
    K = 1.0 / (L ** 2)
np.fill_diagonal(K, 0.0)               # no self-interaction
K[~np.isfinite(K)] = 0.0
np.fill_diagonal(L, 0.0)

# initialize on a circle
theta = np.linspace(0, 2*np.pi, n, endpoint=False)
pos = np.column_stack([np.cos(theta), np.sin(theta)]) * n / 4.0

# gradient descent on Kamada-Kawai energy: sum_{i<j} 0.5*K_ij*(|xi-xj|-L_ij)^2
lr = 0.01
for it in range(2000):
    grad = np.zeros_like(pos)
    for i in range(n):
        diff = pos[i] - pos                      # vectors i->others
        d = np.sqrt((diff ** 2).sum(axis=1))     # current distances
        d[d == 0] = 1e-9
        coeff = K[i] * (1 - L[i] / d)            # per-node force coefficient
        grad[i] = (coeff[:, None] * diff).sum(axis=0)
    pos -= lr * grad

print("\n=== Layout ===")
print("Kamada-Kawai layout computed via gradient descent (from scratch).")

# plot
fig, ax = plt.subplots(figsize=(9, 9))
# draw edges
for u, v in edges:
    x = [pos[u-1, 0], pos[v-1, 0]]
    y = [pos[u-1, 1], pos[v-1, 1]]
    ax.plot(x, y, color="0.7", lw=0.8, zorder=1)
# node size ~ degree, color ~ betweenness
sizes = 80 + degrees * 40
sc = ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=betweenness,
                cmap="viridis", zorder=2, edgecolors="k")
for i in range(n):
    ax.text(pos[i, 0], pos[i, 1], str(i+1), fontsize=7,
            ha="center", va="center", color="white", zorder=3)
plt.colorbar(sc, ax=ax, label="Betweenness centrality")
ax.set_title("Zachary's Karate Club (Kamada-Kawai layout)\n"
             "node size ~ degree, color ~ betweenness")
ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.1.1_s5.png", dpi=130)

# ---------------------------------------------------------------------------
# 7. Separate verification check against known values.
# ---------------------------------------------------------------------------
print("\n=== Verification check ===")
hubs_ok = degrees.max() >= 15 and (degrees >= 12).sum() >= 2
mean_ok = abs(mean_degree - 4.6) < 0.2
diam_ok = (diameter == 5)
clust_ok = abs(mean_clustering - 0.57) < 0.03
top2 = sorted(int(i+1) for i in btw_order[:2])
bridge_ok = (top2 == [1, 34])

print(f"Has a few high-degree hubs (max deg {int(degrees.max())}): {hubs_ok}")
print(f"Mean degree ~4.6 (got {mean_degree:.3f}): {mean_ok}")
print(f"Diameter == 5 (got {diameter}): {diam_ok}")
print(f"Clustering ~0.57 (got {mean_clustering:.3f}): {clust_ok}")
print(f"Two highest-betweenness are nodes 1 & 34 (got {top2}): {bridge_ok}")
print(f"ALL CHECKS PASS: {all([hubs_ok, mean_ok, diam_ok, clust_ok, bridge_ok])}")
print("Explanation: matching the independently-known karate-club statistics "
      "(mean degree ~4.6, diameter 5, clustering ~0.57) and recovering nodes 1 "
      "and 34 as the top bridges confirms our from-scratch routines reproduce "
      "the true network structure and correctly identify the instructor and "
      "president who link the two factions.")
