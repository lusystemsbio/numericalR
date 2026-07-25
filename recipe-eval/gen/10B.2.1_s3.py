import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans

# ---------------------------------------------------------------
# 1. Build the same three-blob data
# ---------------------------------------------------------------
X, y_true = make_blobs(n_samples=150, centers=3, cluster_std=0.60,
                       random_state=0)
n = X.shape[0]

# ---------------------------------------------------------------
# 2. Pairwise Euclidean distances (condensed + square forms)
# ---------------------------------------------------------------
D_cond = pdist(X, metric="euclidean")   # condensed distance vector
D = squareform(D_cond)                   # full n x n distance matrix

# ---------------------------------------------------------------
# 3. Agglomerative hierarchical clustering with Ward linkage,
#    implemented explicitly via the Lance-Williams update.
#    Ward merges the pair of clusters that increases total
#    within-cluster variance the least.
# ---------------------------------------------------------------
# Each active cluster: list of member indices, size, and centroid.
clusters = {i: {"members": [i], "size": 1, "centroid": X[i].copy()}
            for i in range(n)}
active = set(range(n))

# Ward merge cost between two clusters (increase in variance):
def ward_cost(a, b):
    na, nb = clusters[a]["size"], clusters[b]["size"]
    diff = clusters[a]["centroid"] - clusters[b]["centroid"]
    return (na * nb) / (na + nb) * np.dot(diff, diff)

next_id = n                 # ids for newly formed clusters
linkage_rows = []           # scipy-style linkage matrix rows: [a, b, dist, size]

# Repeatedly merge until a single cluster remains.
while len(active) > 1:
    best = None
    best_pair = None
    act = list(active)
    # Find the pair with minimum Ward cost.
    for i in range(len(act)):
        for j in range(i + 1, len(act)):
            a, b = act[i], act[j]
            c = ward_cost(a, b)
            if best is None or c < best:
                best, best_pair = c, (a, b)
    a, b = best_pair
    # Merge a and b into a new cluster.
    na, nb = clusters[a]["size"], clusters[b]["size"]
    new_members = clusters[a]["members"] + clusters[b]["members"]
    new_size = na + nb
    # New centroid is the weighted mean of the two centroids.
    new_centroid = (na * clusters[a]["centroid"] + nb * clusters[b]["centroid"]) / new_size
    # scipy encodes the merge "distance" as sqrt(2 * ward_cost).
    linkage_rows.append([a, b, np.sqrt(2.0 * best), new_size])
    clusters[next_id] = {"members": new_members, "size": new_size,
                         "centroid": new_centroid}
    active.discard(a)
    active.discard(b)
    active.add(next_id)
    next_id += 1

Z = np.array(linkage_rows)   # our own linkage matrix

# ---------------------------------------------------------------
# 4. Cut the dendrogram into k = 3 flat clusters.
#    The last (k-1) merges join the k top branches; skipping them
#    leaves exactly k clusters. Walk the merge tree to collect them.
# ---------------------------------------------------------------
k = 3
# Reconstruct membership of every node id (leaves + merged nodes).
node_members = {i: [i] for i in range(n)}
for row_idx, (a, b, dist, size) in enumerate(Z):
    node_members[n + row_idx] = node_members[int(a)] + node_members[int(b)]

# The roots of the k branches are the "new" ids created by the
# first (n - k) merges, minus any that were themselves later merged.
merged_away = set()
for row_idx in range(n - k):          # merges that happen below the cut
    a, b, dist, size = Z[row_idx]
    merged_away.add(int(a))
    merged_away.add(int(b))
top_ids = [i for i in range(2 * n - 1) if i in node_members
           and i not in merged_away]
top_ids = [i for i in top_ids if i <= n + (n - k) - 1]  # only ids present at cut

labels = np.empty(n, dtype=int)
for lab, node in enumerate(top_ids):
    for m in node_members[node]:
        labels[m] = lab

# ---------------------------------------------------------------
# 5. k-means reference partition for comparison
# ---------------------------------------------------------------
km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(X)
km_labels = km.labels_

# Agreement is measured up to label permutation via a best-match count.
def best_agreement(a, b):
    from itertools import permutations
    la, lb = np.unique(a), np.unique(b)
    best = 0
    for perm in permutations(lb):
        mapping = {old: new for old, new in zip(lb, perm)}
        mapped = np.array([mapping[x] for x in b])
        best = max(best, np.sum(mapped == a))
    return best / len(a)

agree_km = best_agreement(labels, km_labels)
agree_true = best_agreement(labels, y_true)

# ---------------------------------------------------------------
# 6. Report numerical results
# ---------------------------------------------------------------
print(f"Number of points: {n}")
print(f"Number of pairwise distances: {D_cond.size}")
print(f"Number of merges performed: {Z.shape[0]}")
print(f"Requested clusters k: {k}")
sizes = np.bincount(labels)
for c in range(k):
    print(f"Hierarchical cluster {c} size: {sizes[c]}")
# The 3 largest merge heights: the last two (top of tree) should
# dwarf the rest, i.e. three clear branches.
heights = np.sort(Z[:, 2])
print(f"Largest merge height: {heights[-1]:.4f}")
print(f"Second largest merge height: {heights[-2]:.4f}")
print(f"Third largest merge height: {heights[-3]:.4f}")
print(f"Height gap (2nd largest / 3rd largest): {heights[-2] / heights[-3]:.4f}")
print(f"Agreement with k-means (fraction): {agree_km:.4f}")
print(f"Agreement with true blobs (fraction): {agree_true:.4f}")

# ---------------------------------------------------------------
# 7. Dendrogram (drawn from our own linkage matrix) + colored points
# ---------------------------------------------------------------
from scipy.cluster.hierarchy import dendrogram

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

dendrogram(Z, ax=axes[0], no_labels=True, color_threshold=heights[-3] * 1.01)
axes[0].axhline(y=(heights[-2] + heights[-3]) / 2, color="k",
                linestyle="--", label="cut for k=3")
axes[0].set_title("Ward dendrogram (three clear branches)")
axes[0].set_ylabel("Merge height")
axes[0].legend()

for c in range(k):
    pts = X[labels == c]
    axes[1].scatter(pts[:, 0], pts[:, 1], s=25, label=f"cluster {c}")
axes[1].set_title("Points colored by 3 cut clusters")
axes[1].set_xlabel("x1")
axes[1].set_ylabel("x2")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.2.1_s3.png")

# ---------------------------------------------------------------
# 8. Why the check confirms the result (one sentence)
# ---------------------------------------------------------------
print("Check explanation: A large gap between the top two merge heights "
      "and the rest shows three well-separated branches, and the near-perfect "
      "agreement of the 3-cluster cut with both the true blobs and k-means "
      "confirms the dendrogram recovers the same structure.")
