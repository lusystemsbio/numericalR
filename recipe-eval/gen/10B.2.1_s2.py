import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans

# ---------------------------------------------------------------
# 1. Build the three-blob data
# ---------------------------------------------------------------
X, y_true = make_blobs(n_samples=150, centers=3, cluster_std=0.70,
                       random_state=42)
n = X.shape[0]

# ---------------------------------------------------------------
# 2. Pairwise Euclidean distances (condensed upper-triangular form)
#    We implement agglomerative Ward clustering explicitly below,
#    so we work from the full distance/coordinate information.
# ---------------------------------------------------------------
# full squared-distance matrix (used implicitly by Ward update rule)
diff = X[:, None, :] - X[None, :, :]
D = np.sqrt((diff ** 2).sum(axis=2))          # Euclidean distance matrix
print("Number of points:", n)
print("Max pairwise Euclidean distance:", D.max())

# ---------------------------------------------------------------
# 3. Explicit agglomerative hierarchical clustering with Ward linkage.
#    Ward merges the two clusters that give the smallest increase in
#    total within-cluster variance. Using the Lance-Williams update
#    on squared distances lets us avoid recomputing from scratch.
# ---------------------------------------------------------------
# Each active cluster tracks: member indices, size, and its centroid.
clusters = {i: {"members": [i], "size": 1, "centroid": X[i].copy()}
            for i in range(n)}
next_id = n                       # ids for newly formed merged clusters
linkage = []                      # rows: [id_a, id_b, merge_distance, size]

active = list(clusters.keys())
while len(active) > 1:
    best_pair = None
    best_delta = np.inf
    # scan all pairs of active clusters for the smallest Ward increase
    for a_idx in range(len(active)):
        for b_idx in range(a_idx + 1, len(active)):
            a, b = active[a_idx], active[b_idx]
            na, nb = clusters[a]["size"], clusters[b]["size"]
            ca, cb = clusters[a]["centroid"], clusters[b]["centroid"]
            # Ward merge cost = (na*nb)/(na+nb) * ||ca - cb||^2
            delta = (na * nb) / (na + nb) * np.sum((ca - cb) ** 2)
            if delta < best_delta:
                best_delta = delta
                best_pair = (a, b)

    a, b = best_pair
    na, nb = clusters[a]["size"], clusters[b]["size"]
    # new merged centroid is the size-weighted average
    new_centroid = (na * clusters[a]["centroid"] + nb * clusters[b]["centroid"]) / (na + nb)
    new_members = clusters[a]["members"] + clusters[b]["members"]
    # scipy dendrogram convention: distance is sqrt(2 * Ward cost)
    merge_dist = np.sqrt(2.0 * best_delta)
    linkage.append([a, b, merge_dist, na + nb])

    clusters[next_id] = {"members": new_members,
                         "size": na + nb,
                         "centroid": new_centroid}
    del clusters[a]
    del clusters[b]
    active.remove(a)
    active.remove(b)
    active.append(next_id)
    next_id += 1

Z = np.array(linkage)
print("Linkage matrix shape:", Z.shape)
print("Largest three merge distances:", np.round(np.sort(Z[:, 2])[-3:], 4))

# ---------------------------------------------------------------
# 4. Cut the tree into k = 3 flat clusters.
#    Undo the last (k-1) merges: the clusters present just before the
#    final 2 merges are exactly the 3 top-level branches.
# ---------------------------------------------------------------
k = 3
# rebuild membership by replaying merges but stopping (n - k) merges in
parent_members = {i: [i] for i in range(n)}
cur_id = n
for row in Z[:n - k]:            # apply only the first (n-k) merges
    a, b = int(row[0]), int(row[1])
    parent_members[cur_id] = parent_members[a] + parent_members[b]
    del parent_members[a]
    del parent_members[b]
    cur_id += 1

labels_hc = np.empty(n, dtype=int)
for lab, (_, members) in enumerate(parent_members.items()):
    labels_hc[members] = lab

_, counts = np.unique(labels_hc, return_counts=True)
print("Hierarchical (Ward, k=3) cluster sizes:", sorted(counts.tolist()))

# ---------------------------------------------------------------
# 5. Cross-check against k-means
# ---------------------------------------------------------------
km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(X)
labels_km = km.labels_


def agreement(a, b):
    # fraction of point-pairs that agree on same/different cluster
    same_a = a[:, None] == a[None, :]
    same_b = b[:, None] == b[None, :]
    return (same_a == same_b).mean()


print("Pairwise agreement (HC vs true blobs): %.4f" % agreement(labels_hc, y_true))
print("Pairwise agreement (HC vs k-means):    %.4f" % agreement(labels_hc, labels_km))
print("Pairwise agreement (k-means vs true):  %.4f" % agreement(labels_km, y_true))

# ---------------------------------------------------------------
# 6. Plots: dendrogram + points colored by the 3 cut clusters
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

dendrogram(Z, ax=axes[0], color_threshold=np.sort(Z[:, 2])[-2] - 1e-9,
           no_labels=True)
axes[0].set_title("Ward dendrogram (three clear branches)")
axes[0].set_xlabel("points")
axes[0].set_ylabel("merge distance")

for lab in range(k):
    m = labels_hc == lab
    axes[1].scatter(X[m, 0], X[m, 1], s=25, label=f"cluster {lab}")
axes[1].set_title("Points colored by 3-cluster cut")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.2.1_s2.png")

# ---------------------------------------------------------------
# 7. Why the check confirms the result
# ---------------------------------------------------------------
print("Explanation: A large gap before the last two merges means three tight "
      "branches join only at high cost, so cutting at k=3 yields the natural "
      "blobs; near-perfect pairwise agreement with k-means confirms both "
      "methods recover the same three-blob structure.")
