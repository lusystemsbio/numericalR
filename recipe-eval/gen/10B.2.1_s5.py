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
X, y_true = make_blobs(n_samples=150, centers=3, cluster_std=1.0, random_state=5)
n = X.shape[0]
print("Number of points:", n)
print("True number of blobs:", len(np.unique(y_true)))

# ---------------------------------------------------------------
# 2. Pairwise Euclidean distance matrix (the model input)
# ---------------------------------------------------------------
diff = X[:, None, :] - X[None, :, :]          # all pairwise coordinate differences
D = np.sqrt((diff ** 2).sum(axis=2))          # Euclidean distances
print("Distance matrix shape:", D.shape)
print("Max pairwise distance:", D.max())

# ---------------------------------------------------------------
# 3. Agglomerative hierarchical clustering with Ward linkage,
#    implemented explicitly (Lance-Williams update rule).
#    Ward merges the pair of clusters that gives the smallest
#    increase in total within-cluster variance.
# ---------------------------------------------------------------
# active clusters: each starts as a singleton
clusters = {i: [i] for i in range(n)}         # cluster id -> list of member point indices
sizes = {i: 1 for i in range(n)}              # cluster sizes
# store squared distances between clusters (Ward works on squared distances)
Dsq = {(i, j): D[i, j] ** 2 for i in range(n) for j in range(i + 1, n)}

def key(a, b):
    return (a, b) if a < b else (b, a)

next_id = n                                   # ids for newly formed clusters
Z = np.zeros((n - 1, 4))                      # SciPy-style linkage matrix

for step in range(n - 1):
    # 3a. find the closest pair of active clusters (min Ward squared distance)
    (a, b), best = min(Dsq.items(), key=lambda kv: kv[1])
    # 3b. record the merge: distance stored is the actual (non-squared) Ward distance
    Z[step] = [a, b, np.sqrt(best), sizes[a] + sizes[b]]
    # 3c. create the new merged cluster
    new = next_id
    next_id += 1
    na, nb = sizes[a], sizes[b]
    # 3d. Lance-Williams update: distance from merged cluster to every other cluster c
    for c in list(clusters):
        if c == a or c == b:
            continue
        nc = sizes[c]
        d_ac = Dsq[key(a, c)]
        d_bc = Dsq[key(b, c)]
        d_ab = best
        # Ward's Lance-Williams coefficients (on squared distances)
        total = na + nb + nc
        d_new = ((na + nc) * d_ac + (nb + nc) * d_bc - nc * d_ab) / total
        Dsq[key(new, c)] = d_new
    # 3e. remove the two merged clusters from the active set
    for c in list(clusters):
        Dsq.pop(key(a, c), None)
        Dsq.pop(key(b, c), None)
    clusters[new] = clusters[a] + clusters[b]
    sizes[new] = na + nb
    del clusters[a], clusters[b], sizes[a], sizes[b]

print("Linkage matrix shape:", Z.shape)
print("Largest merge distances (last 5):", np.round(Z[-5:, 2], 3).tolist())

# ---------------------------------------------------------------
# 4. Cut the dendrogram into k = 3 flat clusters.
#    The first (n - k) merges are kept; the remaining connected
#    components are the k clusters.
# ---------------------------------------------------------------
k = 3
parent = {i: i for i in range(next_id)}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
# apply only the first n - k merges (i.e. stop before the last k-1 joins)
new_id = n
for step in range(n - k):
    a, b = int(Z[step, 0]), int(Z[step, 1])
    parent[a] = new_id
    parent[b] = new_id
    new_id += 1
# label each original point by the root of its component
roots = {}
labels_hc = np.empty(n, dtype=int)
for i in range(n):
    r = find(i)
    if r not in roots:
        roots[r] = len(roots)
    labels_hc[i] = roots[r]
print("Hierarchical cut cluster sizes:", np.bincount(labels_hc).tolist())

# ---------------------------------------------------------------
# 5. Separate check: compare against k-means with k = 3
# ---------------------------------------------------------------
labels_km = KMeans(n_clusters=3, n_init=10, random_state=0).fit_predict(X)

# The gap between the k=3 cut level and the next merge tells us how
# "clear" the three branches are: a large gap => three well-separated groups.
merge_at_3 = Z[n - k - 1, 2]   # distance of the last merge that still leaves 3 clusters
merge_to_2 = Z[n - k, 2]       # distance of the merge that would join two of the 3 branches
print("Merge distance leaving 3 clusters:", round(merge_at_3, 3))
print("Merge distance dropping to 2 clusters:", round(merge_to_2, 3))
print("Gap ratio (higher = clearer 3 branches):", round(merge_to_2 / merge_at_3, 3))

# Agreement between the hierarchical cut and k-means, invariant to label naming:
def agreement(a, b):
    # fraction of point pairs that are grouped the same way by both labelings
    same_a = a[:, None] == a[None, :]
    same_b = b[:, None] == b[None, :]
    return (same_a == same_b).mean()

agree = agreement(labels_hc, labels_km)
print("Pairwise agreement (Rand index) with k-means:", round(agree, 4))
print("Number of clusters recovered:", len(np.unique(labels_hc)))

# ---------------------------------------------------------------
# 6. Plots: dendrogram + points colored by the 3 cut clusters
# ---------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
dendrogram(Z, ax=ax[0], color_threshold=merge_to_2 - 1e-9, no_labels=True)
ax[0].axhline(y=(merge_at_3 + merge_to_2) / 2, color="k", ls="--",
              label="cut for k = 3")
ax[0].set_title("Ward dendrogram (three clear branches)")
ax[0].set_ylabel("Merge distance")
ax[0].legend()

for c in range(k):
    m = labels_hc == c
    ax[1].scatter(X[m, 0], X[m, 1], s=25, label=f"cluster {c}")
ax[1].set_title("Points colored by the 3-cluster cut")
ax[1].set_xlabel("x1"); ax[1].set_ylabel("x2")
ax[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.2.1_s5.png")

# ---------------------------------------------------------------
# 7. One-sentence explanation of the check
# ---------------------------------------------------------------
print("Explanation: A large gap between the merge distances at k=3 and k=2 shows "
      "the three branches are well separated, and a near-1 pairwise agreement with "
      "k-means confirms the cut recovers the same three blobs, so both independent "
      "methods point to the same true structure.")
