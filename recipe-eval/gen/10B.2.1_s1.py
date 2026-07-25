import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy.cluster.hierarchy import dendrogram

# ---------------------------------------------------------------
# 1. The three-blob data (same input as the k-means test)
# ---------------------------------------------------------------
X, y_true = make_blobs(n_samples=150, centers=3, cluster_std=0.70,
                       random_state=42)
n = X.shape[0]

# ---------------------------------------------------------------
# 2. Pairwise Euclidean distances between all points
# ---------------------------------------------------------------
def euclidean(a, b):
    return np.sqrt(np.sum((a - b) ** 2))

# active[i] holds the current inter-cluster distance dict; we start
# with every point as its own singleton cluster (ids 0..n-1).
dist = {}                                   # dist[(i,j)] with i<j
for i in range(n):
    for j in range(i + 1, n):
        dist[(i, j)] = euclidean(X[i], X[j])

active = list(range(n))                      # currently living cluster ids
size = {i: 1 for i in range(n)}              # number of points per cluster

# ---------------------------------------------------------------
# 3. Agglomerative Ward clustering, done explicitly (no one-step call)
#    Ward's Lance-Williams recurrence updates distances after each merge:
#    d(ij,k) = sqrt( ((ni+nk)dik^2 + (nj+nk)djk^2 - nk*dij^2) / (ni+nj+nk) )
# ---------------------------------------------------------------
def get_d(a, b):
    return dist[(a, b)] if a < b else dist[(b, a)]

Z = []                                       # scipy-style linkage matrix
next_id = n                                  # id given to each newly merged cluster
merge_history = []                           # (idA, idB) actually merged, in order

for _ in range(n - 1):
    # 3a. find the closest pair among the currently active clusters
    a, b, best = None, None, np.inf
    for ii in range(len(active)):
        for jj in range(ii + 1, len(active)):
            ca, cb = active[ii], active[jj]
            d = get_d(ca, cb)
            if d < best:
                best, a, b = d, ca, cb

    # 3b. record this merge (id_a, id_b, distance, combined size)
    new_size = size[a] + size[b]
    Z.append([a, b, best, new_size])
    merge_history.append((a, b))

    # 3c. Ward update: distance from the new cluster to every other cluster
    dab = best
    for c in active:
        if c == a or c == b:
            continue
        dik, djk = get_d(a, c), get_d(b, c)
        ni, nj, nk = size[a], size[b], size[c]
        newd = np.sqrt(((ni + nk) * dik**2 + (nj + nk) * djk**2
                        - nk * dab**2) / (ni + nj + nk))
        lo, hi = (next_id, c) if next_id < c else (c, next_id)
        dist[(lo, hi)] = newd

    # 3d. retire a and b, activate the new cluster
    active.remove(a); active.remove(b); active.append(next_id)
    size[next_id] = new_size
    next_id += 1

Z = np.array(Z)

# ---------------------------------------------------------------
# 4. Cut the tree into k = 3 flat clusters by performing only the
#    first (n - k) merges via a union-find, then reading components.
# ---------------------------------------------------------------
k = 3
parent = list(range(2 * n - 1))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

for step in range(n - k):                    # stop k-1 merges early -> k clusters
    a, b = merge_history[step]
    parent[a] = n + step                     # both children point to their merged id
    parent[b] = n + step

roots = [find(i) for i in range(n)]          # root id of each original point
uniq = {r: idx for idx, r in enumerate(sorted(set(roots)))}
labels = np.array([uniq[r] for r in roots])  # relabel to 0,1,2

# ---------------------------------------------------------------
# 5. k-means reference partition for agreement check
# ---------------------------------------------------------------
km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(X)
labels_km = km.labels_

# ---------------------------------------------------------------
# 6. Check: are there three clear branches? Look at the merge heights.
#    The last two merges (joining the 3 blobs) should tower over the rest.
# ---------------------------------------------------------------
heights = Z[:, 2]
top3 = heights[-3:]                          # heights of the final three merges
gap = top3[2] - top3[1]                      # jump between 3rd-last and 2nd-last

ari_true = adjusted_rand_score(y_true, labels)
ari_km   = adjusted_rand_score(labels_km, labels)

# ---------------------------------------------------------------
# 7. Plots: dendrogram + points colored by the k=3 cut
# ---------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
dendrogram(Z, ax=ax[0], color_threshold=top3[1] + 1e-9, no_labels=True)
ax[0].set_title("Ward dendrogram")
ax[0].set_xlabel("points"); ax[0].set_ylabel("merge height")
ax[1].scatter(X[:, 0], X[:, 1], c=labels, cmap="viridis", s=25)
ax[1].set_title("Points colored by 3-cluster cut")
ax[1].set_xlabel("x1"); ax[1].set_ylabel("x2")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.2.1_s1.png")

# ---------------------------------------------------------------
# 8. Report every numerical result
# ---------------------------------------------------------------
print("Number of points:", n)
print("Number of clusters requested (k):", k)
print("Cluster sizes from Ward cut:", list(np.bincount(labels)))
print("Cluster sizes from k-means:", list(np.bincount(labels_km)))
print("Height of 3rd-from-last merge (largest within-branch join):", top3[0])
print("Height of 2nd-from-last merge (joins two blobs):", top3[1])
print("Height of last merge (joins final branch):", top3[2])
print("Gap between last and 2nd-last merge heights:", gap)
print("Ratio last-merge-height / 3rd-from-last-merge-height:", top3[2] / top3[0])
print("Adjusted Rand index, Ward cut vs true blobs:", ari_true)
print("Adjusted Rand index, Ward cut vs k-means:", ari_km)
print("Three clear branches (last two heights >> the rest):", bool(top3[1] / top3[0] > 1.5))
print("Cut recovers blobs and agrees with k-means:", bool(ari_true > 0.9 and ari_km > 0.9))
# The check confirms the result because a large jump before the final two merges
# means the tree splits cleanly into three high-level branches, and an adjusted
# Rand index near 1 against both the true blobs and k-means shows the 3-way cut
# assigns points to the same groups as the known partition.
print("Explanation: a big height gap before the top two merges plus an adjusted "
      "Rand index near 1 versus the true labels and k-means shows the three "
      "branches are the three blobs, so the hierarchical cut recovers them.")
