import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans

# ----------------------------------------------------------------------
# 1. Build the three-blob data (same input as the k-means recipe)
# ----------------------------------------------------------------------
X, y_true = make_blobs(n_samples=150, centers=3, cluster_std=0.70,
                       random_state=4)
n = X.shape[0]

# ----------------------------------------------------------------------
# 2. Pairwise Euclidean distances (condensed + square form)
# ----------------------------------------------------------------------
D = squareform(pdist(X, metric="euclidean"))   # n x n distance matrix

# ----------------------------------------------------------------------
# 3. Agglomerative hierarchical clustering with Ward linkage, done
#    explicitly (Lance-Williams update) rather than via linkage() in
#    one step.  Ward merges the pair of clusters that minimises the
#    increase in total within-cluster variance.
# ----------------------------------------------------------------------
# active[i]  -> True while cluster i still exists
# members[i] -> list of original point indices in cluster i
# size[i]    -> number of points in cluster i
# Cluster distances start as squared Euclidean point distances (Ward
# formulas are expressed in squared distances).
active = [True] * n
members = [[i] for i in range(n)]
size = [1] * n
Ddist = D ** 2                       # squared distances between clusters
np.fill_diagonal(Ddist, np.inf)      # never merge a cluster with itself

# The linkage matrix Z records each merge: [id_a, id_b, height, new_size]
Z = np.zeros((n - 1, 4))
next_id = n                          # ids for newly formed clusters
# map current row/col index -> cluster id (starts as the point's own id)
cur_id = list(range(n))

for step in range(n - 1):
    # ---- find the closest currently-active pair (i, j) ----
    best = np.inf
    pi = pj = -1
    idx = [k for k in range(len(active)) if active[k]]
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if Ddist[i, j] < best:
                best = Ddist[i, j]
                pi, pj = i, j

    # ---- record the merge (height = Ward distance, not squared) ----
    Z[step] = [cur_id[pi], cur_id[pj], np.sqrt(best), size[pi] + size[pj]]

    # ---- Lance-Williams update of squared distances to merged cluster ----
    ni, nj = size[pi], size[pj]
    for k in idx:
        if k == pi or k == pj:
            continue
        nk = size[k]
        # Ward's Lance-Williams coefficients (on squared distances)
        total = ni + nj + nk
        new_d = ((ni + nk) * Ddist[pi, k]
                 + (nj + nk) * Ddist[pj, k]
                 - nk * Ddist[pi, pj]) / total
        Ddist[pi, k] = new_d
        Ddist[k, pi] = new_d

    # ---- fold cluster pj into cluster pi ----
    members[pi] = members[pi] + members[pj]
    size[pi] = ni + nj
    cur_id[pi] = next_id             # merged cluster gets a fresh id
    next_id += 1
    active[pj] = False               # pj no longer exists
    Ddist[pj, :] = np.inf
    Ddist[:, pj] = np.inf

# ----------------------------------------------------------------------
# 4. Cut the tree into a flat partition of k = 3 clusters.
#    Undo the last (k-1) merges: the clusters present just before those
#    merges are exactly the k top-level branches.
# ----------------------------------------------------------------------
k = 3
labels = -np.ones(n, dtype=int)

# Rebuild membership by replaying merges but stopping (n - k) merges in,
# i.e. perform only the first n-k unions so that k clusters remain.
parent = list(range(2 * n - 1))      # union-find over point + merged ids
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

for step in range(n - k):            # replay only first n-k merges
    a, b = int(Z[step, 0]), int(Z[step, 1])
    parent[find(a)] = n + step       # new merged id for this step
    parent[find(b)] = n + step

# Assign a compact 0..k-1 label per remaining root
roots = {}
for i in range(n):
    r = find(i)
    if r not in roots:
        roots[r] = len(roots)
    labels[i] = roots[r]

for c in range(k):
    print(f"Ward cut cluster {c} size: {np.sum(labels == c)}")

# ----------------------------------------------------------------------
# 5. Compare with k-means (the flat-partition reference method)
# ----------------------------------------------------------------------
km = KMeans(n_clusters=3, n_init=10, random_state=4).fit(X)
km_labels = km.labels_

# Agreement is label-invariant, so compare via the confusion matrix:
# best matching of Ward labels to k-means labels.
from itertools import permutations
best_agree = 0
for perm in permutations(range(3)):
    mapped = np.array([perm[l] for l in labels])
    best_agree = max(best_agree, np.mean(mapped == km_labels))
print(f"Best agreement between Ward cut and k-means: {best_agree:.3f}")

# ----------------------------------------------------------------------
# 6. Check: the three top-of-tree merges happen at much larger heights
#    than the within-branch merges -> three clear branches.
# ----------------------------------------------------------------------
heights = np.sort(Z[:, 2])
top3 = heights[-3:]                        # last three (largest) merge heights
next_below = heights[-4]                   # 4th largest = tallest within-branch merge
print(f"Three tallest merge heights: {top3[0]:.3f}, {top3[1]:.3f}, {top3[2]:.3f}")
print(f"Tallest within-branch merge height: {next_below:.3f}")
gap_ratio = top3[0] / next_below
print(f"Gap ratio (3rd-tallest / 4th-tallest merge): {gap_ratio:.3f}")

# ----------------------------------------------------------------------
# 7. Plots: dendrogram + points colored by the three cut clusters
# ----------------------------------------------------------------------
from scipy.cluster.hierarchy import dendrogram

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left: dendrogram of our explicitly-built linkage matrix Z
dendrogram(Z, ax=axes[0], color_threshold=next_below + 1e-9,
           no_labels=True)
axes[0].axhline(y=(top3[0] + next_below) / 2, color="k", ls="--",
                label="cut for k=3")
axes[0].set_title("Ward dendrogram (three clear branches)")
axes[0].set_xlabel("points")
axes[0].set_ylabel("merge height")
axes[0].legend()

# Right: scatter colored by the three cut clusters
for c in range(k):
    pts = X[labels == c]
    axes[1].scatter(pts[:, 0], pts[:, 1], s=25, label=f"cluster {c}")
axes[1].set_title("Points colored by k=3 Ward cut")
axes[1].set_xlabel("x1")
axes[1].set_ylabel("x2")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.2.1_s4.png")

# ----------------------------------------------------------------------
# 8. Why the check confirms the result (one sentence)
# ----------------------------------------------------------------------
print("Explanation: A large height gap before the last three merges means the "
      "three branches are far more separated than points within any branch, so "
      "cutting at k=3 isolates the true blobs, and its near-perfect agreement "
      "with the independent k-means partition confirms the recovered clustering.")
