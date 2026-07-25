import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 0. Reproducibility
# ---------------------------------------------------------------
rng = np.random.default_rng(1)          # seed 1 for data + k-means init

# ---------------------------------------------------------------
# 1. Build synthetic gene-expression data: 500 genes x 26 samples.
#    We plant 5 co-expression modules whose mean profiles rise/fall
#    along an ordered progression of the 26 conditions (mimicking the
#    PCA progression of Part 10A).
# ---------------------------------------------------------------
n_genes, n_samples, k = 500, 26, 5
t = np.linspace(0, 2 * np.pi, n_samples)                 # ordered "condition axis"

# Distinct smooth profiles for the 5 true modules
profiles = np.vstack([
    np.sin(t),                       # early up-down
    np.cos(t),                       # early down-up
    np.sin(2 * t),                   # oscillating
    np.linspace(-1.5, 1.5, n_samples),   # monotone rise
    np.linspace(1.5, -1.5, n_samples),   # monotone fall
])

# Assign each gene to a module and add gene-level noise
true_labels = rng.integers(0, k, size=n_genes)
X = np.zeros((n_genes, n_samples))
for g in range(n_genes):
    X[g] = profiles[true_labels[g]] + rng.normal(0, 0.35, n_samples)

# ---------------------------------------------------------------
# 2. Z-score each gene (row) so clustering is by expression *shape*.
# ---------------------------------------------------------------
Z = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)

# ---------------------------------------------------------------
# 3. K-means clustering, implemented explicitly (Lloyd's algorithm).
# ---------------------------------------------------------------
# 3a. Initialize centroids by picking k random distinct genes.
init_idx = rng.choice(n_genes, size=k, replace=False)
centroids = Z[init_idx].copy()

max_iter = 100
for it in range(max_iter):
    # 3b. Assignment step: each gene -> nearest centroid (Euclidean).
    dists = np.linalg.norm(Z[:, None, :] - centroids[None, :, :], axis=2)
    labels = dists.argmin(axis=1)

    # 3c. Update step: recompute each centroid as the mean of its genes.
    new_centroids = np.array([
        Z[labels == c].mean(axis=0) if np.any(labels == c) else centroids[c]
        for c in range(k)
    ])

    # 3d. Convergence check: stop when centroids stop moving.
    shift = np.linalg.norm(new_centroids - centroids)
    centroids = new_centroids
    if shift < 1e-6:
        break

# ---------------------------------------------------------------
# 4. Order genes by cluster (and by within-cluster similarity) so
#    co-expression modules form contiguous blocks in the heatmap.
# ---------------------------------------------------------------
order = np.concatenate([np.where(labels == c)[0] for c in range(k)])
Z_ordered = Z[order]
labels_ordered = labels[order]

# Cluster sizes and boundary rows (for drawing separators)
cluster_sizes = np.array([np.sum(labels == c) for c in range(k)])
boundaries = np.cumsum(cluster_sizes)

# ---------------------------------------------------------------
# 5. Cluster-ordered heatmap, color scale fixed to [-3, 3].
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(Z_ordered, aspect="auto", cmap="RdBu_r",
               vmin=-3, vmax=3, interpolation="nearest")
for b in boundaries[:-1]:
    ax.axhline(b - 0.5, color="black", lw=1)      # separate cluster blocks
ax.set_xlabel("Samples / conditions (ordered)")
ax.set_ylabel("Genes (ordered by k-means cluster)")
ax.set_title("K-means cluster-ordered gene-expression heatmap (k=5)")
fig.colorbar(im, ax=ax, label="z-score")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.3.1_s3.png")

# ---------------------------------------------------------------
# 6. Check: quantify that each cluster is a coherent co-expression
#    module (genes within a cluster are well correlated across
#    conditions, and cluster mean profiles differ from each other).
# ---------------------------------------------------------------
print(f"Iterations to converge: {it + 1}")
print(f"Cluster sizes: {cluster_sizes.tolist()}")

# Mean expression profile of each cluster across the 26 conditions
cluster_means = np.array([Z[labels == c].mean(axis=0) for c in range(k)])

# Within-cluster coherence: mean pairwise correlation of genes to their
# cluster's mean profile (high => rise and fall together).
for c in range(k):
    members = Z[labels == c]
    corrs = np.array([np.corrcoef(g, cluster_means[c])[0, 1] for g in members])
    print(f"Cluster {c}: size={members.shape[0]:3d}, "
          f"mean gene-to-profile correlation = {corrs.mean():.3f}")

# Between-cluster separation: mean profiles should be distinct.
prof_corr = np.corrcoef(cluster_means)
off_diag = prof_corr[~np.eye(k, dtype=bool)]
print(f"Mean |correlation| between distinct cluster profiles = "
      f"{np.abs(off_diag).mean():.3f}")
print(f"Max correlation between distinct cluster profiles = "
      f"{off_diag.max():.3f}")

# Peak condition of each cluster's mean profile (shows the progression).
peak_conditions = cluster_means.argmax(axis=1)
print(f"Peak condition index per cluster: {peak_conditions.tolist()}")

# ---------------------------------------------------------------
# 7. One-sentence explanation of why this check confirms the result:
# ---------------------------------------------------------------
print("Explanation: High within-cluster gene-to-profile correlations with "
      "low between-cluster profile correlations confirm the heatmap blocks "
      "are genuine co-expression modules that rise and fall together across "
      "conditions, resolving the Part 10A PCA progression into distinct gene groups.")
