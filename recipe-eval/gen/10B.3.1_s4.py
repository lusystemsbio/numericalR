import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 1. Build the model input: gene-expression data, 500 genes x 26 samples.
#    We simulate a PCA-style progression (as in Part 10A): the 26 samples
#    form a smooth "condition" axis, and genes belong to latent co-expression
#    modules that rise/fall along that axis at different phases.
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)          # seed 1 for reproducibility
n_genes, n_samples = 500, 26

# condition axis (the ordered progression across the 26 samples)
t = np.linspace(0, 2 * np.pi, n_samples)

# each gene gets a random phase -> defines which module it belongs to
phases = rng.uniform(0, 2 * np.pi, n_genes)
# true underlying signal: gene expression = cos(condition - gene_phase)
signal = np.cos(t[None, :] - phases[:, None])          # 500 x 26
noise = rng.normal(0, 0.35, size=(n_genes, n_samples)) # measurement noise
expr = signal + noise                                  # raw expression

# ----------------------------------------------------------------------
# 2. Z-score each gene across its 26 samples (per-gene standardization).
# ----------------------------------------------------------------------
gene_mean = expr.mean(axis=1, keepdims=True)
gene_std = expr.std(axis=1, ddof=0, keepdims=True)
Z = (expr - gene_mean) / gene_std                      # 500 x 26, z-scored

print("Z-scored data shape (genes x samples):", Z.shape)
print("Per-gene mean (should be ~0):", np.round(Z.mean(axis=1).mean(), 6))
print("Per-gene std  (should be ~1):", np.round(Z.std(axis=1).mean(), 6))

# ----------------------------------------------------------------------
# 3. K-means clustering of the z-scored genes, implemented explicitly.
#    Each gene is a point in 26-dimensional (sample) space.
# ----------------------------------------------------------------------
k = 5
max_iter = 100
np.random.seed(1)                        # seed 1 for the k-means init

# --- initialize centroids by picking k random distinct genes ---
init_idx = rng.choice(n_genes, size=k, replace=False)
centroids = Z[init_idx].copy()           # k x 26

labels = np.zeros(n_genes, dtype=int)
for it in range(max_iter):
    # --- assignment step: each gene -> nearest centroid (Euclidean) ---
    # distances: n_genes x k
    dists = np.linalg.norm(Z[:, None, :] - centroids[None, :, :], axis=2)
    new_labels = dists.argmin(axis=1)

    # --- stop if assignments no longer change ---
    if it > 0 and np.array_equal(new_labels, labels):
        labels = new_labels
        print(f"K-means converged after {it} iterations.")
        break
    labels = new_labels

    # --- update step: recompute each centroid as the mean of its genes ---
    for c in range(k):
        members = Z[labels == c]
        if len(members) > 0:
            centroids[c] = members.mean(axis=0)
        else:
            # re-seed an empty cluster to a random gene
            centroids[c] = Z[rng.integers(n_genes)]
else:
    print(f"K-means reached max_iter = {max_iter}.")

# report cluster sizes
for c in range(k):
    print(f"Cluster {c}: {(labels == c).sum()} genes")

# --- within-cluster sum of squares (a quality summary) ---
wcss = sum(((Z[labels == c] - centroids[c]) ** 2).sum() for c in range(k))
print("Total within-cluster sum of squares (WCSS):", round(float(wcss), 4))

# ----------------------------------------------------------------------
# 4. Order genes by cluster so co-expression modules become contiguous
#    blocks. Within a cluster, order clusters by the sample (condition)
#    at which the module peaks, so the PCA progression is visible.
# ----------------------------------------------------------------------
peak_sample = centroids.argmax(axis=1)            # condition of each module's peak
cluster_order = np.argsort(peak_sample)           # order modules along progression

row_order = np.concatenate([np.where(labels == c)[0] for c in cluster_order])
Z_ordered = Z[row_order]

# boundaries between clusters (for drawing separators)
sizes = [(labels == c).sum() for c in cluster_order]
boundaries = np.cumsum(sizes)[:-1]

print("Cluster peak sample index (per original cluster id):", peak_sample.tolist())
print("Cluster display order (by peak condition):", cluster_order.tolist())

# ----------------------------------------------------------------------
# 5. Cluster-ordered expression heatmap, color scale -3 to 3.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(Z_ordered, aspect="auto", cmap="RdBu_r",
               vmin=-3, vmax=3, interpolation="nearest")

# horizontal lines separating the clusters
for b in boundaries:
    ax.axhline(b - 0.5, color="black", linewidth=1.0)

ax.set_xlabel("Sample (condition progression)")
ax.set_ylabel("Genes (ordered by k-means cluster)")
ax.set_title(f"Cluster-ordered gene-expression heatmap (k={k})")
cbar = fig.colorbar(im, ax=ax, label="z-scored expression")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.3.1_s4.png",
            dpi=150)

# ----------------------------------------------------------------------
# 6. Check: confirm the heatmap splits genes into co-expression modules
#    that rise and fall together across conditions. We quantify this by
#    the mean within-cluster correlation between each gene and its cluster
#    centroid: high values mean genes in a cluster truly co-vary.
# ----------------------------------------------------------------------
mean_corrs = []
for c in range(k):
    members = Z[labels == c]
    cen = centroids[c]
    # correlation of each member gene's profile with its centroid profile
    corrs = [np.corrcoef(g, cen)[0, 1] for g in members]
    m = float(np.mean(corrs))
    mean_corrs.append(m)
    print(f"Cluster {c}: mean gene-to-centroid correlation = {round(m, 4)}")

print("Overall mean within-cluster correlation:", round(float(np.mean(mean_corrs)), 4))

# also: correlation between the module centroids should span the progression
# (adjacent modules positively correlated, opposite-phase modules negatively)
centroid_corr = np.corrcoef(centroids[cluster_order])
print("Module-centroid correlation matrix (display order):")
print(np.round(centroid_corr, 2))

# ----------------------------------------------------------------------
# Why the check confirms the result (one sentence):
# A high mean within-cluster gene-to-centroid correlation shows that each
# k-means cluster is a genuine co-expression module whose genes rise and
# fall together across the 26 conditions, so the contiguous colored blocks
# in the heatmap resolve the continuous PCA progression of Part 10A into
# discrete, phase-ordered gene groups.
# ----------------------------------------------------------------------
print("Check passed:" ,
      np.mean(mean_corrs) > 0.7,
      "- clusters are coherent co-expression modules tracing the progression.")
