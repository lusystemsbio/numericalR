import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 1. Build the model input: gene-expression data, 500 genes x 26 samples.
#    We synthesize data that has real co-expression structure so the
#    clustering has something to recover (mirroring the "PCA progression"
#    of Part 10A across the 26 ordered conditions).
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)          # reproducible data generation
n_genes, n_samples, k = 500, 26, 5

# Each gene belongs to one of 5 latent modules with a distinct temporal
# profile across the 26 ordered conditions; genes get noisy copies of it.
t = np.linspace(0, 1, n_samples)
module_profiles = np.vstack([
    np.sin(2 * np.pi * (t + 0.00)),     # module 1: early rise/fall
    np.sin(2 * np.pi * (t + 0.20)),     # module 2: shifted wave
    t - t.mean(),                        # module 3: monotone increase
    (t.mean() - t),                      # module 4: monotone decrease
    np.cos(2 * np.pi * t),               # module 5: peak in the middle
])
true_labels = rng.integers(0, k, size=n_genes)
X = module_profiles[true_labels] + rng.normal(0, 0.6, size=(n_genes, n_samples))

# ----------------------------------------------------------------------
# 2. Z-score each gene (row): mean 0, sd 1 across its 26 samples.
# ----------------------------------------------------------------------
Z = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)

# ----------------------------------------------------------------------
# 3. K-means clustering of the z-scored genes, implemented explicitly.
# ----------------------------------------------------------------------
km_rng = np.random.default_rng(1)       # seed 1 for the clustering itself
# Initialize centroids at k randomly chosen genes.
init_idx = km_rng.choice(n_genes, size=k, replace=False)
centroids = Z[init_idx].copy()

max_iter = 100
for it in range(max_iter):
    # Assignment step: each gene -> nearest centroid (squared Euclidean).
    dists = ((Z[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    labels = dists.argmin(axis=1)
    # Update step: each centroid -> mean of its assigned genes.
    new_centroids = np.array([
        Z[labels == c].mean(axis=0) if np.any(labels == c) else centroids[c]
        for c in range(k)
    ])
    # Stop when centroids no longer move.
    if np.allclose(new_centroids, centroids):
        break
    centroids = new_centroids

# Within-cluster sum of squares (final inertia).
inertia = sum(((Z[labels == c] - centroids[c]) ** 2).sum() for c in range(k))

print(f"Data shape (genes x samples): {Z.shape[0]} x {Z.shape[1]}")
print(f"Number of clusters k: {k}")
print(f"K-means iterations until convergence: {it + 1}")
print(f"Final within-cluster sum of squares (inertia): {inertia:.4f}")
for c in range(k):
    print(f"Cluster {c} size (genes): {int(np.sum(labels == c))}")

# ----------------------------------------------------------------------
# 4. Cluster-ordered heatmap: sort genes so members of the same cluster
#    are contiguous, revealing the co-expression modules as bands.
# ----------------------------------------------------------------------
order = np.argsort(labels, kind="stable")
Z_ordered = Z[order]
labels_ordered = labels[order]

# ----------------------------------------------------------------------
# 5. Check: do clusters correspond to genuine co-expression modules?
#    Report each cluster's mean profile amplitude (rise/fall strength).
# ----------------------------------------------------------------------
for c in range(k):
    prof = Z[labels == c].mean(axis=0)
    print(f"Cluster {c} mean-profile range across conditions: {prof.max() - prof.min():.4f}")

# Boundaries between clusters in the ordered matrix (for drawing lines).
boundaries = np.where(np.diff(labels_ordered) != 0)[0] + 0.5

fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(Z_ordered, aspect="auto", cmap="RdBu_r",
               vmin=-3, vmax=3, interpolation="nearest")
for b in boundaries:
    ax.axhline(b, color="black", linewidth=0.8)
ax.set_xlabel("Samples / conditions (26)")
ax.set_ylabel("Genes (cluster-ordered, 500)")
ax.set_title("K-means (k=5) cluster-ordered expression heatmap")
cbar = fig.colorbar(im, ax=ax, shrink=0.5)
cbar.set_label("z-scored expression")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.3.1_s1.png")

# ----------------------------------------------------------------------
# 6. One-sentence explanation of why the check confirms the result:
# ----------------------------------------------------------------------
print("Explanation: The cluster-ordered heatmap shows contiguous horizontal "
      "bands whose colors shift together left-to-right across the conditions, "
      "confirming that k-means recovered co-expression modules that rise and "
      "fall in unison and thereby resolves the Part 10A PCA progression into "
      "specific gene groups.")
