import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 0. Build the model input: gene-expression data (500 genes x 26 samples).
#    The 26 samples represent an ordered progression of conditions (as in
#    the PCA of Part 10A). We synthesize several latent "programs" that
#    rise and fall along that progression, then assign each gene to one
#    program plus noise -- this gives real co-expression modules to recover.
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)          # reproducible synthetic data
n_genes, n_samples = 500, 26
t = np.linspace(0.0, 1.0, n_samples)    # the ordered condition axis

# Five smooth latent programs across the progression (rise/fall shapes).
programs = np.vstack([
    np.sin(np.pi * t),                       # up then down (mid peak)
    t,                                       # monotonic rise
    1.0 - t,                                 # monotonic fall
    np.sin(2 * np.pi * t),                   # early up, late down (biphasic)
    np.exp(-((t - 0.75) ** 2) / 0.02),       # late sharp burst
])
n_programs = programs.shape[0]

# Assign each gene to one program; expression = program signal + noise.
gene_program = rng.integers(0, n_programs, size=n_genes)
amplitude = rng.uniform(1.5, 3.0, size=n_genes)[:, None]
baseline = rng.uniform(-1.0, 1.0, size=n_genes)[:, None]
noise = rng.normal(0.0, 0.6, size=(n_genes, n_samples))
X = baseline + amplitude * programs[gene_program] + noise   # (500 x 26)

# ----------------------------------------------------------------------
# 1. Z-score each gene across samples (mean 0, sd 1 per row).
# ----------------------------------------------------------------------
mu = X.mean(axis=1, keepdims=True)
sd = X.std(axis=1, keepdims=True)
sd[sd == 0] = 1.0                       # guard against constant genes
Z = (X - mu) / sd
print(f"Z-scored matrix shape: {Z.shape[0]} genes x {Z.shape[1]} samples")
print(f"Per-gene mean (should be ~0): {np.abs(Z.mean(axis=1)).max():.3e} max abs")
print(f"Per-gene sd   (should be ~1): {Z.std(axis=1).mean():.4f} average")

# ----------------------------------------------------------------------
# 2. K-means clustering of the z-scored genes, implemented explicitly.
# ----------------------------------------------------------------------
k = 5
seed_rng = np.random.default_rng(1)     # seed 1 for the clustering itself

# Initialize centroids by picking k distinct genes at random.
init_idx = seed_rng.choice(n_genes, size=k, replace=False)
centroids = Z[init_idx].copy()

labels = np.zeros(n_genes, dtype=int)
for iteration in range(100):            # Lloyd's algorithm iterations
    # --- Assignment step: each gene -> nearest centroid (Euclidean). ---
    # dists[g, c] = squared distance from gene g to centroid c.
    dists = ((Z[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    new_labels = dists.argmin(axis=1)

    # --- Update step: recompute each centroid as its members' mean. ---
    new_centroids = np.zeros_like(centroids)
    for c in range(k):
        members = Z[new_labels == c]
        if len(members) > 0:
            new_centroids[c] = members.mean(axis=0)
        else:                           # re-seed an empty cluster
            new_centroids[c] = Z[seed_rng.integers(n_genes)]

    # --- Convergence check: stop when assignments stop changing. ---
    if np.array_equal(new_labels, labels) and iteration > 0:
        labels = new_labels
        centroids = new_centroids
        print(f"K-means converged after {iteration} iterations")
        break
    labels, centroids = new_labels, new_centroids

# Report cluster sizes.
for c in range(k):
    print(f"Cluster {c}: {(labels == c).sum()} genes")

# Compute total within-cluster sum of squares (clustering objective).
wcss = sum(((Z[labels == c] - centroids[c]) ** 2).sum() for c in range(k))
print(f"Within-cluster sum of squares (WCSS): {wcss:.4f}")

# ----------------------------------------------------------------------
# 3. Order genes so that cluster members are contiguous (cluster-ordered).
#    Order clusters by the sample at which their centroid peaks, so the
#    heatmap reads as a progression left-to-right and top-to-bottom.
# ----------------------------------------------------------------------
cluster_peak = centroids.argmax(axis=1)             # peak sample per cluster
cluster_order = np.argsort(cluster_peak)            # clusters along progression
row_order = np.concatenate([np.where(labels == c)[0] for c in cluster_order])
Z_ordered = Z[row_order]

# ----------------------------------------------------------------------
# 4. Cluster-ordered expression heatmap, color scale -3 to 3.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(Z_ordered, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)

# Draw white lines at cluster boundaries for readability.
boundary = 0
for c in cluster_order:
    boundary += (labels == c).sum()
    if boundary < n_genes:
        ax.axhline(boundary - 0.5, color="white", linewidth=1.2)

ax.set_xlabel("Sample (ordered condition progression)")
ax.set_ylabel("Genes (grouped by k-means cluster)")
ax.set_title(f"Cluster-ordered gene-expression heatmap (k={k}, seed=1)")
fig.colorbar(im, ax=ax, label="z-score")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.3.1_s2.png")

# ----------------------------------------------------------------------
# 5. Separate check: confirm each cluster is a co-expression module whose
#    genes rise and fall together across conditions. We measure, per
#    cluster, the average correlation of each gene's profile with its
#    cluster centroid (module coherence). High coherence + distinct
#    centroid peak samples => the PCA progression resolves into groups.
# ----------------------------------------------------------------------
print("\n--- Co-expression module check ---")
overall_coherence = []
for c in cluster_order:
    members = Z[labels == c]
    cen = centroids[c]
    # Correlation of each member gene with the cluster's mean profile.
    corrs = [np.corrcoef(g, cen)[0, 1] for g in members]
    coherence = np.mean(corrs)
    overall_coherence.extend(corrs)
    print(f"Cluster {c}: mean gene-vs-centroid correlation = {coherence:.3f}, "
          f"centroid peaks at sample {cen.argmax()}")

print(f"Overall mean within-cluster coherence: {np.mean(overall_coherence):.3f}")
print(f"Distinct centroid peak samples across clusters: "
      f"{sorted(int(centroids[c].argmax()) for c in cluster_order)}")

# One-sentence explanation of why this check confirms the result:
print("\nWhy the check confirms the result: high within-cluster correlation "
      "with the centroid means each cluster's genes move up and down in lockstep "
      "across the conditions, and the distinct centroid peak samples show these "
      "modules occupy different points along the PCA progression -- so the "
      "clustering has resolved that continuous progression into concrete "
      "co-expression gene groups.")
