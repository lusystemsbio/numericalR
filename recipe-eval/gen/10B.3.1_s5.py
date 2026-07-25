import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 0. Build a synthetic gene-expression matrix (500 genes x 26 samples).
#    We embed several latent "programs" that rise/fall along an ordered
#    progression of the 26 conditions -- this mimics the PCA progression
#    of Part 10A, which we will resolve into gene groups below.
# ---------------------------------------------------------------
rng = np.random.default_rng(1)          # seed 1 (data generation)
n_genes, n_samples = 500, 26
t = np.linspace(0, 1, n_samples)        # ordered "progression" axis

# latent temporal programs across the progression (co-expression modules)
programs = np.vstack([
    np.sin(2 * np.pi * t),              # early up / late down
    np.cos(2 * np.pi * t),              # up at ends, down in middle
    t,                                  # monotone rising
    1 - t,                              # monotone falling
    np.exp(-((t - 0.5) ** 2) / 0.02),   # transient burst in the middle
])
n_prog = programs.shape[0]

# assign each gene mostly to one program, add small mixing + noise
X = np.zeros((n_genes, n_samples))
gene_program = rng.integers(0, n_prog, size=n_genes)
for g in range(n_genes):
    w = np.zeros(n_prog)
    w[gene_program[g]] = 1.0
    w += rng.normal(0, 0.15, n_prog)                 # slight cross-talk
    X[g] = w @ programs + rng.normal(0, 0.4, n_samples)

# ---------------------------------------------------------------
# 1. z-score each gene (row) -> mean 0, sd 1 across its 26 samples
# ---------------------------------------------------------------
Z = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
print("Data shape (genes x samples):", Z.shape)
print("Per-gene mean after z-score (should be ~0):", np.round(Z.mean(axis=1).mean(), 6))
print("Per-gene std  after z-score (should be ~1):", np.round(Z.std(axis=1).mean(), 6))

# ---------------------------------------------------------------
# 2. k-means clustering of the z-scored genes, implemented explicitly
# ---------------------------------------------------------------
k = 5
km_rng = np.random.default_rng(1)       # seed 1 (clustering)

# initialize centroids by picking k distinct genes at random
init_idx = km_rng.choice(n_genes, size=k, replace=False)
centroids = Z[init_idx].copy()

labels = np.zeros(n_genes, dtype=int)
for iteration in range(100):
    # -- assignment step: each gene -> nearest centroid (Euclidean) --
    # squared distances from every gene to every centroid
    d2 = ((Z[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    new_labels = d2.argmin(axis=1)

    # -- update step: recompute each centroid as mean of its members --
    new_centroids = np.zeros_like(centroids)
    for c in range(k):
        members = Z[new_labels == c]
        if len(members) == 0:
            # empty cluster: reseed to a random gene
            new_centroids[c] = Z[km_rng.integers(n_genes)]
        else:
            new_centroids[c] = members.mean(axis=0)

    # -- convergence check --
    if np.array_equal(new_labels, labels) and np.allclose(new_centroids, centroids):
        labels = new_labels
        centroids = new_centroids
        print("k-means converged at iteration:", iteration)
        break
    labels, centroids = new_labels, new_centroids

# within-cluster sum of squares (inertia)
inertia = sum(((Z[labels == c] - centroids[c]) ** 2).sum() for c in range(k))
print("Number of clusters k:", k)
for c in range(k):
    print(f"Cluster {c} size (genes):", int((labels == c).sum()))
print("Total within-cluster sum of squares (inertia):", round(float(inertia), 4))

# ---------------------------------------------------------------
# 3. order genes by cluster so co-expressed genes are contiguous
# ---------------------------------------------------------------
order = np.argsort(labels, kind="stable")
Z_ordered = Z[order]
labels_ordered = labels[order]

# ---------------------------------------------------------------
# 4. cluster-ordered expression heatmap (color scale -3..3)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(Z_ordered, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_xlabel("Samples (ordered progression, 26 conditions)")
ax.set_ylabel("Genes (ordered by k-means cluster)")
ax.set_title(f"Cluster-ordered gene-expression heatmap (k={k})")

# draw white lines at the boundaries between clusters
boundaries = np.where(np.diff(labels_ordered) != 0)[0] + 0.5
for b in boundaries:
    ax.axhline(b, color="black", linewidth=0.8)
fig.colorbar(im, ax=ax, label="z-score")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.3.1_s5.png", dpi=120)

# ---------------------------------------------------------------
# 5. CHECK: do clusters behave as co-expression modules that rise/fall
#    together across conditions? Summarize each cluster's mean profile
#    and its internal coherence (mean pairwise correlation of members).
# ---------------------------------------------------------------
print("\n--- Co-expression module check ---")
for c in range(k):
    members = Z[labels == c]
    profile = members.mean(axis=0)
    # peak / trough condition of the cluster's mean profile
    print(f"Cluster {c}: peak at sample {int(profile.argmax())}, "
          f"trough at sample {int(profile.argmin())}, "
          f"profile range = {round(float(profile.ptp()), 3)}")
    # internal coherence: mean off-diagonal correlation among member genes
    if len(members) > 1:
        corr = np.corrcoef(members)
        off_diag = corr[np.triu_indices_from(corr, k=1)]
        print(f"Cluster {c}: mean within-cluster gene-gene correlation = "
              f"{round(float(off_diag.mean()), 3)}")

# average within-cluster correlation vs. overall correlation
all_corr = np.corrcoef(Z)
overall = all_corr[np.triu_indices_from(all_corr, k=1)].mean()
within_vals = []
for c in range(k):
    m = Z[labels == c]
    if len(m) > 1:
        cc = np.corrcoef(m)
        within_vals.append(cc[np.triu_indices_from(cc, k=1)].mean())
print("Overall mean gene-gene correlation:", round(float(overall), 3))
print("Mean within-cluster gene-gene correlation:", round(float(np.mean(within_vals)), 3))

# Explanation (one sentence):
print("\nWhy this confirms the result: because within-cluster genes are far "
      "more correlated than average and each cluster's mean profile has a "
      "distinct peak/trough along the ordered conditions, the heatmap's blocks "
      "are genuine co-expression modules that rise and fall together, resolving "
      "the smooth PCA progression of Part 10A into specific gene groups.")
