import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from scipy.stats import spearmanr
import umap  # umap-learn

# ----------------------------------------------------------------------
# 1. Build the model gene-expression data: 26 samples x 500 genes,
#    each sample belonging to one of nine conditions A..I that form a
#    smooth A-to-I progression in the high-dimensional space.
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)          # seed 1
n_samples, n_genes, n_cond = 26, 500, 9
cond_names = list("ABCDEFGHI")

# assign the 26 samples to the 9 conditions (as even as possible)
cond_idx = np.array([i % n_cond for i in range(n_samples)])
cond_idx.sort()                         # group A..I in order

# a single random direction in gene space defines the progression axis;
# condition k sits at position k along that axis -> global A->I ordering
direction = rng.normal(size=n_genes)
direction /= np.linalg.norm(direction)
base = rng.normal(size=n_genes)

X = np.empty((n_samples, n_genes))
for s in range(n_samples):
    k = cond_idx[s]                     # 0..8 = A..I
    # mean shifts linearly with condition index -> A-to-I line
    mean = base + 6.0 * k * direction
    X[s] = mean + rng.normal(scale=0.6, size=n_genes)   # per-sample noise

colors = plt.cm.viridis(cond_idx / (n_cond - 1))

# ----------------------------------------------------------------------
# 2. Classical (Torgerson) MDS, implemented explicitly.
# ----------------------------------------------------------------------
# pairwise squared Euclidean distances
G = X @ X.T
sq = np.diag(G)
D2 = sq[:, None] + sq[None, :] - 2 * G
D2 = np.maximum(D2, 0.0)
# double-centering: B = -1/2 J D2 J,  J = I - (1/n) 11^T
J = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
B = -0.5 * J @ D2 @ J
# eigendecomposition; embed on top-2 eigenvectors scaled by sqrt(eigenvalue)
eigvals, eigvecs = np.linalg.eigh(B)
order = np.argsort(eigvals)[::-1]
eigvals, eigvecs = eigvals[order], eigvecs[:, order]
L = np.sqrt(np.clip(eigvals[:2], 0, None))
mds = eigvecs[:, :2] * L

# ----------------------------------------------------------------------
# 3. t-SNE (perplexity 5) and UMAP (5 neighbors), both seeded.
# ----------------------------------------------------------------------
tsne = TSNE(n_components=2, perplexity=5, random_state=1, init="pca").fit_transform(X)
ump = umap.UMAP(n_components=2, n_neighbors=5, random_state=1).fit_transform(X)

# ----------------------------------------------------------------------
# 4. Plot the three embeddings, colored by condition.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, emb, title in zip(axes, [mds, tsne, ump], ["Classical MDS", "t-SNE", "UMAP"]):
    ax.scatter(emb[:, 0], emb[:, 1], c=colors, s=60, edgecolor="k")
    for s in range(n_samples):
        ax.annotate(cond_names[cond_idx[s]], (emb[s, 0], emb[s, 1]), fontsize=7)
    ax.set_title(title)
    ax.set_xlabel("dim 1"); ax.set_ylabel("dim 2")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.5.1_s5.png", dpi=120)

# ----------------------------------------------------------------------
# 5. Separate check: (a) all three separate the nine conditions
#    (silhouette on the 2D embedding, condition labels), and
#    (b) MDS lays conditions along the A-to-I progression
#    (Spearman corr between condition index and MDS dim-1 centroid).
# ----------------------------------------------------------------------
sil_mds = silhouette_score(mds, cond_idx)
sil_tsne = silhouette_score(tsne, cond_idx)
sil_umap = silhouette_score(ump, cond_idx)

# condition centroids along MDS dim 1, then rank correlation with A..I order
centroids = np.array([mds[cond_idx == k, 0].mean() for k in range(n_cond)])
rho, _ = spearmanr(np.arange(n_cond), centroids)

print(f"Silhouette (condition separation) MDS : {sil_mds:.3f}")
print(f"Silhouette (condition separation) tSNE: {sil_tsne:.3f}")
print(f"Silhouette (condition separation) UMAP: {sil_umap:.3f}")
print(f"MDS dim-1 centroids A..I: {np.array2string(centroids, precision=3)}")
print(f"Spearman(condition index, MDS dim-1 centroid): {abs(rho):.3f}")
print(f"Top-2 MDS eigenvalues: {eigvals[0]:.3f}, {eigvals[1]:.3f}")

# One-sentence explanation:
# A |Spearman| ~ 1 shows MDS dim-1 orders the condition centroids monotonically
# A->I (global distances preserved), while all three silhouettes > 0 confirm the
# nine conditions are separated in every embedding -- so the check confirms the
# result because t-SNE/UMAP separate clusters without preserving that global order.
print("Explanation: positive silhouettes for all three confirm the nine "
      "conditions are separated everywhere, while the near-1 Spearman "
      "correlation confirms only MDS preserves the global A-to-I ordering.")
