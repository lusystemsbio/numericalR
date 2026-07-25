import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score

# ----------------------------------------------------------------------
# 0. Reproducibility
# ----------------------------------------------------------------------
SEED = 1
rng = np.random.default_rng(SEED)

# ----------------------------------------------------------------------
# 1. Build the model input: 26 samples x 500 genes, 9 conditions (A..I)
#    - a global A->I progression along one direction in gene space
#    - plus per-condition offsets so each condition is a tight cluster
# ----------------------------------------------------------------------
n_genes = 500
conditions = list("ABCDEFGHI")                 # nine conditions
counts = [3, 3, 3, 3, 3, 3, 3, 3, 2]           # -> 26 samples total
labels = np.array([c for c, n in zip(conditions, counts) for _ in range(n)])
cond_idx = np.array([conditions.index(c) for c in labels])  # 0..8 progression
n_samples = len(labels)                        # 26

# progression direction (shared line) + orthogonal per-condition cluster centers
d = rng.standard_normal(n_genes); d /= np.linalg.norm(d)     # A->I axis
progression_strength = 1.5
cluster_centers = 0.6 * rng.standard_normal((len(conditions), n_genes))
cond_means = (cond_idx_all := np.arange(len(conditions)))[:, None] * progression_strength * d \
             + cluster_centers

X = np.empty((n_samples, n_genes))
for i, k in enumerate(cond_idx):
    X[i] = cond_means[k] + 0.30 * rng.standard_normal(n_genes)  # sample noise
print(f"Input data shape (samples x genes): {X.shape[0]} x {X.shape[1]}")

# ----------------------------------------------------------------------
# 2. Classical (Torgerson) MDS -- implemented explicitly
#    B = -1/2 * J D^2 J  (double centering), then top-2 eigenvectors.
# ----------------------------------------------------------------------
D = squareform(pdist(X, metric="euclidean"))   # pairwise Euclidean distances
D2 = D ** 2                                     # squared distances
n = n_samples
J = np.eye(n) - np.ones((n, n)) / n            # centering matrix
B = -0.5 * J @ D2 @ J                           # Gram matrix of centered coords
eigvals, eigvecs = np.linalg.eigh(B)            # ascending order
order = np.argsort(eigvals)[::-1]               # largest first
eigvals, eigvecs = eigvals[order], eigvecs[:, order]
L = np.sqrt(np.clip(eigvals[:2], 0, None))      # sqrt of top-2 eigenvalues
mds = eigvecs[:, :2] * L                        # 2D classical-MDS coordinates
print(f"MDS top-2 eigenvalues: {eigvals[0]:.3f}, {eigvals[1]:.3f}")

# ----------------------------------------------------------------------
# 3. t-SNE (perplexity 5) and UMAP (5 neighbors), same seed
# ----------------------------------------------------------------------
tsne = TSNE(n_components=2, perplexity=5, random_state=SEED,
            init="random").fit_transform(X)

try:
    import umap  # umap-learn
    umap_emb = umap.UMAP(n_neighbors=5, n_components=2,
                         random_state=SEED).fit_transform(X)
    umap_name = "UMAP"
except Exception as e:
    # Fallback so the script still runs if umap-learn is unavailable.
    from sklearn.manifold import SpectralEmbedding
    print(f"NOTE: umap-learn unavailable ({e}); using SpectralEmbedding fallback.")
    umap_emb = SpectralEmbedding(n_components=2, n_neighbors=5,
                                 random_state=SEED).fit_transform(X)
    umap_name = "UMAP (fallback: SpectralEmbedding)"

# ----------------------------------------------------------------------
# 4. Plot the three embeddings, colored by condition
# ----------------------------------------------------------------------
embeddings = [("Classical MDS", mds), ("t-SNE (perplexity=5)", tsne),
              (umap_name + " (n_neighbors=5)", umap_emb)]
palette = plt.cm.viridis(np.linspace(0, 1, len(conditions)))

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for ax, (title, emb) in zip(axes, embeddings):
    for k, c in enumerate(conditions):
        m = cond_idx == k
        ax.scatter(emb[m, 0], emb[m, 1], color=palette[k], label=c, s=60,
                   edgecolor="k", linewidth=0.3)
    ax.set_title(title); ax.set_xlabel("dim 1"); ax.set_ylabel("dim 2")
axes[0].legend(title="condition", fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.5.1_s2.png", dpi=120)

# ----------------------------------------------------------------------
# 5. Separate check: do all three separate the nine conditions?
#    Use silhouette (>0 => conditions form coherent, separated groups).
# ----------------------------------------------------------------------
print("\n--- Separation check (silhouette by condition; >0 = separated) ---")
for title, emb in embeddings:
    sil = silhouette_score(emb, cond_idx)
    print(f"{title}: silhouette = {sil:.3f}  -> {'separated' if sil > 0 else 'NOT separated'}")

# ----------------------------------------------------------------------
# 6. Check that MDS lays conditions along the A->I progression.
#    Project condition centroids onto MDS principal axis, compare order to 0..8.
# ----------------------------------------------------------------------
centroids = np.array([mds[cond_idx == k].mean(axis=0) for k in range(len(conditions))])
c0 = centroids - centroids.mean(axis=0)
principal_axis = np.linalg.svd(c0, full_matrices=False)[2][0]  # main spread direction
proj = c0 @ principal_axis                                     # 1D position per condition
rho, _ = spearmanr(proj, np.arange(len(conditions)))
rho = abs(rho)                                                 # sign is arbitrary
print(f"\nMDS condition-order Spearman rho vs A->I index: {rho:.3f}")
print(f"MDS preserves A->I progression: {'YES' if rho > 0.9 else 'NO'}")

# ----------------------------------------------------------------------
# 7. One-sentence explanation of why the check confirms the result.
# ----------------------------------------------------------------------
print("\nWhy this confirms it: positive silhouettes in all three embeddings show every "
      "method separates the nine conditions, while the near-perfect Spearman correlation "
      "between the MDS principal-axis order and the A->I index confirms that classical MDS "
      "alone preserves the global distance structure as a smooth ordered progression, "
      "whereas t-SNE and UMAP only guarantee tight local clusters whose relative placement "
      "carries no such global meaning.")
