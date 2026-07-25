import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 0. Build the model data: 26 samples x 500 genes, 9 conditions A..I.
#    Conditions sit on a 1D latent "progression" (A=0 ... I=8); gene
#    means slide along a few loading directions with that progression,
#    plus per-sample noise. This gives a real A-to-I global ordering
#    plus tight per-condition local clusters.
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)          # seed 1
n_genes = 500
cond_names = list("ABCDEFGHI")          # nine conditions
sizes = [3, 3, 3, 3, 3, 3, 3, 3, 2]     # replicate counts, sum = 26
labels = np.array([c for c, s in zip(cond_names, sizes) for _ in range(s)])
cond_idx = np.array([cond_names.index(c) for c in labels])   # 0..8 progression
n = len(labels)                         # 26 samples

# a handful of smooth loading directions that vary with the progression
n_dir = 5
loadings = rng.normal(size=(n_dir, n_genes))
t = cond_idx.astype(float)              # progression coordinate per sample
# nonlinear progression signal along each direction
signal = np.stack([t, t**1.5, np.sin(0.6 * t), np.cos(0.4 * t), t * 0.5], axis=1)
X = signal @ loadings                   # (26, 500) smooth mean part
X += rng.normal(scale=0.6, size=(n, n_genes))   # per-sample noise -> tight clusters
X = X - X.mean(axis=0, keepdims=True)   # center genes

# distinct color per condition
palette = plt.cm.viridis(np.linspace(0, 1, len(cond_names)))
colors = np.array([palette[i] for i in cond_idx])


def pairwise_sq_dists(Y):
    """Squared Euclidean distance matrix."""
    s = np.sum(Y * Y, axis=1)
    D2 = s[:, None] + s[None, :] - 2.0 * (Y @ Y.T)
    return np.maximum(D2, 0.0)


# ----------------------------------------------------------------------
# 1. Classical (Torgerson) MDS -- explicit double-centering + eigendecomp.
#    Preserves GLOBAL distances, so it recovers the A-to-I ordering.
# ----------------------------------------------------------------------
D2 = pairwise_sq_dists(X)                        # squared distances
J = np.eye(n) - np.ones((n, n)) / n              # centering operator
B = -0.5 * J @ D2 @ J                            # double-centered Gram matrix
evals, evecs = np.linalg.eigh(B)                 # symmetric eigendecomposition
order = np.argsort(evals)[::-1]                  # largest eigenvalues first
evals, evecs = evals[order], evecs[:, order]
L = np.sqrt(np.maximum(evals[:2], 0.0))          # top-2 positive eigenvalues
mds = evecs[:, :2] * L                           # coordinates = V * sqrt(lambda)


# ----------------------------------------------------------------------
# 2. t-SNE -- explicit: perplexity-calibrated P, Student-t Q, KL descent.
# ----------------------------------------------------------------------
def hi_dim_affinities(X, perplexity):
    """Row-normalized Gaussian affinities P_{j|i} matched to a perplexity."""
    D2 = pairwise_sq_dists(X)
    P = np.zeros((n, n))
    target = np.log(perplexity)                  # target entropy (nats)
    for i in range(n):
        d = D2[i].copy()
        d[i] = np.inf                            # exclude self
        beta_lo, beta_hi, beta = -np.inf, np.inf, 1.0   # beta = 1/(2 sigma^2)
        for _ in range(50):                      # binary search on beta
            w = np.exp(-d * beta)
            w[i] = 0.0
            sw = w.sum()
            if sw == 0:
                sw = 1e-12
            p = w / sw
            H = -np.sum(p[p > 0] * np.log(p[p > 0]))     # Shannon entropy
            if H < target:                       # too little spread -> lower beta
                beta_hi = beta
                beta = (beta + beta_lo) / 2 if beta_lo > -np.inf else beta / 2
            else:                                # too much spread -> raise beta
                beta_lo = beta
                beta = (beta + beta_hi) / 2 if beta_hi < np.inf else beta * 2
        P[i] = p
    P = (P + P.T) / (2 * n)                       # symmetrize + normalize
    return np.maximum(P, 1e-12)


P = hi_dim_affinities(X, perplexity=5)
Y = 1e-4 * rng.normal(size=(n, 2))               # small random init
Yprev = Y.copy()
lr, n_iter = 200.0, 1000
for it in range(n_iter):
    momentum = 0.5 if it < 250 else 0.8
    Pcur = P * 12.0 if it < 250 else P           # early exaggeration
    num = 1.0 / (1.0 + pairwise_sq_dists(Y))     # Student-t kernel
    np.fill_diagonal(num, 0.0)
    Q = np.maximum(num / num.sum(), 1e-12)       # low-dim affinities
    # gradient of KL divergence
    PQ = (Pcur - Q) * num
    grad = 4.0 * (np.diag(PQ.sum(axis=1)) - PQ) @ Y
    Ynew = Y - lr * grad + momentum * (Y - Yprev)
    Yprev, Y = Y, Ynew
    Y -= Y.mean(axis=0)                          # recenter
tsne = Y


# ----------------------------------------------------------------------
# 3. UMAP -- explicit: fuzzy simplicial set + attractive/repulsive SGD.
# ----------------------------------------------------------------------
k = 5
D = np.sqrt(pairwise_sq_dists(X))                # Euclidean distances
knn = np.argsort(D, axis=1)[:, 1:k + 1]          # k nearest neighbors (no self)
rho = np.array([D[i, knn[i]].min() for i in range(n)])   # dist to nearest nb
target = np.log2(k)                              # desired fuzzy set size
sigmas = np.ones(n)
for i in range(n):                               # binary search per-point sigma
    lo, hi, sig = 0.0, np.inf, 1.0
    for _ in range(50):
        val = np.sum(np.exp(-(np.maximum(D[i, knn[i]] - rho[i], 0.0)) / sig))
        if val > target:
            hi = sig; sig = (lo + hi) / 2
        else:
            lo = sig; sig = (lo + hi) / 2 if hi < np.inf else sig * 2
    sigmas[i] = sig
A = np.zeros((n, n))                             # directed membership strengths
for i in range(n):
    A[i, knn[i]] = np.exp(-(np.maximum(D[i, knn[i]] - rho[i], 0.0)) / sigmas[i])
G = A + A.T - A * A.T                            # probabilistic t-conorm (symmetric)

a, b = 1.929, 0.7915                             # low-dim curve params (min_dist~0.1)
Y = rng.normal(scale=1.0, size=(n, 2))           # random init
edges = [(i, j, G[i, j]) for i in range(n) for j in range(n) if i != j and G[i, j] > 1e-3]
n_epochs, alpha0 = 500, 1.0
for ep in range(n_epochs):
    alpha = alpha0 * (1 - ep / n_epochs)         # decaying learning rate
    for i, j, w in edges:
        if rng.random() > w:                     # sample edge by its strength
            continue
        d = Y[i] - Y[j]
        d2 = np.dot(d, d) + 1e-6
        # attractive gradient between connected points
        coef = (-2 * a * b * d2 ** (b - 1)) / (1 + a * d2 ** b)
        Y[i] += alpha * np.clip(coef * d, -4, 4)
        Y[j] -= alpha * np.clip(coef * d, -4, 4)
        for _ in range(5):                       # negative sampling (repulsion)
            kk = rng.integers(n)
            if kk == i:
                continue
            dn = Y[i] - Y[kk]
            dn2 = np.dot(dn, dn) + 1e-6
            coefn = (2 * b) / ((1e-3 + dn2) * (1 + a * dn2 ** b))
            Y[i] += alpha * np.clip(coefn * dn, -4, 4)
umap_emb = Y - Y.mean(axis=0)


# ----------------------------------------------------------------------
# 4. Plot the three embeddings, colored by condition.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, emb, title in zip(axes, [mds, tsne, umap_emb],
                          ["Classical MDS", "t-SNE (perp=5)", "UMAP (k=5)"]):
    ax.scatter(emb[:, 0], emb[:, 1], c=colors, s=60, edgecolor="k", linewidth=0.4)
    for i, name in enumerate(cond_names):        # label each condition centroid
        c = emb[cond_idx == i].mean(axis=0)
        ax.annotate(name, c, fontsize=11, fontweight="bold")
    ax.set_title(title)
    ax.set_xlabel("dim 1"); ax.set_ylabel("dim 2")
handles = [plt.Line2D([0], [0], marker='o', ls='', markerfacecolor=palette[i],
                      markeredgecolor='k', label=cond_names[i]) for i in range(9)]
fig.legend(handles=handles, loc="upper center", ncol=9, title="condition")
fig.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.5.1_s4.png", dpi=130)


# ----------------------------------------------------------------------
# 5. Separate check: do all three separate the nine conditions, and does
#    MDS follow the A-to-I progression while t-SNE/UMAP form tight local
#    clusters (whose between-cluster distances carry no global meaning)?
# ----------------------------------------------------------------------
def neighbor_purity(emb, kk=3):
    """Fraction of each point's kk nearest neighbors sharing its condition."""
    Dd = np.sqrt(pairwise_sq_dists(emb))
    np.fill_diagonal(Dd, np.inf)
    nn = np.argsort(Dd, axis=1)[:, :kk]
    return np.mean([np.mean(cond_idx[nn[i]] == cond_idx[i]) for i in range(n)])


def spearman(x, y):
    """Spearman rank correlation via Pearson on ranks."""
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return np.corrcoef(rx, ry)[0, 1]


print("=== Separation: local neighbor purity (1.0 = every neighbor same condition) ===")
for emb, name in zip([mds, tsne, umap_emb], ["MDS ", "tSNE", "UMAP"]):
    print(f"{name} neighbor purity (k=3): {neighbor_purity(emb):.3f}")

# MDS global structure: does an MDS axis track the A-to-I condition order?
rho1 = spearman(mds[:, 0], cond_idx)
rho2 = spearman(mds[:, 1], cond_idx)
best_axis = 1 if abs(rho2) > abs(rho1) else 0
print("\n=== MDS global progression (Spearman corr of MDS axis vs A..I index) ===")
print(f"MDS dim1 vs condition order: {rho1:.3f}")
print(f"MDS dim2 vs condition order: {rho2:.3f}")
print(f"Best MDS axis |Spearman| vs A-to-I: {max(abs(rho1), abs(rho2)):.3f} (axis {best_axis + 1})")

# t-SNE / UMAP: between-cluster distances are NOT meaningful -> weak global order
print("\n=== t-SNE / UMAP global order (should be weaker than MDS) ===")
for emb, name in zip([tsne, umap_emb], ["tSNE", "UMAP"]):
    r1, r2 = spearman(emb[:, 0], cond_idx), spearman(emb[:, 1], cond_idx)
    print(f"{name} best axis |Spearman| vs A-to-I: {max(abs(r1), abs(r2)):.3f}")

print("\n=== Interpretation ===")
print("Explanation: high neighbor purity in all three embeddings confirms every")
print("method separates the nine conditions, while a near-1 MDS-axis correlation")
print("with the A-to-I index (vs. a lower t-SNE/UMAP correlation) confirms that")
print("only MDS preserves the global progression -- t-SNE and UMAP keep tight local")
print("clusters whose between-cluster distances are not meaningful.")
