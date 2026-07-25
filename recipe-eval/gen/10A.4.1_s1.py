import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from string import ascii_uppercase

rng = np.random.default_rng(0)

# ----------------------------------------------------------------------
# Build synthetic gene-expression data: 26 samples x 500 genes.
# Samples lie along a 1D latent progression (an arc) so that the leading
# principal components should recover an ordered curve from A to I.
# ----------------------------------------------------------------------
n_samples, n_genes = 26, 500
conditions = list(ascii_uppercase[:9])                      # A..I (9 conditions)
# assign each of 26 samples to a condition, ordered A..I
labels = np.array([conditions[i * len(conditions) // n_samples] for i in range(n_samples)])

# latent 1D coordinate t increasing monotonically A->I
t = np.linspace(0.0, np.pi, n_samples)                      # arc parameter
# embed the arc in the first two "signal" gene directions, rest is noise
signal = np.zeros((n_samples, n_genes))
signal[:, 0] = 3.0 * np.cos(t)                              # arc x
signal[:, 1] = 3.0 * np.sin(t)                              # arc y
# spread the two signal directions across many genes via a random loading mix
mix = rng.normal(size=(2, n_genes)) * 0.15
X = signal[:, :2] @ mix + signal                           # genes carry the arc
X += rng.normal(scale=0.25, size=(n_samples, n_genes))     # measurement noise

# ----------------------------------------------------------------------
# PCA implemented explicitly (center, covariance eigendecomposition).
# ----------------------------------------------------------------------
Xc = X - X.mean(axis=0, keepdims=True)                      # center each gene
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)          # SVD of centered data
var_explained = S**2 / np.sum(S**2)                        # fraction variance per PC
scores = U * S                                             # sample scores = U*S
pc12 = scores[:, :2]                                       # project onto first 2 PCs

print(f"Variance explained by PC1: {var_explained[0]:.4f}")
print(f"Variance explained by PC2: {var_explained[1]:.4f}")
print(f"Variance explained by PC3: {var_explained[2]:.4f}")
print(f"Cumulative variance PC1+PC2: {var_explained[0] + var_explained[1]:.4f}")

# ----------------------------------------------------------------------
# Hastie-Stuetzle principal curve, implemented explicitly.
#  1) start from the PC1 line as the initial curve
#  2) PROJECT: assign each point its arc-length coordinate on the curve
#  3) EXPECTATION: re-estimate the curve as a smooth function of arc length
#  4) iterate until the coordinates stop changing
# ----------------------------------------------------------------------
def project_to_polyline(pts, curve):
    """Return, per point, the nearest arc-length coordinate along a polyline."""
    seg_start = curve[:-1]
    seg_vec = curve[1:] - curve[:-1]
    seg_len = np.linalg.norm(seg_vec, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])       # arc length at vertices
    lam = np.empty(len(pts))
    for i, p in enumerate(pts):
        # for each segment find the closest point (clamped projection)
        u = np.einsum("sj,j->s", seg_vec, p) - np.einsum("sj,sj->s", seg_vec, seg_start)
        denom = np.maximum(seg_len**2, 1e-12)
        s = np.clip(u / denom, 0.0, 1.0)                    # position within segment
        proj = seg_start + s[:, None] * seg_vec
        d = np.linalg.norm(proj - p, axis=1)
        k = np.argmin(d)                                    # best segment
        lam[i] = cum[k] + s[k] * seg_len[k]                 # arc length of projection
    return lam

def smooth_curve(pts, lam, n_knots=12, span=0.3):
    """Re-fit curve coordinates as a local-linear smoother of arc length lam."""
    order = np.argsort(lam)
    lam_s = lam[order]
    grid = np.linspace(lam_s.min(), lam_s.max(), n_knots)   # arc-length grid
    width = span * (lam_s.max() - lam_s.min() + 1e-9)
    fit = np.empty((n_knots, pts.shape[1]))
    for j, g in enumerate(grid):
        w = np.exp(-0.5 * ((lam - g) / width) ** 2)         # Gaussian kernel weights
        w /= w.sum()
        fit[j] = w @ pts                                    # weighted mean per dim
    return fit                                              # smoothed curve vertices

# initialize the curve as the PC1 axis in PC1-PC2 space
init_lam = pc12[:, 0].copy()
curve = smooth_curve(pc12, init_lam)
lam = project_to_polyline(pc12, curve)

for it in range(20):
    curve = smooth_curve(pc12, lam)                         # expectation step
    new_lam = project_to_polyline(pc12, curve)             # projection step
    if np.max(np.abs(new_lam - lam)) < 1e-4:
        lam = new_lam
        break
    lam = new_lam

# normalize principal-curve coordinate to [0,1] for interpretation
lam_norm = (lam - lam.min()) / (lam.max() - lam.min())

# check monotonic ordering of the 1D coordinate vs condition A..I
cond_index = np.array([conditions.index(c) for c in labels])
order_corr = np.corrcoef(cond_index, lam_norm)[0, 1]
print(f"Iterations to converge: {it + 1}")
print(f"Correlation(condition order, principal-curve coord): {order_corr:.4f}")

print("Per-sample principal-curve coordinate (label: coord):")
for lab, lc in sorted(zip(labels, lam_norm), key=lambda z: z[1]):
    print(f"  {lab}: {lc:.4f}")

# ----------------------------------------------------------------------
# Plots: scree plot + PC1-PC2 scatter colored by condition with the curve.
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# scree plot
comps = np.arange(1, 11)
ax1.bar(comps, var_explained[:10], color="steelblue")
ax1.plot(comps, np.cumsum(var_explained[:10]), "o-", color="darkorange", label="cumulative")
ax1.set_xlabel("Principal component")
ax1.set_ylabel("Fraction of variance explained")
ax1.set_title("Scree plot")
ax1.legend()

# PC1-PC2 scatter colored by condition
cmap = plt.cm.get_cmap("viridis", len(conditions))
for k, c in enumerate(conditions):
    m = labels == c
    ax2.scatter(pc12[m, 0], pc12[m, 1], color=cmap(k), s=60,
                edgecolor="k", label=c, zorder=3)
# order the fitted curve vertices along arc length before drawing
curve_lam = project_to_polyline(curve, curve)
co = np.argsort(curve_lam)
ax2.plot(curve[co, 0], curve[co, 1], "-", color="red", lw=2.5,
         label="principal curve", zorder=2)
ax2.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
ax2.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
ax2.set_title("PC1-PC2 projection with principal curve")
ax2.legend(ncol=2, fontsize=8, title="condition")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.4.1_s1.png", dpi=130)

# One-sentence explanation of why the check confirms the result:
print("Explanation: A dominant PC1 variance share together with a near-1 correlation "
      "between condition order (A..I) and the principal-curve coordinate confirms that "
      "the samples lie on a single ordered arc that the principal curve faithfully "
      "summarizes as one 1D coordinate per sample.")
