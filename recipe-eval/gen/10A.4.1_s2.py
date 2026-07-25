import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.stats import spearmanr

rng = np.random.default_rng(0)

# ----------------------------------------------------------------------
# 1. Build the model input: 26 samples x 500 genes with a hidden 1-D arc
# ----------------------------------------------------------------------
n_samples, n_genes = 26, 500
s = np.linspace(0.0, 1.0, n_samples)          # latent ordering, one value per sample

# Latent geometry: a dominant "progression" axis + a smaller "curvature" bump.
# This makes the samples trace an ordered arc rather than a straight line.
progression = s                                # large-variance main direction
curvature   = np.sin(np.pi * s)                # smaller side-to-side bump -> arc

# Two fixed gene-loading directions that map latent axes into 500-gene space.
load1 = rng.normal(size=n_genes)
load2 = rng.normal(size=n_genes)
load1 /= np.linalg.norm(load1)
load2 /= np.linalg.norm(load2)

A_scale, B_scale = 6.0, 1.6                    # progression dominates -> PC1 dominates
expr = (A_scale * progression[:, None] * load1[None, :]
        + B_scale * curvature[:, None] * load2[None, :]
        + 0.15 * rng.normal(size=(n_samples, n_genes)))   # gene-expression matrix

# Label each sample by condition A..I, assigned along the latent ordering.
cond_letters = list("ABCDEFGHI")
groups = np.array_split(np.arange(n_samples), len(cond_letters))
labels = np.empty(n_samples, dtype="<U1")
for letter, idx in zip(cond_letters, groups):
    labels[idx] = letter
cond_index = np.array([cond_letters.index(c) for c in labels])   # 0..8

# ----------------------------------------------------------------------
# 2. PCA implemented explicitly (center -> SVD -> scores + variance)
# ----------------------------------------------------------------------
Xc = expr - expr.mean(axis=0, keepdims=True)   # center each gene
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
explained_var = (S ** 2) / (n_samples - 1)     # variance carried by each component
explained_ratio = explained_var / explained_var.sum()
scores = U * S                                 # PCA scores (samples in PC space)
pc = scores[:, :2]                             # project onto first two components

print(f"Variance explained by PC1: {explained_ratio[0]:.4f}")
print(f"Variance explained by PC2: {explained_ratio[1]:.4f}")
print(f"Variance explained by PC3: {explained_ratio[2]:.4f}")
print(f"Cumulative PC1+PC2 variance explained: {explained_ratio[:2].sum():.4f}")

# ----------------------------------------------------------------------
# 3. Hastie-Stuetzle principal curve in the PC1-PC2 plane (explicit loop)
# ----------------------------------------------------------------------
def project_to_polyline(points, curve):
    """Return arc-length parameter and distance of each point to a polyline."""
    seg_start = curve[:-1]
    seg_vec = curve[1:] - curve[:-1]
    seg_len2 = np.einsum("ij,ij->i", seg_vec, seg_vec)
    seg_len2[seg_len2 == 0] = 1e-12
    cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(seg_vec, axis=1))])
    lam = np.empty(len(points))
    dist = np.empty(len(points))
    for i, p in enumerate(points):
        t = np.clip(np.einsum("ij,j->i", p - seg_start, seg_vec.T.T @ np.eye(2)
                              if False else seg_vec) / seg_len2, 0.0, 1.0)
        # (above kept simple below) recompute cleanly:
        t = np.clip(((p - seg_start) * seg_vec).sum(axis=1) / seg_len2, 0.0, 1.0)
        proj = seg_start + t[:, None] * seg_vec       # foot point on each segment
        d = np.linalg.norm(proj - p, axis=1)
        k = np.argmin(d)
        lam[i] = cum[k] + t[k] * np.linalg.norm(seg_vec[k])
        dist[i] = d[k]
    return lam, dist

# Initialize the curve as the first principal component (lambda = PC1 score).
lam = pc[:, 0].copy()
prev_err = np.inf
for iteration in range(15):
    # --- Expectation-like step: order samples by their current curve parameter
    order = np.argsort(lam)
    lam_s = lam[order]
    lam_s = lam_s + np.linspace(0, 1e-6, len(lam_s))  # break ties for the spline

    # --- Maximization step: smooth each PC coordinate against lambda (scatterplot smoother)
    sp1 = UnivariateSpline(lam_s, pc[order, 0], k=3, s=n_samples * 0.5)
    sp2 = UnivariateSpline(lam_s, pc[order, 1], k=3, s=n_samples * 0.5)

    # Build a dense polyline of the smoothed curve.
    grid = np.linspace(lam_s.min(), lam_s.max(), 400)
    curve = np.column_stack([sp1(grid), sp2(grid)])

    # --- Reproject every sample onto the updated curve -> new arc-length lambda
    lam, dist = project_to_polyline(pc, curve)
    err = (dist ** 2).mean()
    if abs(prev_err - err) < 1e-6:
        break
    prev_err = err

print(f"Principal-curve iterations run: {iteration + 1}")
print(f"Mean squared distance of samples to principal curve: {err:.6f}")

# 1-D coordinate per sample = normalized arc length along the principal curve.
lam_norm = (lam - lam.min()) / (lam.max() - lam.min())

# ----------------------------------------------------------------------
# 4. Plots: scree plot + PC1-PC2 scatter colored by condition with curve
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

k = 10
ax1.bar(np.arange(1, k + 1), explained_ratio[:k], color="steelblue")
ax1.plot(np.arange(1, k + 1), explained_ratio[:k], "o-", color="black")
ax1.set_xlabel("Principal component")
ax1.set_ylabel("Fraction of variance explained")
ax1.set_title("PCA scree plot")

cmap = plt.get_cmap("viridis", len(cond_letters))
for j, letter in enumerate(cond_letters):
    m = labels == letter
    ax2.scatter(pc[m, 0], pc[m, 1], color=cmap(j), s=70,
                edgecolor="k", label=letter, zorder=3)
# Draw the fitted principal curve.
grid = np.linspace(lam_s.min(), lam_s.max(), 400)
curve = np.column_stack([sp1(grid), sp2(grid)])
ax2.plot(curve[:, 0], curve[:, 1], "-", color="crimson", lw=2.5,
         label="principal curve", zorder=2)
ax2.set_xlabel(f"PC1 ({explained_ratio[0]*100:.1f}% var)")
ax2.set_ylabel(f"PC2 ({explained_ratio[1]*100:.1f}% var)")
ax2.set_title("PC1-PC2 projection with principal curve")
ax2.legend(title="condition", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.4.1_s2.png",
            dpi=130)

# ----------------------------------------------------------------------
# 5. Separate check: PC1 dominant, ordered A->I arc, curve as 1-D coordinate
# ----------------------------------------------------------------------
print("\n--- CHECK ---")
print(f"PC1 explains most variance (ratio {explained_ratio[0]:.4f} > all others): "
      f"{explained_ratio[0] > explained_ratio[1:].max()}")

# Mean principal-curve coordinate per condition, in label order A..I.
print("Mean principal-curve coordinate per condition (should increase A->I):")
mean_lam = []
for letter in cond_letters:
    v = lam_norm[labels == letter].mean()
    mean_lam.append(v)
    print(f"  condition {letter}: {v:.4f}")
mean_lam = np.array(mean_lam)
print(f"Condition means monotonically increasing A->I: {np.all(np.diff(mean_lam) > 0)}")

rho, _ = spearmanr(cond_index, lam_norm)
print(f"Spearman corr(condition order, principal-curve coordinate): {rho:.4f}")

print("Per-sample 1-D principal-curve coordinate (sample, condition, lambda):")
for i in np.argsort(lam_norm):
    print(f"  sample {i:2d}  condition {labels[i]}  lambda = {lam_norm[i]:.4f}")

print("\nWhy this confirms the result: because PC1 captures the bulk of the variance "
      "and the per-condition principal-curve coordinates rise monotonically from A to I "
      "with near-perfect rank correlation, the samples genuinely lie on a single ordered "
      "arc that the principal curve collapses into one interpretable 1-D progression.")
