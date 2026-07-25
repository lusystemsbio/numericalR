import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# ----------------------------------------------------------------------
# 1. Build synthetic gene-expression data: 26 samples x 500 genes.
#    The samples are placed along a curved 1-D latent trajectory (an arc)
#    so that a principal curve should recover their ordering A -> I.
# ----------------------------------------------------------------------
rng = np.random.default_rng(0)
n_samples, n_genes = 26, 500

# latent progression parameter t in [0, 1], one value per sample (ordered)
t = np.linspace(0.0, 1.0, n_samples)

# two smooth "programs" that trace an arc in a 2-D latent space
u1 = np.cos(np.pi * t)          # bends from +1 -> -1
u2 = np.sin(np.pi * t)          # rises then falls -> arc shape
latent = np.column_stack([u1, u2])          # (26, 2) arc

# random loadings map the 2-D latent arc into 500 genes, plus noise
loadings = rng.normal(size=(2, n_genes))
X = latent @ loadings + 0.15 * rng.normal(size=(n_samples, n_genes))

# label the 26 samples by condition A..I (9 conditions along the arc)
conds = np.array([chr(ord('A') + int(i * 9 / n_samples)) for i in range(n_samples)])

# ----------------------------------------------------------------------
# 2. PCA implemented explicitly via SVD of the centered data matrix.
# ----------------------------------------------------------------------
Xc = X - X.mean(axis=0)                       # center each gene
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
scores = U * S                                # sample coordinates in PC space
var_explained = (S**2) / np.sum(S**2)         # fraction of variance per PC

scores2d = scores[:, :2]                       # project onto PC1-PC2
print("Variance explained by PC1:", var_explained[0])
print("Variance explained by PC2:", var_explained[1])
print("Variance explained by PC3:", var_explained[2])
print("Cumulative variance PC1+PC2:", var_explained[0] + var_explained[1])

# ----------------------------------------------------------------------
# 3. Hastie-Stuetzle principal curve, implemented explicitly.
#    Alternate: (a) projection of points onto the current curve to get a
#    1-D arc-length coordinate lambda, then (b) smoothing each coordinate
#    against lambda to update the curve. Repeat to convergence.
# ----------------------------------------------------------------------
def arc_length(curve_pts):
    # cumulative Euclidean distance along an ordered polyline
    d = np.sqrt(np.sum(np.diff(curve_pts, axis=0)**2, axis=1))
    return np.concatenate([[0.0], np.cumsum(d)])

def project_to_polyline(pts, curve_pts, curve_lam):
    # for each data point find nearest point on the polyline; return its
    # arc-length coordinate lambda and the projected coordinates
    lam = np.zeros(len(pts))
    proj = np.zeros_like(pts)
    for i, p in enumerate(pts):
        best_d2, best_lam, best_proj = np.inf, 0.0, curve_pts[0]
        for j in range(len(curve_pts) - 1):
            a, b = curve_pts[j], curve_pts[j + 1]
            ab = b - a
            L2 = ab @ ab
            s = 0.0 if L2 == 0 else np.clip((p - a) @ ab / L2, 0.0, 1.0)
            foot = a + s * ab                  # closest point on this segment
            d2 = (p - foot) @ (p - foot)
            if d2 < best_d2:
                best_d2 = d2
                best_lam = curve_lam[j] + s * (curve_lam[j + 1] - curve_lam[j])
                best_proj = foot
        lam[i] = best_lam
        proj[i] = best_proj
    return lam, proj

# (a) initialize the curve as the PC1 line through the data
direction = np.array([1.0, 0.0])              # PC1 axis in PC space
lam = scores2d @ direction                     # initial 1-D coordinate
lam = lam - lam.min()

prev_criterion = np.inf
for iteration in range(15):
    order = np.argsort(lam)                     # order samples along the curve
    lam_s = lam[order]
    # (b) smoothing step: each PC coordinate as a smooth function of lambda
    sx = UnivariateSpline(lam_s, scores2d[order, 0], k=3, s=len(lam) * 0.02)
    sy = UnivariateSpline(lam_s, scores2d[order, 1], k=3, s=len(lam) * 0.02)
    # rebuild a fine polyline representation of the smoothed curve
    grid = np.linspace(lam_s.min(), lam_s.max(), 200)
    curve_pts = np.column_stack([sx(grid), sy(grid)])
    curve_lam = arc_length(curve_pts)          # reparameterize by arc length
    # (a) projection step: reproject data onto the updated curve
    lam, proj = project_to_polyline(scores2d, curve_pts, curve_lam)
    criterion = np.mean(np.sum((scores2d - proj)**2, axis=1))  # residual var
    if abs(prev_criterion - criterion) < 1e-6:
        break
    prev_criterion = criterion

print("Principal curve iterations run:", iteration + 1)
print("Principal curve residual variance:", criterion)

# per-sample 1-D principal-curve coordinate (normalized to [0,1])
pc_coord = (lam - lam.min()) / (lam.max() - lam.min())

# ----------------------------------------------------------------------
# 4. Check: does the 1-D principal-curve coordinate follow A -> I order?
#    Compare curve ordering against the true condition ordering.
# ----------------------------------------------------------------------
cond_rank = np.array([ord(c) - ord('A') for c in conds])
# Spearman-style rank correlation between condition order and curve coord
r = np.corrcoef(np.argsort(np.argsort(pc_coord)),
                np.argsort(np.argsort(cond_rank)))[0, 1]
print("Rank correlation (curve coordinate vs condition A..I):", abs(r))
for c in sorted(set(conds)):
    print(f"Mean principal-curve coordinate for condition {c}:",
          pc_coord[conds == c].mean())

# ----------------------------------------------------------------------
# 5. Plots: scree plot + PC1-PC2 scatter colored by condition with curve.
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

ax1.bar(np.arange(1, 11), var_explained[:10], color="steelblue")
ax1.set_xlabel("Principal component")
ax1.set_ylabel("Fraction of variance explained")
ax1.set_title("PCA scree plot")

uniq = sorted(set(conds))
cmap = plt.get_cmap("viridis", len(uniq))
for k, c in enumerate(uniq):
    m = conds == c
    ax2.scatter(scores2d[m, 0], scores2d[m, 1], color=cmap(k), s=60,
                edgecolor="k", label=c, zorder=3)
# draw the fitted principal curve, ordered by arc length
o = np.argsort(curve_lam)
ax2.plot(curve_pts[o, 0], curve_pts[o, 1], "r-", lw=2.5,
         label="principal curve", zorder=2)
ax2.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}% var)")
ax2.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}% var)")
ax2.set_title("PC1-PC2 projection with principal curve")
ax2.legend(title="Condition", fontsize=8, ncol=2)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.4.1_s5.png", dpi=130)

# One-sentence explanation of why the check confirms the result:
print("Explanation: A high rank correlation between the principal-curve "
      "coordinate and the A->I condition order, together with PC1 holding "
      "most of the variance, confirms the samples lie on an ordered 1-D arc "
      "that the principal curve faithfully summarizes as a single coordinate.")
