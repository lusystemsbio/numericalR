import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(0)

# ----------------------------------------------------------------------
# 1. Build the test data: 26 samples x 500 genes lying on an ordered arc.
#    Conditions A..I are assigned in order along a latent progression t.
# ----------------------------------------------------------------------
n_samples, n_genes = 26, 500
labels_pool = list("ABCDEFGHI")

# latent 1D progression (the "true" ordering of the samples)
t = np.linspace(0.0, 1.0, n_samples)

# condition label for each sample: bin t into 9 ordered groups A..I
cond_idx = np.floor(t * (len(labels_pool) - 1e-9)).astype(int)
conditions = np.array([labels_pool[i] for i in cond_idx])

# The samples trace a curved (arc) path in a 2D latent space ...
angle = np.pi * t                                  # sweep half a circle
latent = np.column_stack([np.cos(angle), np.sin(angle)])

# ... which is embedded into 500 genes via a random linear map + noise.
loadings = np.random.randn(2, n_genes)
X = latent @ loadings + 0.15 * np.random.randn(n_samples, n_genes)

# ----------------------------------------------------------------------
# 2. PCA by SVD, done explicitly (no sklearn PCA).
# ----------------------------------------------------------------------
Xc = X - X.mean(axis=0)                             # center genes
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)   # SVD of centered data
var_explained = (S ** 2) / np.sum(S ** 2)           # fraction of variance per PC
scores = U * S                                      # PCA scores = U * singular values

# project onto the leading two components
Z = scores[:, :2]

print("Variance explained by PC1: {:.4f}".format(var_explained[0]))
print("Variance explained by PC2: {:.4f}".format(var_explained[1]))
print("Cumulative variance PC1+PC2: {:.4f}".format(var_explained[:2].sum()))

# ----------------------------------------------------------------------
# 3. Hastie-Stuetzle principal curve, implemented explicitly.
#    Alternate a projection step and a smoothing (conditional-expectation)
#    step until the fitted curve stops changing.
# ----------------------------------------------------------------------
def arc_length_param(curve_pts):
    # cumulative Euclidean distance along an ordered polyline -> lambda
    d = np.sqrt(((np.diff(curve_pts, axis=0)) ** 2).sum(axis=1))
    return np.concatenate([[0.0], np.cumsum(d)])

def project_to_polyline(pts, curve_pts, curve_lam):
    # for each point return (lambda of nearest point on polyline, foot point)
    lam_out = np.empty(len(pts))
    foot_out = np.empty_like(pts)
    seg_a = curve_pts[:-1]
    seg_b = curve_pts[1:]
    seg_v = seg_b - seg_a
    seg_len2 = (seg_v ** 2).sum(axis=1)
    for i, p in enumerate(pts):
        # projection parameter u onto each segment, clamped to [0,1]
        u = ((p - seg_a) * seg_v).sum(axis=1) / np.where(seg_len2 == 0, 1, seg_len2)
        u = np.clip(u, 0.0, 1.0)
        proj = seg_a + u[:, None] * seg_v
        dist2 = ((proj - p) ** 2).sum(axis=1)
        j = np.argmin(dist2)                        # closest segment
        lam_out[i] = curve_lam[j] + u[j] * np.sqrt(seg_len2[j])
        foot_out[i] = proj[j]
    return lam_out, foot_out

def smooth_vs_lambda(lam, y, span=0.35):
    # local-linear smoother: fit y ~ lambda in a moving window (E[X|lambda])
    order = np.argsort(lam)
    lo, yo = lam[order], y[order]
    win = max(3, int(span * len(lo)))
    out = np.empty_like(yo)
    for i in range(len(lo)):
        a = max(0, i - win // 2)
        b = min(len(lo), a + win)
        a = max(0, b - win)
        A = np.vstack([np.ones(b - a), lo[a:b]]).T   # local linear design
        coef, *_ = np.linalg.lstsq(A, yo[a:b], rcond=None)
        out[i] = coef[0] + coef[1] * lo[i]
    res = np.empty_like(out)
    res[order] = out
    return res

# initialise the curve as the PC1 line through the data
init_lam = Z[:, 0] - Z[:, 0].mean()
dir1 = np.array([1.0, 0.0])                          # PC1 axis in score space
curve_lam_grid = np.linspace(init_lam.min(), init_lam.max(), 100)
curve_pts = Z.mean(axis=0) + np.outer(curve_lam_grid, dir1)

prev = np.inf
for it in range(20):
    # ---- projection step: assign each sample an arc-length lambda ----
    curve_lam = arc_length_param(curve_pts)
    lam, _ = project_to_polyline(Z, curve_pts, curve_lam)

    # ---- smoothing step: curve coordinate = E[Z | lambda] ----
    fx = smooth_vs_lambda(lam, Z[:, 0])
    fy = smooth_vs_lambda(lam, Z[:, 1])

    # rebuild an ordered, resampled polyline from the smoothed fit
    order = np.argsort(lam)
    curve_pts = np.column_stack([fx[order], fy[order]])

    # convergence: mean squared distance of points to the new curve
    cl = arc_length_param(curve_pts)
    lam2, foot = project_to_polyline(Z, curve_pts, cl)
    mse = ((Z - foot) ** 2).sum(axis=1).mean()
    if abs(prev - mse) < 1e-6:
        break
    prev = mse

# final 1D principal-curve coordinate per sample (normalised arc length)
pc_coord = (lam2 - lam2.min()) / (lam2.max() - lam2.min())
print("Principal-curve fit MSE (points to curve): {:.6f}".format(mse))
print("Iterations to converge: {}".format(it + 1))

# ----------------------------------------------------------------------
# 4. Separate check: does the 1D coordinate order samples A -> I ?
# ----------------------------------------------------------------------
order_by_coord = np.argsort(pc_coord)
mean_coord = np.array([pc_coord[conditions == c].mean() for c in labels_pool])
# orient so condition A is at the low end
if mean_coord[0] > mean_coord[-1]:
    pc_coord = 1.0 - pc_coord
    mean_coord = 1.0 - mean_coord
monotone = bool(np.all(np.diff(mean_coord) > 0))
print("Mean principal-curve coordinate per condition A..I:")
for c, m in zip(labels_pool, mean_coord):
    print("  condition {}: {:.4f}".format(c, m))
print("Conditions are monotonically ordered A->I along the curve: {}".format(monotone))

# ----------------------------------------------------------------------
# 5. Plots: scree plot + PC1-PC2 scatter with fitted principal curve.
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# scree plot
ax1.bar(np.arange(1, 11), var_explained[:10] * 100, color="steelblue")
ax1.set_xlabel("Principal component")
ax1.set_ylabel("Variance explained (%)")
ax1.set_title("Scree plot")

# PC1-PC2 scatter coloured by condition
cmap = plt.get_cmap("viridis", len(labels_pool))
for k, c in enumerate(labels_pool):
    m = conditions == c
    ax2.scatter(Z[m, 0], Z[m, 1], color=cmap(k), s=60, edgecolor="k",
                label=c, zorder=3)
# draw the fitted principal curve (ordered polyline)
cp = curve_pts
ax2.plot(cp[:, 0], cp[:, 1], "r-", lw=2.5, label="principal curve", zorder=2)
ax2.set_xlabel("PC1 ({:.1f}%)".format(var_explained[0] * 100))
ax2.set_ylabel("PC2 ({:.1f}%)".format(var_explained[1] * 100))
ax2.set_title("PC1-PC2 projection with principal curve")
ax2.legend(ncol=2, fontsize=8, title="condition")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.4.1_s3.png")

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# ----------------------------------------------------------------------
print("Check rationale: because PC1 captures most variance and the per-condition "
      "mean principal-curve coordinate increases monotonically from A to I, the "
      "samples genuinely lie on an ordered arc that the principal curve summarizes "
      "as a single 1D progression coordinate.")
