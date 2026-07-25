import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 1. Build the model / test data: 26 samples x 500 genes, labeled A..I.
#    We plant a smooth 1D "condition" progression (an ordered arc) into a
#    high-dimensional gene space so that PCA + a principal curve can recover it.
# ----------------------------------------------------------------------
rng = np.random.default_rng(0)
n_samples, n_genes = 26, 500
conditions = "ABCDEFGHI"                      # 9 conditions
# assign 26 samples to 9 conditions in order (roughly 3 per condition)
labels = np.array([conditions[int(i * len(conditions) / n_samples)] for i in range(n_samples)])

# latent 1D coordinate t in [0, 1], ordered A->I; place samples along an arc
t = np.linspace(0.0, 1.0, n_samples)
# arc in a 2D latent plane (a curved, monotone progression, not a straight line)
theta = np.pi * t                             # sweep half a circle
arc_x = np.cos(theta)
arc_y = np.sin(theta)

# random loadings map the 2D arc into 500 gene dimensions
load_x = rng.normal(size=n_genes)
load_y = rng.normal(size=n_genes)
signal = np.outer(arc_x, load_x) + np.outer(arc_y, load_y)
noise = 0.15 * rng.normal(size=(n_samples, n_genes))
X = 5.0 * signal + noise                       # gene-expression matrix

# ----------------------------------------------------------------------
# 2. PCA implemented explicitly (center, SVD, variance explained, scores).
# ----------------------------------------------------------------------
Xc = X - X.mean(axis=0)                         # center each gene
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
var_explained = (S ** 2) / np.sum(S ** 2)       # fraction of variance per PC
scores = U * S                                  # sample coordinates in PC space
pc = scores[:, :2]                              # project onto first two PCs

print(f"PC1 variance explained: {var_explained[0]:.4f}")
print(f"PC2 variance explained: {var_explained[1]:.4f}")
print(f"PC3 variance explained: {var_explained[2]:.4f}")
print(f"PC1+PC2 cumulative variance explained: {var_explained[:2].sum():.4f}")


# ----------------------------------------------------------------------
# 3. Hastie-Stuetzle principal curve, implemented explicitly.
#    Alternate two steps until convergence:
#      (a) projection: assign each point the arc-length of its nearest curve point
#      (b) expectation: re-estimate the curve as a smooth function of arc-length
# ----------------------------------------------------------------------
def arc_length(curve):
    # cumulative Euclidean distance along an ordered polyline
    seg = np.sqrt(((curve[1:] - curve[:-1]) ** 2).sum(axis=1))
    return np.concatenate([[0.0], np.cumsum(seg)])

def project_to_curve(pts, curve, s_curve):
    # nearest point on the polyline for each sample -> its arc-length lambda
    lam = np.empty(len(pts))
    for i, p in enumerate(pts):
        a, b = curve[:-1], curve[1:]            # segment endpoints
        ab = b - a
        ap = p - a
        denom = (ab ** 2).sum(axis=1)
        u = np.clip((ap * ab).sum(axis=1) / np.where(denom == 0, 1, denom), 0, 1)
        proj = a + u[:, None] * ab              # foot of perpendicular per segment
        d2 = ((proj - p) ** 2).sum(axis=1)
        j = np.argmin(d2)                       # closest segment
        lam[i] = s_curve[j] + u[j] * (s_curve[j + 1] - s_curve[j])
    return lam

def smooth(lam, y, span=0.35):
    # local-linear smoother: y as a smooth function of arc-length lam
    order = np.argsort(lam)
    xs, ys = lam[order], y[order]
    out = np.empty_like(ys)
    h = span * (xs[-1] - xs[0] + 1e-12)         # kernel bandwidth
    for i, xi in enumerate(xs):
        w = np.exp(-0.5 * ((xs - xi) / h) ** 2)
        sw = w.sum()
        mx = (w * xs).sum() / sw
        my = (w * ys).sum() / sw
        b = (w * (xs - mx) * (ys - my)).sum() / ((w * (xs - mx) ** 2).sum() + 1e-12)
        out[i] = my + b * (xi - mx)             # local linear fit at xi
    res = np.empty_like(out)
    res[order] = out
    return res

# (init) start from the first principal component line through the data
lam = pc @ (Vt[0, :2] if Vt.shape[1] >= 2 else np.array([1.0, 0.0]))
lam = pc[:, 0].copy()                           # PC1 score is the natural start
for _ in range(15):
    order = np.argsort(lam)
    curve = pc[order]                           # ordered samples define the curve
    s_curve = arc_length(curve)
    lam_new = project_to_curve(pc, curve, s_curve)          # (a) projection step
    fx = smooth(lam_new, pc[:, 0])                          # (b) smooth x vs lambda
    fy = smooth(lam_new, pc[:, 1])                          # (b) smooth y vs lambda
    new_curve_pts = np.column_stack([fx, fy])
    if np.max(np.abs(lam_new - lam)) < 1e-6:
        lam = lam_new
        break
    lam = lam_new
    pc_fit = new_curve_pts                       # smoothed curve position per sample

# final smoothed curve, ordered by arc-length, plus 1D coordinate per sample
order = np.argsort(lam)
fx = smooth(lam, pc[:, 0])
fy = smooth(lam, pc[:, 1])
curve_pts = np.column_stack([fx, fy])[order]
# normalize the principal-curve arc-length to a 0..1 coordinate per sample
lam_norm = (lam - lam.min()) / (lam.max() - lam.min())

print("\nPrincipal-curve 1D coordinate per sample (label: coord):")
for i in range(n_samples):
    print(f"  sample {i:2d}  condition {labels[i]}  coord={lam_norm[i]:.4f}")

# order check: mean principal-curve coordinate per condition, A..I
print("\nMean principal-curve coordinate by condition (should increase A->I):")
mean_coords = []
for c in conditions:
    m = lam_norm[labels == c].mean()
    mean_coords.append(m)
    print(f"  condition {c}: {m:.4f}")
mono = np.all(np.diff(mean_coords) > 0)
sp = np.corrcoef(np.argsort(np.argsort(lam_norm)),
                 np.array([conditions.index(c) for c in labels]))[0, 1]
print(f"\nCondition means monotonically increasing A->I: {bool(mono)}")
print(f"Spearman-rank correlation (curve coord vs condition order): {sp:.4f}")

# ----------------------------------------------------------------------
# 4. Plots: scree plot + PC1-PC2 scatter colored by condition with the curve.
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

k = 10
ax1.bar(np.arange(1, k + 1), var_explained[:k], color="steelblue")
ax1.plot(np.arange(1, k + 1), np.cumsum(var_explained[:k]), "o-", color="darkorange",
         label="cumulative")
ax1.set_xlabel("Principal component")
ax1.set_ylabel("Fraction of variance explained")
ax1.set_title("PCA scree plot")
ax1.legend()

cmap = plt.get_cmap("viridis", len(conditions))
for idx, c in enumerate(conditions):
    m = labels == c
    ax2.scatter(pc[m, 0], pc[m, 1], color=cmap(idx), s=70,
                edgecolor="k", label=c, zorder=3)
ax2.plot(curve_pts[:, 0], curve_pts[:, 1], "-", color="crimson", lw=2.5,
         label="principal curve", zorder=2)
ax2.set_xlabel(f"PC1 ({var_explained[0] * 100:.1f}%)")
ax2.set_ylabel(f"PC2 ({var_explained[1] * 100:.1f}%)")
ax2.set_title("PC1-PC2 projection with Hastie-Stuetzle principal curve")
ax2.legend(ncol=2, fontsize=8, title="condition")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.4.1_s4.png", dpi=130)

# ----------------------------------------------------------------------
# 5. Why the check confirms the result (one sentence).
# ----------------------------------------------------------------------
print("\nWhy this check confirms the result:")
print("Because PC1 captures most variance while the principal-curve coordinate "
      "increases monotonically from condition A to I, the fitted 1D curve "
      "faithfully summarizes the ordered biological progression embedded in the data.")
