import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Reconstruct the 100x2 "strip" data from 10A.1 ---
# Points scattered along a line: x roughly uniform, y = line + small noise.
rng = np.random.default_rng(0)
n = 100
x = rng.uniform(0, 10, n)                 # spread along the strip
y = 2.0 + 0.5 * x + rng.normal(0, 0.3, n)  # a line plus small vertical scatter
X = np.column_stack([x, y])               # samples in rows, features in columns

# --- Center only (no scaling): subtract the column means ---
means = X.mean(axis=0)
Xc = X - means

# === Route 1: eigen-decomposition of the covariance matrix ===
# Covariance with (n-1) denominator; features are along columns.
C = (Xc.T @ Xc) / (n - 1)
eigvals, eigvecs = np.linalg.eigh(C)   # eigh: symmetric, ascending eigenvalues
order = np.argsort(eigvals)[::-1]      # sort descending so PC1 is first
eig_var = eigvals[order]               # component variances from eigen route
eig_dirs = eigvecs[:, order]           # columns are the PC directions

# === Route 2: singular value decomposition of the centered data ===
# Xc = U S Vt; the rows of Vt are the principal directions.
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
svd_var = (S ** 2) / (n - 1)           # component variances from SVD route
svd_dirs = Vt.T                        # columns are the PC directions

# --- Align sign of directions so the two routes are comparable ---
for k in range(eig_dirs.shape[1]):
    if np.dot(eig_dirs[:, k], svd_dirs[:, k]) < 0:
        svd_dirs[:, k] *= -1

# PC1 direction (from SVD route) and its variance share
pc1 = svd_dirs[:, 0]
var_share_pc1 = svd_var[0] / svd_var.sum()

# --- Print numerical results ---
print("Component variances (eigen route) PC1:", eig_var[0])
print("Component variances (eigen route) PC2:", eig_var[1])
print("Component variances (SVD route)   PC1:", svd_var[0])
print("Component variances (SVD route)   PC2:", svd_var[1])
print("Max abs difference between routes' variances:", np.max(np.abs(eig_var - svd_var)))
print("Routes agree (allclose):", np.allclose(eig_var, svd_var))
print("PC1 direction (SVD route):", pc1[0], pc1[1])
print("PC1 slope (dy/dx):", pc1[1] / pc1[0])
print("PC1 fraction of total variance:", var_share_pc1)
print("Regression slope of 10A.1 (y on x) for comparison:",
      np.polyfit(x, y, 1)[0])

# --- Plot: strip data with PC1 direction drawn through the mean ---
scale = 3 * np.sqrt(svd_var[0])   # length ~3 std devs along PC1
plt.figure(figsize=(6, 6))
plt.scatter(X[:, 0], X[:, 1], s=18, alpha=0.6, label="strip data")
plt.plot([means[0] - scale * pc1[0], means[0] + scale * pc1[0]],
         [means[1] - scale * pc1[1], means[1] + scale * pc1[1]],
         "r-", lw=2, label="PC1 direction")
plt.scatter([means[0]], [means[1]], c="k", s=40, marker="x", label="mean")
plt.gca().set_aspect("equal", adjustable="datalim")
plt.xlabel("feature 1 (x)")
plt.ylabel("feature 2 (y)")
plt.title("Strip data with PC1 direction")
plt.legend()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the covariance matrix equals (Xc^T Xc)/(n-1) and the "
      "SVD gives Xc = U S V^T, the squared singular values over (n-1) are exactly "
      "the covariance eigenvalues, so matching variances from both routes confirms "
      "the from-scratch PCA is correct.")
