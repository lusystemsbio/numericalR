import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Recreate the 100x2 "strip" data from 10A.1: points scattered
# tightly along a line so that they form a thin strip.
# ---------------------------------------------------------------
rng = np.random.default_rng(0)
n = 100
x = np.linspace(0, 10, n)                 # coordinate along the strip
y = 0.5 * x + rng.normal(0, 0.3, n)       # line (slope 0.5) + small noise
X = np.column_stack([x, y])               # samples in rows, features in columns

# ---------------------------------------------------------------
# Center only (subtract the column means); no scaling.
# ---------------------------------------------------------------
mean = X.mean(axis=0)
Xc = X - mean

# ---------------------------------------------------------------
# Route 1: eigen-decomposition of the covariance matrix.
# ---------------------------------------------------------------
# Covariance matrix of the centered data (divide by n-1).
cov = (Xc.T @ Xc) / (n - 1)
# Eigenvalues (component variances) and eigenvectors (directions).
eigvals, eigvecs = np.linalg.eigh(cov)
# eigh returns ascending order; sort descending so PC1 is first.
order = np.argsort(eigvals)[::-1]
eig_variances = eigvals[order]
eig_vectors = eigvecs[:, order]

# ---------------------------------------------------------------
# Route 2: singular value decomposition of the centered data.
# ---------------------------------------------------------------
# Xc = U S Vt ; columns of V (rows of Vt) are the principal directions.
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
# The variance along each PC is s^2 / (n-1).
svd_variances = (S ** 2) / (n - 1)
svd_vectors = Vt.T

# ---------------------------------------------------------------
# Report the component variances from both routes.
# ---------------------------------------------------------------
print("Eigen route component variances (PC1, PC2):", eig_variances[0], eig_variances[1])
print("SVD route component variances   (PC1, PC2):", svd_variances[0], svd_variances[1])

# ---------------------------------------------------------------
# Check: the two routes return identical variances.
# ---------------------------------------------------------------
max_var_diff = np.max(np.abs(eig_variances - svd_variances))
print("Max absolute difference in variances between routes:", max_var_diff)
print("Routes agree on variances (within 1e-8):", bool(max_var_diff < 1e-8))

# ---------------------------------------------------------------
# Check: PC1 carries almost all the variance.
# ---------------------------------------------------------------
total_var = eig_variances.sum()
pc1_fraction = eig_variances[0] / total_var
print("Total variance:", total_var)
print("Fraction of variance carried by PC1:", pc1_fraction)

# ---------------------------------------------------------------
# Check: PC1 direction recovers the regression direction of 10A.1.
# Sign of an eigenvector is arbitrary; orient PC1 to have positive x.
# ---------------------------------------------------------------
pc1 = eig_vectors[:, 0]
if pc1[0] < 0:
    pc1 = -pc1
pc1_slope = pc1[1] / pc1[0]               # slope implied by PC1 direction
print("PC1 direction (unit vector):", pc1[0], pc1[1])
print("Slope implied by PC1:", pc1_slope)
print("True generating slope (10A.1 regression direction):", 0.5)

# ---------------------------------------------------------------
# Plot: strip data with the PC1 direction drawn through the mean.
# ---------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.scatter(X[:, 0], X[:, 1], s=15, alpha=0.6, label="strip data")
# Scale the arrow length by PC1 standard deviation for visibility.
scale = 3.0 * np.sqrt(eig_variances[0])
plt.plot([mean[0] - scale * pc1[0], mean[0] + scale * pc1[0]],
         [mean[1] - scale * pc1[1], mean[1] + scale * pc1[1]],
         color="red", linewidth=2, label="PC1 direction")
plt.scatter([mean[0]], [mean[1]], color="black", zorder=5, label="mean")
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.title("Strip data with PC1 direction")
plt.legend()
plt.axis("equal")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.2.1_s2.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("Explanation: Because the covariance eigen-decomposition and the SVD of the "
      "centered data are algebraically the same factorization (cov = Vt.T diag(s^2/(n-1)) Vt), "
      "their identical variances plus a PC1 that captures nearly all variance and matches the "
      "regression slope confirm the PCA is computed correctly.")
