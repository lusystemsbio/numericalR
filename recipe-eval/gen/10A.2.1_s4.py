import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Recreate the 100x2 "strip" data from 10A.1: points scattered
# along a line (a regression-style cloud stretched along a direction).
# ---------------------------------------------------------------
rng = np.random.default_rng(0)
n = 100
x = rng.uniform(0, 10, n)                 # predictor spread along the strip
y = 2.0 + 0.8 * x + rng.normal(0, 0.5, n) # narrow scatter about a line
X = np.column_stack([x, y])               # samples in rows, features in columns

# ---------------------------------------------------------------
# Center only (subtract each column's mean); no scaling.
# ---------------------------------------------------------------
mean = X.mean(axis=0)
Xc = X - mean

# ---------------------------------------------------------------
# Route 1: eigen-decomposition of the covariance matrix.
# ---------------------------------------------------------------
# Sample covariance with the (n-1) denominator.
cov = (Xc.T @ Xc) / (n - 1)
# Symmetric matrix -> use eigh; eigenvalues ascending.
eigvals, eigvecs = np.linalg.eigh(cov)
# Sort descending so PC1 is first.
order = np.argsort(eigvals)[::-1]
eigvals = eigvals[order]
eigvecs = eigvecs[:, order]
var_eig = eigvals                         # component variances = eigenvalues
pc1_eig = eigvecs[:, 0]

# ---------------------------------------------------------------
# Route 2: SVD of the centered data. Xc = U S V^T.
# Columns of V are the principal directions; variances = s^2/(n-1).
# ---------------------------------------------------------------
U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
var_svd = (s ** 2) / (n - 1)              # component variances from singular values
pc1_svd = Vt[0, :]

# Fix sign so the two PC1 vectors point the same way (sign is arbitrary).
if np.dot(pc1_eig, pc1_svd) < 0:
    pc1_svd = -pc1_svd

# ---------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------
variances_match = np.allclose(var_eig, var_svd)
frac_pc1 = var_eig[0] / var_eig.sum()     # fraction of variance carried by PC1
# Regression direction of 10A.1: slope of y on x -> direction (1, slope).
slope = np.polyfit(x, y, 1)[0]
reg_dir = np.array([1.0, slope])
reg_dir = reg_dir / np.linalg.norm(reg_dir)
pc1_use = pc1_eig / np.linalg.norm(pc1_eig)
if np.dot(pc1_use, reg_dir) < 0:
    pc1_use = -pc1_use
# Angle between PC1 and the regression direction.
angle_deg = np.degrees(np.arccos(np.clip(np.dot(pc1_use, reg_dir), -1, 1)))

# ---------------------------------------------------------------
# Print every numerical result.
# ---------------------------------------------------------------
print(f"Covariance-route variances (PC1, PC2): {var_eig[0]:.6f}, {var_eig[1]:.6f}")
print(f"SVD-route variances       (PC1, PC2): {var_svd[0]:.6f}, {var_svd[1]:.6f}")
print(f"Routes return identical variances (allclose): {variances_match}")
print(f"Max abs variance difference between routes: {np.max(np.abs(var_eig - var_svd)):.3e}")
print(f"PC1 direction (covariance route): [{pc1_eig[0]:.6f}, {pc1_eig[1]:.6f}]")
print(f"PC1 direction (SVD route):        [{pc1_svd[0]:.6f}, {pc1_svd[1]:.6f}]")
print(f"Fraction of total variance carried by PC1: {frac_pc1:.6f}")
print(f"Regression direction (10A.1), unit: [{reg_dir[0]:.6f}, {reg_dir[1]:.6f}]")
print(f"Angle between PC1 and regression direction (deg): {angle_deg:.4f}")

# ---------------------------------------------------------------
# Plot: strip data with PC1 direction drawn through the mean.
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.scatter(X[:, 0], X[:, 1], s=15, alpha=0.6, label="strip data")
scale = 2.0 * np.sqrt(var_eig[0])         # length ~ 2 std devs along PC1
t = np.array([-scale, scale])
plt.plot(mean[0] + t * pc1_use[0], mean[1] + t * pc1_use[1],
         "r-", lw=2, label="PC1 direction")
plt.scatter([mean[0]], [mean[1]], c="k", marker="x", s=80, label="mean")
plt.xlabel("feature 1 (x)")
plt.ylabel("feature 2 (y)")
plt.title("Strip data with PC1 direction")
plt.legend()
plt.axis("equal")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.2.1_s4.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# Because covariance eigenvalues equal s^2/(n-1) from the SVD of the same
# centered data, identical variances from both routes verify the
# implementation is algebraically correct, and PC1 holding nearly all the
# variance while aligning with the y-on-x slope shows PCA recovers the
# strip's regression direction.
# ---------------------------------------------------------------
print("Check meaning: identical variances from two independent algebraic routes "
      "confirm the from-scratch PCA is correct, and PC1 carrying ~all variance "
      "along the regression direction confirms it captures the strip's structure.")
