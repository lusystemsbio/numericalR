import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Recreate the 100x2 "strip" data from 10A.1: points scattered along a line ---
rng = np.random.default_rng(0)
n = 100
x = rng.uniform(0, 10, n)              # spread along the strip
y = 2.0 * x + 1.0 + rng.normal(0, 1.0, n)  # a line (slope 2) plus small noise
X = np.column_stack([x, y])            # samples in rows, features in columns

# --- Center the data (centering only, no scaling) ---
mean = X.mean(axis=0)                  # per-feature mean
Xc = X - mean                          # centered data matrix

# === Route 1: eigen-decomposition of the covariance matrix ===
# Covariance uses (n-1) denominator (sample covariance)
C = (Xc.T @ Xc) / (n - 1)             # 2x2 covariance matrix
eigvals, eigvecs = np.linalg.eigh(C)   # eigh: symmetric, ascending eigenvalues
order = np.argsort(eigvals)[::-1]      # sort descending -> PC1 first
eig_variances = eigvals[order]         # component variances
eig_directions = eigvecs[:, order]     # columns are PC directions

# === Route 2: singular value decomposition of the centered data ===
# Xc = U S Vt ; variances = S^2 / (n-1); PC directions are rows of Vt
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
svd_variances = (S ** 2) / (n - 1)     # component variances from singular values
svd_directions = Vt.T                  # columns are PC directions

# --- Report component variances from both routes ---
print(f"Route 1 (covariance eigen) variances: PC1={eig_variances[0]:.6f}, PC2={eig_variances[1]:.6f}")
print(f"Route 2 (SVD of centered)   variances: PC1={svd_variances[0]:.6f}, PC2={svd_variances[1]:.6f}")

# --- Check: the two routes return identical variances ---
max_var_diff = np.max(np.abs(eig_variances - svd_variances))
print(f"Max absolute difference in variances between routes: {max_var_diff:.3e}")
print(f"Two routes agree on variances (within tolerance): {np.allclose(eig_variances, svd_variances)}")

# --- Check: PC1 carries almost all the variance ---
frac_pc1 = eig_variances[0] / eig_variances.sum()
print(f"Fraction of total variance carried by PC1: {frac_pc1:.6f}")

# --- Check: PC1 recovers the regression direction (slope) of 10A.1 ---
pc1 = eig_directions[:, 0]
if pc1[0] < 0:                         # fix sign so it points to +x for readability
    pc1 = -pc1
pc1_slope = pc1[1] / pc1[0]            # slope implied by PC1 direction
print(f"PC1 direction vector: [{pc1[0]:.6f}, {pc1[1]:.6f}]")
print(f"Slope implied by PC1 direction: {pc1_slope:.6f}")
print(f"True line slope used to generate strip data: {2.0:.6f}")

# --- Plot: strip data with PC1 direction drawn through the mean ---
plt.figure(figsize=(7, 6))
plt.scatter(X[:, 0], X[:, 1], s=15, alpha=0.6, label="strip data")
# draw PC1 as a line segment through the mean, scaled by PC1 std for visibility
scale = 3.0 * np.sqrt(eig_variances[0])
p0 = mean - scale * pc1
p1 = mean + scale * pc1
plt.plot([p0[0], p1[0]], [p0[1], p1[1]], "r-", lw=2, label="PC1 direction")
plt.scatter(*mean, color="k", zorder=5, label="mean")
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.title("Strip data with PC1 direction")
plt.legend()
plt.axis("equal")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.2.1_s5.png")

# --- One-sentence explanation of why the check confirms the result ---
print("Explanation: Because covariance = Xc^T Xc/(n-1) and SVD gives Xc = U S Vt with "
      "eigenvectors V and eigenvalues S^2, the two routes are algebraically the same "
      "decomposition, so identical variances plus a dominant PC1 aligned with the "
      "generating slope confirm PCA correctly recovers the strip's regression direction.")
