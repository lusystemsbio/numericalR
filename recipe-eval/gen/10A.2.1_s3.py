import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Recreate the 100x2 "strip" data from 10A.1: points scattered along a line ---
rng = np.random.default_rng(0)
n = 100
x = rng.uniform(-3, 3, n)              # spread along the strip
y = 2.0 * x + rng.normal(0, 0.5, n)    # a line (slope 2) with small scatter
X = np.column_stack([x, y])            # samples in rows, features in columns

# --- Center the data (centering only, no scaling) ---
mean = X.mean(axis=0)                  # per-feature mean
Xc = X - mean                          # centered data matrix

# --- Route 1: eigen-decomposition of the covariance matrix ---
# Covariance uses (n-1) in the denominator (sample covariance)
C = (Xc.T @ Xc) / (n - 1)             # 2x2 covariance matrix
eigvals, eigvecs = np.linalg.eigh(C)  # eigh: symmetric matrix, ascending eigenvalues
order = np.argsort(eigvals)[::-1]      # sort descending so PC1 is first
eig_var = eigvals[order]               # component variances from covariance route
eig_dirs = eigvecs[:, order]           # eigenvectors = principal directions

# --- Route 2: SVD of the centered data ---
# Xc = U S Vt ; columns of V are principal directions, variances = s^2/(n-1)
U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
svd_var = (s ** 2) / (n - 1)           # component variances from SVD route
svd_dirs = Vt.T                        # principal directions (columns)

# --- PC1 direction (fix sign so both routes point the same way for reporting) ---
pc1_eig = eig_dirs[:, 0]
pc1_svd = svd_dirs[:, 0]
if np.dot(pc1_eig, pc1_svd) < 0:
    pc1_svd = -pc1_svd
# orient PC1 to have positive x-component for a stable slope reading
if pc1_eig[0] < 0:
    pc1_eig = -pc1_eig
    pc1_svd = -pc1_svd

# --- Regression direction from 10A.1 for comparison (slope of y on x) ---
reg_slope = np.polyfit(x, y, 1)[0]
pc1_slope = pc1_eig[1] / pc1_eig[0]    # slope implied by PC1 direction

# --- Fraction of variance carried by PC1 ---
frac_pc1_eig = eig_var[0] / eig_var.sum()
frac_pc1_svd = svd_var[0] / svd_var.sum()

# --- Check: do the two routes agree on variances? ---
variances_agree = np.allclose(eig_var, svd_var)

# --- Print every numerical result ---
print("Covariance-route component variances (PC1, PC2):", eig_var[0], eig_var[1])
print("SVD-route component variances (PC1, PC2):       ", svd_var[0], svd_var[1])
print("Max abs difference between routes' variances:   ", np.max(np.abs(eig_var - svd_var)))
print("Two routes return identical variances (allclose):", variances_agree)
print("PC1 direction (covariance route):                ", pc1_eig[0], pc1_eig[1])
print("PC1 direction (SVD route):                       ", pc1_svd[0], pc1_svd[1])
print("PC1 fraction of total variance (covariance):     ", frac_pc1_eig)
print("PC1 fraction of total variance (SVD):            ", frac_pc1_svd)
print("PC1 implied slope (dy/dx):                       ", pc1_slope)
print("Regression slope from 10A.1 (y on x):            ", reg_slope)

# --- Plot: strip data with PC1 direction drawn through the mean ---
t = np.linspace(-4, 4, 2)
line = mean[:, None] + pc1_eig[:, None] * (t[None, :] * s[0] / np.sqrt(n - 1))
plt.figure(figsize=(6, 6))
plt.scatter(X[:, 0], X[:, 1], s=15, alpha=0.6, label="strip data")
plt.plot(line[0], line[1], "r-", lw=2, label="PC1 direction")
plt.scatter([mean[0]], [mean[1]], color="k", zorder=5, label="mean")
plt.axis("equal")
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.title("Strip data with PC1 direction")
plt.legend()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10A.2.1_s3.png")

# One sentence: matching variances confirm the result because the eigen-decomposition of
# the covariance matrix and the SVD of the centered data are mathematically the same
# factorization (C = Xc^T Xc /(n-1) has eigenvalues s^2/(n-1)), so identical outputs
# verify the from-scratch implementation is correct, with PC1's near-total variance and
# strip-aligned slope showing it recovers the 10A.1 regression direction.
print("Check explanation: identical variances from both routes confirm correctness because "
      "eigen-decomposition of Xc^T Xc/(n-1) and SVD of Xc are the same factorization; "
      "PC1 carrying nearly all variance and matching the regression slope shows it points along the strip.")
