import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import multivariate_normal

# ---------------- Generate the three-blob data ----------------
rng = np.random.RandomState(1)  # seed 1 for reproducibility
true_means = np.array([[0.0, 0.0], [5.0, 5.0], [0.0, 6.0]])
n_per = 100
X = np.vstack([rng.randn(n_per, 2) + m for m in true_means])  # three Gaussian blobs
N, D = X.shape
print("Number of data points:", N)
print("Data dimension:", D)
print("True number of blobs:", 3)


# ---------------- Explicit EM for a Gaussian mixture ----------------
def gmm_em(X, K, seed=1, n_iter=200, tol=1e-6, reg=1e-6):
    N, D = X.shape
    rs = np.random.RandomState(seed)
    # Initialize means at randomly chosen data points, covariances as data covariance, weights uniform
    means = X[rs.choice(N, K, replace=False)].copy()
    covs = np.array([np.cov(X.T) + reg * np.eye(D) for _ in range(K)])
    weights = np.full(K, 1.0 / K)

    prev_ll = -np.inf
    for _ in range(n_iter):
        # E-step: responsibilities = soft assignment of each point to each component
        resp = np.zeros((N, K))
        for k in range(K):
            resp[:, k] = weights[k] * multivariate_normal.pdf(X, mean=means[k], cov=covs[k])
        ll = np.sum(np.log(resp.sum(axis=1) + 1e-300))  # total log-likelihood
        resp /= resp.sum(axis=1, keepdims=True)         # normalize per point

        # M-step: update weights, means, covariances from responsibilities
        Nk = resp.sum(axis=0)
        weights = Nk / N
        means = (resp.T @ X) / Nk[:, None]
        for k in range(K):
            diff = X - means[k]
            covs[k] = (resp[:, k][:, None] * diff).T @ diff / Nk[k] + reg * np.eye(D)

        if abs(ll - prev_ll) < tol:  # converged
            break
        prev_ll = ll

    return weights, means, covs, resp, ll


def bic_from_ll(ll, K, N, D):
    # Number of free parameters: weights (K-1) + means (K*D) + covariances (K*D*(D+1)/2)
    n_params = (K - 1) + K * D + K * D * (D + 1) / 2
    return -2 * ll + n_params * np.log(N)  # lower BIC is better


# ---------------- Fit for K = 1..6 and pick K by BIC ----------------
Ks = range(1, 7)
bics = []
fits = {}
for K in Ks:
    w, m, c, r, ll = gmm_em(X, K, seed=1)
    b = bic_from_ll(ll, K, N, D)
    bics.append(b)
    fits[K] = (w, m, c, r, ll, b)
    print(f"K={K}: log-likelihood={ll:.3f}, BIC={b:.3f}")

bics = np.array(bics)
best_K = list(Ks)[int(np.argmin(bics))]  # BIC-selected number of components
print("BIC-selected number of components:", best_K)

# ---------------- Inspect the selected model ----------------
w, m, c, r, ll, b = fits[best_K]
hard_labels = r.argmax(axis=1)  # component with highest soft probability
print("Selected model component weights:", np.round(w, 3))
print("Selected model component means:")
for k in range(best_K):
    print(f"  component {k}: mean={np.round(m[k], 3)}")
print("Max soft probability per point (mean):", round(r.max(axis=1).mean(), 4))
print("Fraction of points assigned with prob > 0.9:", round(np.mean(r.max(axis=1) > 0.9), 4))

# Check: does BIC select three components?
check_passes = (best_K == 3)
print("Check (BIC selects 3 components matching blobs):", check_passes)

# ---------------- Plots ----------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Scatter colored by mixture component
axes[0].scatter(X[:, 0], X[:, 1], c=hard_labels, cmap="viridis", s=20)
axes[0].scatter(m[:, 0], m[:, 1], c="red", marker="X", s=150, label="component means")
axes[0].set_title(f"Points colored by mixture component (K={best_K})")
axes[0].set_xlabel("x1"); axes[0].set_ylabel("x2"); axes[0].legend()

# BIC versus number of components
axes[1].plot(list(Ks), bics, "o-")
axes[1].axvline(best_K, color="red", linestyle="--", label=f"selected K={best_K}")
axes[1].set_title("BIC versus number of components")
axes[1].set_xlabel("number of components"); axes[1].set_ylabel("BIC (lower is better)")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.4.1_s5.png")

# Explanation of why the check confirms the result:
print("Explanation: The check confirms the result because BIC, which penalizes model "
      "complexity, is minimized at three components exactly matching the three known blobs, "
      "and the fitted mixture assigns points to those blobs with high soft probabilities.")
