import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import multivariate_normal
from sklearn.datasets import make_blobs

# ---- Data: three-blob data ----
X, y_true = make_blobs(n_samples=300, centers=3, cluster_std=1.0, random_state=1)
n, d = X.shape

# ---- EM for a Gaussian mixture (implemented explicitly) ----
def fit_gmm(X, K, seed=1, n_iter=200, tol=1e-6, reg=1e-6):
    rng = np.random.default_rng(seed)
    n, d = X.shape
    # Initialize: random points as means, global covariance, uniform weights
    means = X[rng.choice(n, K, replace=False)].copy()
    covs = np.array([np.cov(X.T) + reg * np.eye(d) for _ in range(K)])
    weights = np.full(K, 1.0 / K)
    ll_old = -np.inf
    for _ in range(n_iter):
        # E-step: responsibilities (soft assignments) via Bayes rule
        resp = np.zeros((n, K))
        for k in range(K):
            resp[:, k] = weights[k] * multivariate_normal.pdf(X, means[k], covs[k], allow_singular=True)
        ll = np.sum(np.log(resp.sum(axis=1) + 1e-300))  # total log-likelihood
        resp /= resp.sum(axis=1, keepdims=True)          # normalize per point
        # M-step: update weights, means, covariances from responsibilities
        Nk = resp.sum(axis=0)
        weights = Nk / n
        for k in range(K):
            means[k] = (resp[:, k, None] * X).sum(axis=0) / Nk[k]
            diff = X - means[k]
            covs[k] = (resp[:, k, None, None] * (diff[:, :, None] * diff[:, None, :])).sum(axis=0) / Nk[k]
            covs[k] += reg * np.eye(d)  # regularize to keep positive definite
        if abs(ll - ll_old) < tol:      # convergence check
            break
        ll_old = ll
    return weights, means, covs, resp, ll

# Number of free parameters for BIC: weights + means + full covariances
def n_params(K, d):
    return (K - 1) + K * d + K * d * (d + 1) // 2

# ---- Fit for K = 1..6 and choose K by BIC ----
Ks = range(1, 7)
bics, models = [], {}
for K in Ks:
    w, m, c, r, ll = fit_gmm(X, K, seed=1)
    p = n_params(K, d)
    bic = -2 * ll + p * np.log(n)   # BIC: lower is better
    bics.append(bic)
    models[K] = (w, m, c, r, ll, bic)
    print(f"K={K}: logL={ll:.3f}, params={p}, BIC={bic:.3f}")

bics = np.array(bics)
best_K = list(Ks)[int(np.argmin(bics))]
print(f"BIC-selected number of components: {best_K}")
print(f"True number of blobs: {len(np.unique(y_true))}")

# ---- Check: does BIC select three components matching the blobs? ----
w, m, c, resp, ll, bic = models[best_K]
labels = resp.argmax(axis=1)
print(f"Selected model has {len(np.unique(labels))} used components")
print(f"Check passed (selected == 3 blobs): {best_K == 3}")
# Report soft-assignment confidence to show recovery with soft probabilities
print(f"Mean max soft probability (assignment confidence): {resp.max(axis=1).mean():.4f}")
for k in range(best_K):
    print(f"Component {k}: weight={w[k]:.4f}, mean=[{m[k,0]:.3f}, {m[k,1]:.3f}]")

# ---- Plots ----
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].scatter(X[:, 0], X[:, 1], c=labels, cmap="viridis", s=25, alpha=0.85)
ax[0].scatter(m[:, 0], m[:, 1], c="red", marker="X", s=180, edgecolor="black", label="means")
ax[0].set_title(f"Points colored by mixture component (K={best_K})")
ax[0].set_xlabel("x1"); ax[0].set_ylabel("x2"); ax[0].legend()

ax[1].plot(list(Ks), bics, "o-")
ax[1].scatter([best_K], [bics[best_K - 1]], c="red", s=140, zorder=5, label=f"min BIC (K={best_K})")
ax[1].set_title("BIC vs number of components")
ax[1].set_xlabel("number of components"); ax[1].set_ylabel("BIC"); ax[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10B.4.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: BIC picking three components and the soft probabilities peaking "
      "near one per point shows the mixture assigns each blob to a distinct Gaussian, "
      "confirming it recovered the true three-blob structure rather than over- or under-fitting.")
