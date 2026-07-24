import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------
rng = np.random.default_rng(123)

# ---------------------------------------------------------------
# Model: three Gaussian blobs in 2D, 50 points each,
# means at 0, 1.5, 3 (in both coords), std = 0.5
# ---------------------------------------------------------------
means = [0.0, 1.5, 3.0]
std = 0.5
n_per = 50
blobs = [rng.normal(loc=m, scale=std, size=(n_per, 2)) for m in means]
X = np.vstack(blobs)  # n x 2 point matrix
n = X.shape[0]


# ---------------------------------------------------------------
# k-means (Lloyd's algorithm) implemented from scratch
# ---------------------------------------------------------------
def assign_clusters(X, centroids):
    # squared Euclidean distance from every point to every centroid
    d2 = ((X[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    labels = d2.argmin(axis=1)          # nearest-centroid assignment
    wcss = d2[np.arange(len(X)), labels].sum()  # within-cluster sum of squares
    return labels, wcss


def update_centroids(X, labels, k, prev):
    new = prev.copy()
    for j in range(k):
        pts = X[labels == j]
        if len(pts) > 0:                # mean of assigned points
            new[j] = pts.mean(axis=0)
        # empty cluster: keep previous centroid
    return new


def kmeans_once(X, k, rng, max_iter=100):
    # random initial centroids: pick k distinct points
    idx = rng.choice(len(X), size=k, replace=False)
    centroids = X[idx].astype(float)
    labels, wcss = assign_clusters(X, centroids)
    for _ in range(max_iter):
        centroids = update_centroids(X, labels, k, centroids)   # mean update
        new_labels, wcss = assign_clusters(X, centroids)        # reassignment
        if np.array_equal(new_labels, labels):  # converged
            labels = new_labels
            break
        labels = new_labels
    return centroids, labels, wcss


def kmeans(X, k, rng, n_restarts=10, max_iter=100):
    best = None  # keep the run with the lowest WCSS across restarts
    for _ in range(n_restarts):
        centroids, labels, wcss = kmeans_once(X, k, rng, max_iter)
        if best is None or wcss < best[2]:
            best = (centroids, labels, wcss)
    return best


# ---------------------------------------------------------------
# Run: k = 3, 10 restarts, up to 100 iterations
# ---------------------------------------------------------------
k = 3
centroids, labels, wcss = kmeans(X, k, rng, n_restarts=10, max_iter=100)

print(f"Best within-cluster sum of squares (WCSS): {wcss:.6f}")
for j in range(k):
    print(f"Centroid {j}: ({centroids[j, 0]:.4f}, {centroids[j, 1]:.4f})  size = {int((labels == j).sum())}")

# ---------------------------------------------------------------
# Check: does k-means recover the three blobs?
# Match each recovered centroid to its nearest true blob mean.
# ---------------------------------------------------------------
true_means = np.array([[m, m] for m in means])  # (0,0), (1.5,1.5), (3,3)
order = np.argsort(centroids[:, 0])             # sort recovered centroids by x
sorted_centroids = centroids[order]
print("\nRecovery check (recovered centroid vs. true blob mean):")
max_err = 0.0
for tc, tm in zip(sorted_centroids, true_means):
    err = np.sqrt(((tc - tm) ** 2).sum())
    max_err = max(max_err, err)
    print(f"true ({tm[0]:.2f}, {tm[1]:.2f}) -> recovered ({tc[0]:.4f}, {tc[1]:.4f}), distance = {err:.4f}")
print(f"Maximum centroid-to-true-mean distance: {max_err:.4f}")
recovered = max_err < std  # each centroid within one std of a blob center
print(f"Recovered all three blob centers (within one std = {std}): {recovered}")

# Compare against a deliberately unlucky single run to show restarts help
worst_wcss = -np.inf
bad_rng = np.random.default_rng(999)
for _ in range(50):
    _, _, w = kmeans_once(X, k, bad_rng, max_iter=100)
    worst_wcss = max(worst_wcss, w)
print(f"\nWorst single-run WCSS observed (no restart guard): {worst_wcss:.6f}")
print(f"Restart-selected WCSS is <= worst single run: {wcss <= worst_wcss}")

# ---------------------------------------------------------------
# Explanation of why the check confirms the result
# ---------------------------------------------------------------
print(
    "\nWhy this confirms it: because each of the three recovered centroids lands within one "
    "standard deviation of a distinct true blob mean (and the kept solution has the lowest WCSS "
    "across restarts), we know k-means placed one centroid in the middle of each blob rather than "
    "collapsing to a poor split from an unlucky initialization."
)

# ---------------------------------------------------------------
# Plot: points colored by cluster, centroids marked
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
colors = ["tab:red", "tab:green", "tab:blue"]
for j in range(k):
    pts = X[labels == j]
    ax.scatter(pts[:, 0], pts[:, 1], c=colors[j], s=25, alpha=0.7, label=f"cluster {j}")
ax.scatter(centroids[:, 0], centroids[:, 1], c="black", marker="X", s=250,
           edgecolors="white", linewidths=1.5, label="centroids", zorder=5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("k-means (from scratch): 3 Gaussian blobs, k=3")
ax.legend()
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10B.1.1_s2.png")
print("\nSaved figure to 10B.1.1_s2.png")
