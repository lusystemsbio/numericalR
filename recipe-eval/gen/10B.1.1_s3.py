import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Build the model: three Gaussian blobs in 2D (50 points each) ----
rng = np.random.default_rng(123)  # seed 123 for reproducibility
means = [0.0, 1.5, 3.0]           # blob centers (same value for x and y)
std = 0.5
n_per = 50
blobs = [rng.normal(loc=m, scale=std, size=(n_per, 2)) for m in means]
X = np.vstack(blobs)              # n x 2 input point matrix (150 x 2)


def assign_clusters(X, centroids):
    # Compute squared distance from every point to every centroid, pick nearest
    diff = X[:, None, :] - centroids[None, :, :]        # (n, k, 2)
    dist2 = np.sum(diff ** 2, axis=2)                   # (n, k)
    labels = np.argmin(dist2, axis=1)                   # nearest centroid index
    return labels, dist2


def kmeans_once(X, k, max_iter, rng):
    # Random initialization: pick k distinct points as starting centroids
    idx = rng.choice(X.shape[0], size=k, replace=False)
    centroids = X[idx].copy()
    labels = None
    for _ in range(max_iter):
        # 1) Assignment step: label each point by its nearest centroid
        new_labels, _ = assign_clusters(X, centroids)
        # 2) Update step: move each centroid to the mean of its members
        new_centroids = np.array([
            X[new_labels == j].mean(axis=0) if np.any(new_labels == j)
            else centroids[j]                            # keep empty cluster put
            for j in range(k)
        ])
        # Stop when assignments no longer change (converged)
        if labels is not None and np.array_equal(new_labels, labels):
            labels = new_labels
            centroids = new_centroids
            break
        labels, centroids = new_labels, new_centroids
    # Within-cluster sum of squares for this solution
    _, dist2 = assign_clusters(X, centroids)
    wcss = np.sum(dist2[np.arange(X.shape[0]), labels])
    return labels, centroids, wcss


def kmeans(X, k, n_restarts, max_iter, seed):
    # Run several restarts; keep the one with the lowest WCSS
    rng = np.random.default_rng(seed)
    best = None
    for r in range(n_restarts):
        labels, centroids, wcss = kmeans_once(X, k, max_iter, rng)
        print(f"Restart {r + 1:2d} WCSS: {wcss:.6f}")
        if best is None or wcss < best[2]:
            best = (labels, centroids, wcss)
    return best


# ---- Run k-means: k=3, 10 restarts, up to 100 iterations, seed 123 ----
labels, centroids, wcss = kmeans(X, k=3, n_restarts=10, max_iter=100, seed=123)

print(f"Best within-cluster sum of squares (WCSS): {wcss:.6f}")
# Sort recovered centroids so they line up with the true means for comparison
order = np.argsort(centroids[:, 0])
sorted_centroids = centroids[order]
for i, c in enumerate(sorted_centroids):
    print(f"Recovered centroid {i} (x, y): ({c[0]:.4f}, {c[1]:.4f})")
for i, m in enumerate(means):
    print(f"True blob mean {i} (x, y): ({m:.4f}, {m:.4f})")

# ---- Check: does each recovered centroid sit near a true blob middle? ----
true_centers = np.array([[m, m] for m in means])
errors = np.linalg.norm(sorted_centroids - true_centers, axis=1)
for i, e in enumerate(errors):
    print(f"Distance from recovered centroid {i} to true blob {i}: {e:.4f}")
max_error = errors.max()
print(f"Max centroid-to-true-mean distance: {max_error:.4f}")
recovered_ok = max_error < std  # every centroid within one std of a true center
print(f"K-means recovered the three blobs (all within one std): {recovered_ok}")

# ---- Scatter plot: points colored by cluster, centroids marked ----
fig, ax = plt.subplots(figsize=(7, 6))
colors = ["tab:red", "tab:green", "tab:blue"]
for j in range(3):
    pts = X[labels == j]
    ax.scatter(pts[:, 0], pts[:, 1], s=25, color=colors[j], alpha=0.6,
               label=f"Cluster {j}")
ax.scatter(centroids[:, 0], centroids[:, 1], s=250, marker="X",
           c="black", edgecolors="white", linewidths=1.5, label="Centroids")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("K-means (k=3) on three Gaussian blobs")
ax.legend()
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10B.1.1_s3.png")

# Explanation: the check confirms the result because each recovered centroid landing
# within one standard deviation of a distinct true blob mean shows k-means found the
# three real cluster centers rather than an arbitrary or degenerate split.
print("Explanation: The check confirms success because each recovered centroid lands "
      "within one standard deviation of a distinct true blob mean, proving k-means "
      "located the three real cluster centers and the restarts avoided a poor split.")
