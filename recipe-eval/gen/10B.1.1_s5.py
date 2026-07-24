import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Reproducibility ----
rng = np.random.default_rng(123)

# ---- Model: three Gaussian blobs in 2D (50 points each) ----
means = [0.0, 1.5, 3.0]          # blob center coordinate (used for both x and y)
std = 0.5
pts_per_blob = 50
blobs = [rng.normal(loc=m, scale=std, size=(pts_per_blob, 2)) for m in means]
X = np.vstack(blobs)             # n x 2 point matrix
n = X.shape[0]


def assign_clusters(X, centroids):
    # For each point, compute distance to every centroid and pick the nearest.
    # dists has shape (n_points, k)
    dists = np.linalg.norm(X[:, None, :] - centroids[None, :, :], axis=2)
    return np.argmin(dists, axis=1)


def compute_wcss(X, centroids, labels):
    # Within-cluster sum of squares: sum of squared distances to assigned centroid.
    return np.sum((X - centroids[labels]) ** 2)


def kmeans_single(X, k, max_iter, rng):
    # One run of Lloyd's algorithm from a random initialization.
    # Initialize centroids as k random points drawn from the data.
    init_idx = rng.choice(X.shape[0], size=k, replace=False)
    centroids = X[init_idx].copy()

    labels = assign_clusters(X, centroids)
    for _ in range(max_iter):
        # --- Update step: move each centroid to the mean of its members ---
        new_centroids = centroids.copy()
        for j in range(k):
            members = X[labels == j]
            if len(members) > 0:
                new_centroids[j] = members.mean(axis=0)
            else:
                # Empty cluster: re-seed to a random point to avoid dead centroids.
                new_centroids[j] = X[rng.integers(X.shape[0])]

        # --- Assignment step: nearest-centroid labels ---
        new_labels = assign_clusters(X, new_centroids)

        centroids = new_centroids
        # Converged when assignments stop changing.
        if np.array_equal(new_labels, labels):
            labels = new_labels
            break
        labels = new_labels

    wcss = compute_wcss(X, centroids, labels)
    return centroids, labels, wcss


def kmeans(X, k, n_restarts, max_iter, rng):
    # Keep the solution with the lowest WCSS over several random restarts.
    best = None
    for r in range(n_restarts):
        centroids, labels, wcss = kmeans_single(X, k, max_iter, rng)
        print(f"Restart {r + 1} WCSS: {wcss:.4f}")
        if best is None or wcss < best[2]:
            best = (centroids, labels, wcss)
    return best


# ---- Run k-means: k = 3, 10 restarts, up to 100 iterations ----
k = 3
centroids, labels, best_wcss = kmeans(X, k, n_restarts=10, max_iter=100, rng=rng)

print(f"Best WCSS over restarts: {best_wcss:.4f}")
for j in range(k):
    print(f"Centroid {j}: ({centroids[j, 0]:.4f}, {centroids[j, 1]:.4f})")
for j in range(k):
    print(f"Cluster {j} size: {int(np.sum(labels == j))}")

# ---- Check: did k-means recover the three blobs? ----
# Match each recovered centroid to the nearest true blob mean (2,2) etc.
true_centers = np.array([[m, m] for m in means])
matched = []
for tc in true_centers:
    d = np.linalg.norm(centroids - tc, axis=1)
    matched.append(np.argmin(d))
    print(f"True center ({tc[0]:.1f}, {tc[1]:.1f}) -> nearest centroid {np.argmin(d)}, distance {d.min():.4f}")

recovered = (len(set(matched)) == k) and all(
    np.linalg.norm(centroids[matched[i]] - true_centers[i]) < std for i in range(k)
)
print(f"Recovered three distinct blobs (each centroid within one std of a true mean): {recovered}")
# This check confirms the result because each recovered centroid landing within one
# standard deviation of a distinct true blob mean means k-means placed exactly one
# center in the middle of each blob, and the best-of-restarts WCSS shows the guard
# against an unlucky initialization that would otherwise split one blob and merge two.

# ---- Plot: points colored by cluster, centroids marked ----
plt.figure(figsize=(7, 6))
colors = ["tab:red", "tab:green", "tab:blue"]
for j in range(k):
    cluster_pts = X[labels == j]
    plt.scatter(cluster_pts[:, 0], cluster_pts[:, 1], c=colors[j], s=25, alpha=0.6, label=f"Cluster {j}")
plt.scatter(centroids[:, 0], centroids[:, 1], c="black", marker="X", s=250,
            edgecolors="white", linewidths=1.5, label="Centroids")
plt.title("K-means clustering of three Gaussian blobs (k=3)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10B.1.1_s5.png")
