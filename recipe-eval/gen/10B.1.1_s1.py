import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def assign_clusters(X, centroids):
    # Compute squared distance from each point to each centroid, pick nearest.
    diffs = X[:, None, :] - centroids[None, :, :]      # (n, k, 2)
    dist2 = np.sum(diffs ** 2, axis=2)                 # (n, k)
    return np.argmin(dist2, axis=1)                    # nearest-centroid label per point


def update_centroids(X, labels, k):
    # Recompute each centroid as the mean of the points currently assigned to it.
    centroids = np.zeros((k, X.shape[1]))
    for j in range(k):
        members = X[labels == j]
        if len(members) > 0:
            centroids[j] = members.mean(axis=0)
        else:
            centroids[j] = X[np.random.randint(len(X))]  # reseed an empty cluster
    return centroids


def wcss(X, labels, centroids):
    # Within-cluster sum of squares: total squared distance of points to their centroid.
    return float(np.sum((X - centroids[labels]) ** 2))


def kmeans(X, k, n_restarts=10, max_iter=100):
    best_labels, best_centroids, best_wcss = None, None, np.inf
    for restart in range(n_restarts):
        # Random initialization: pick k distinct points as starting centroids.
        centroids = X[np.random.choice(len(X), k, replace=False)].copy()
        for _ in range(max_iter):
            labels = assign_clusters(X, centroids)          # assignment step
            new_centroids = update_centroids(X, labels, k)  # update step
            if np.allclose(new_centroids, centroids):       # converged
                break
            centroids = new_centroids
        labels = assign_clusters(X, centroids)
        cur_wcss = wcss(X, labels, centroids)
        # Keep the restart with the lowest WCSS (guards against bad initializations).
        if cur_wcss < best_wcss:
            best_wcss, best_labels, best_centroids = cur_wcss, labels, centroids
    return best_labels, best_centroids, best_wcss


# --- Build the model: three Gaussian blobs of 50 points each in 2D ---
np.random.seed(123)
means = [0.0, 1.5, 3.0]
std = 0.5
n_per = 50
blobs = [np.random.normal(m, std, size=(n_per, 2)) for m in means]
X = np.vstack(blobs)  # n x 2 point matrix

# --- Run k-means ---
k = 3
labels, centroids, final_wcss = kmeans(X, k=k, n_restarts=10, max_iter=100)

# --- Report numerical results ---
print("Final within-cluster sum of squares:", final_wcss)
# Sort recovered centroids for a clean comparison to the true blob means.
order = np.argsort(centroids[:, 0])
sorted_centroids = centroids[order]
for i, c in enumerate(sorted_centroids):
    print(f"Recovered centroid {i}: x={c[0]:.4f}, y={c[1]:.4f}")
true_centers = np.array([[m, m] for m in means])
for i, c in enumerate(true_centers):
    print(f"True blob center {i}: x={c[0]:.4f}, y={c[1]:.4f}")
errors = np.linalg.norm(sorted_centroids - true_centers, axis=1)
for i, e in enumerate(errors):
    print(f"Centroid {i} distance to true center: {e:.4f}")
print("Max centroid recovery error:", float(errors.max()))
counts = np.array([np.sum(labels == j) for j in range(k)])
for j in range(k):
    print(f"Cluster {j} size: {counts[j]}")

# --- Plot: points colored by cluster, centroids marked ---
plt.figure(figsize=(7, 6))
for j in range(k):
    pts = X[labels == j]
    plt.scatter(pts[:, 0], pts[:, 1], s=25, label=f"cluster {j}")
plt.scatter(centroids[:, 0], centroids[:, 1], c="black", marker="X",
            s=250, edgecolors="white", linewidths=1.5, label="centroids")
plt.title("K-means clustering of three Gaussian blobs (k=3)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10B.1.1_s1.png")

# The check confirms the result because each recovered centroid lands within a
# small distance (well under the blob separation of 1.5) of a distinct true
# blob mean, which can only happen if k-means split the data into the three
# actual blobs rather than an unlucky poor partition.
print("Check passed: all centroids within blob separation of true means?",
      bool(errors.max() < 0.75))
