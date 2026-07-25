import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import random
from collections import Counter

# ---- reproducibility ----
random.seed(1)
np.random.seed(1)

n = 1000  # number of vertices

# =========================================================
# Erdos-Renyi G(n, p): connect each unordered vertex pair
# independently with probability p, built from scratch.
# =========================================================
def erdos_renyi(n, p):
    degrees = [0] * n
    # iterate over every unordered pair (i, j), i < j
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < p:   # independent coin flip per pair
                degrees[i] += 1
                degrees[j] += 1
    return degrees

# =========================================================
# Barabasi-Albert: grow the graph by adding one vertex at a
# time, attaching it to m existing vertices with probability
# proportional to their current degree (preferential
# attachment), built from scratch.
# =========================================================
def barabasi_albert(n, m):
    # 'repeated_nodes' holds each node once per edge endpoint;
    # sampling from it uniformly gives probability proportional to degree.
    repeated_nodes = []
    degrees = [0] * n

    # seed: start with m fully-connected initial vertices
    for i in range(m):
        for j in range(i + 1, m):
            degrees[i] += 1
            degrees[j] += 1
            repeated_nodes.extend([i, j])
    # if the seed produced no edges yet, seed the pool with the initial nodes
    if not repeated_nodes:
        repeated_nodes = list(range(m))

    # add remaining vertices one by one
    for new in range(m, n):
        targets = set()
        # pick m distinct targets, favoring high-degree nodes
        while len(targets) < m:
            targets.add(random.choice(repeated_nodes))
        # connect new vertex to each chosen target
        for t in targets:
            degrees[new] += 1
            degrees[t] += 1
            repeated_nodes.extend([new, t])  # update pool so future picks stay proportional
    return degrees

# ---- generate the two networks ----
er_degrees = erdos_renyi(n, 0.006)
ba_degrees = barabasi_albert(n, 3)

# ---- summary statistics ----
er_mean = np.mean(er_degrees)
er_std = np.std(er_degrees)
ba_mean = np.mean(ba_degrees)
ba_std = np.std(ba_degrees)

print(f"Erdos-Renyi (n=1000, p=0.006) mean degree: {er_mean:.4f}")
print(f"Erdos-Renyi (n=1000, p=0.006) std degree: {er_std:.4f}")
print(f"Erdos-Renyi expected mean degree (n-1)*p: {(n - 1) * 0.006:.4f}")
print(f"Erdos-Renyi coefficient of variation (std/mean): {er_std / er_mean:.4f}")
print(f"Barabasi-Albert (n=1000, m=3) mean degree: {ba_mean:.4f}")
print(f"Barabasi-Albert (n=1000, m=3) std degree: {ba_std:.4f}")
print(f"Barabasi-Albert expected mean degree ~2m: {2 * 3:.4f}")
print(f"Barabasi-Albert coefficient of variation (std/mean): {ba_std / ba_mean:.4f}")
print(f"Erdos-Renyi max degree: {max(er_degrees)}")
print(f"Barabasi-Albert max degree: {max(ba_degrees)}")

# ---- degree distribution for the log-log power-law check ----
ba_counts = Counter(ba_degrees)
ba_k = np.array(sorted(k for k in ba_counts if k > 0))
ba_pk = np.array([ba_counts[k] for k in ba_k], dtype=float) / n

# fit a straight line to log10(P(k)) vs log10(k): slope = power-law exponent
log_k = np.log10(ba_k)
log_pk = np.log10(ba_pk)
slope, intercept = np.polyfit(log_k, log_pk, 1)
print(f"Barabasi-Albert log-log fitted slope (power-law exponent -gamma): {slope:.4f}")
print(f"Barabasi-Albert log-log fitted intercept: {intercept:.4f}")

# ---- plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Erdos-Renyi: histogram (degrees cluster tightly around the mean)
ax1.hist(er_degrees, bins=range(min(er_degrees), max(er_degrees) + 2),
         color="steelblue", edgecolor="black", align="left")
ax1.axvline(er_mean, color="red", linestyle="--", label=f"mean = {er_mean:.2f}")
ax1.set_title("Erdos-Renyi degree distribution (histogram)")
ax1.set_xlabel("degree k")
ax1.set_ylabel("count")
ax1.legend()

# Barabasi-Albert: log-log scatter (straight line = scale-free power law)
ax2.scatter(ba_k, ba_pk, color="darkorange", edgecolor="black", label="P(k)")
ax2.plot(ba_k, 10 ** (intercept + slope * log_k), color="red",
         linestyle="--", label=f"fit slope = {slope:.2f}")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_title("Barabasi-Albert degree distribution (log-log)")
ax2.set_xlabel("degree k")
ax2.set_ylabel("P(k)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.2.1_s1.png")

# ---- one-sentence explanation of why the check confirms the result ----
print("Explanation: A small coefficient of variation with a bell-shaped histogram "
      "confirms the Erdos-Renyi degrees concentrate near the mean, while a straight "
      "line on log-log axes confirms P(k) ~ k^(-gamma), the scale-free power law that "
      "produces hubs seen in real biological networks.")
