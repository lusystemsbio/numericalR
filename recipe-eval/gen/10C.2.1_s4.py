import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Reproducibility: single RNG seeded once, shared by both models.
# ---------------------------------------------------------------
rng = np.random.default_rng(1)

# Parameters
n = 1000        # number of vertices
p = 0.006       # Erdos-Renyi edge probability
m = 3           # Barabasi-Albert edges added per new vertex


# ---------------------------------------------------------------
# Erdos-Renyi G(n, p), implemented from scratch.
# Every unordered pair (i, j) is connected independently with prob p.
# ---------------------------------------------------------------
def erdos_renyi(n, p, rng):
    degrees = np.zeros(n, dtype=int)   # track degree of each vertex
    # loop over all unique vertex pairs i < j
    for i in range(n):
        # draw a Bernoulli(p) outcome for every j > i at once
        connect = rng.random(n - i - 1) < p
        for offset, is_edge in enumerate(connect):
            if is_edge:
                j = i + 1 + offset
                degrees[i] += 1        # edge contributes to both endpoints
                degrees[j] += 1
    return degrees


# ---------------------------------------------------------------
# Barabasi-Albert, implemented from scratch (preferential attachment).
# Start from a small connected seed of m+1 vertices, then add each new
# vertex with m edges chosen proportional to existing vertex degrees.
# ---------------------------------------------------------------
def barabasi_albert(n, m, rng):
    degrees = np.zeros(n, dtype=int)
    # "targets" is a list of vertex ids; each id appears once per incident
    # edge, so drawing uniformly from it gives probability proportional
    # to degree (classic preferential-attachment trick).
    targets = []

    # Seed: a small clique of the first m+1 vertices so every vertex
    # starts with nonzero degree and can be chosen.
    seed = m + 1
    for i in range(seed):
        for j in range(i + 1, seed):
            degrees[i] += 1
            degrees[j] += 1
            targets.append(i)
            targets.append(j)

    # Grow the network one new vertex at a time.
    for new in range(seed, n):
        chosen = set()
        # keep sampling until we have m distinct targets for this vertex
        while len(chosen) < m:
            t = targets[rng.integers(len(targets))]  # prob proportional to degree
            chosen.add(t)
        for t in chosen:
            degrees[new] += 1
            degrees[t] += 1
            targets.append(new)   # record both endpoints of the new edge
            targets.append(t)
    return degrees


# ---------------------------------------------------------------
# Generate the two networks.
# ---------------------------------------------------------------
er_deg = erdos_renyi(n, p, rng)
ba_deg = barabasi_albert(n, m, rng)

# ---------------------------------------------------------------
# Summary statistics.
# ---------------------------------------------------------------
er_mean = er_deg.mean()
er_std = er_deg.std()
ba_mean = ba_deg.mean()
ba_std = ba_deg.std()

print(f"Erdos-Renyi mean degree:            {er_mean:.4f}")
print(f"Erdos-Renyi std degree:             {er_std:.4f}")
print(f"Erdos-Renyi expected mean (n-1)*p:  {(n - 1) * p:.4f}")
print(f"Erdos-Renyi coeff. of variation:    {er_std / er_mean:.4f}")
print(f"Barabasi-Albert mean degree:        {ba_mean:.4f}")
print(f"Barabasi-Albert std degree:         {ba_std:.4f}")
print(f"Barabasi-Albert expected mean (2m):  {2 * m:.4f}")
print(f"Barabasi-Albert coeff. of variation: {ba_std / ba_mean:.4f}")
print(f"Barabasi-Albert max degree (hub):   {ba_deg.max()}")
print(f"Erdos-Renyi max degree:             {er_deg.max()}")

# ---------------------------------------------------------------
# Quantitative check on the BA log-log power law:
# fit a straight line to log(k) vs log(count) and report the slope,
# which estimates the (negative) power-law exponent.
# ---------------------------------------------------------------
ba_k, ba_counts = np.unique(ba_deg, return_counts=True)
positive = ba_counts > 0
logk = np.log10(ba_k[positive])
logc = np.log10(ba_counts[positive])
slope, intercept = np.polyfit(logk, logc, 1)
print(f"Barabasi-Albert log-log fit slope (~ -gamma): {slope:.4f}")

# ---------------------------------------------------------------
# Plots: ER histogram (left), BA log-log scatter (right).
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(er_deg, bins=range(er_deg.min(), er_deg.max() + 2),
             color="steelblue", edgecolor="black", align="left")
axes[0].axvline(er_mean, color="red", linestyle="--", label=f"mean={er_mean:.1f}")
axes[0].set_title("Erdos-Renyi degree distribution")
axes[0].set_xlabel("degree k")
axes[0].set_ylabel("number of vertices")
axes[0].legend()

axes[1].scatter(ba_k[positive], ba_counts[positive], color="darkorange",
                edgecolor="black", zorder=3)
axes[1].plot(10**logk, 10**(intercept + slope * logk), "r--",
             label=f"slope={slope:.2f}")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_title("Barabasi-Albert degree distribution (log-log)")
axes[1].set_xlabel("degree k")
axes[1].set_ylabel("count")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.2.1_s4.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("Explanation: The ER degrees have a small coefficient of variation and "
      "cluster tightly in a single-peaked histogram around (n-1)*p, whereas the "
      "BA degrees fall on a straight line in log-log axes, and a straight log-log "
      "line means count ~ k^slope, i.e. a scale-free power law with heavy-tailed "
      "hubs like those seen in biological networks.")
