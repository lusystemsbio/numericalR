import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random
from collections import Counter
import math

# ---------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------
n = 1000          # number of vertices
p = 0.006         # Erdos-Renyi edge probability
m = 3             # Barabasi-Albert edges added per new vertex
seed = 1
random.seed(seed)

# ---------------------------------------------------------------
# Erdos-Renyi G(n, p) from scratch
# For every unordered pair (i, j) with i < j, connect them
# independently with probability p. We only track degrees.
# ---------------------------------------------------------------
er_degree = [0] * n
for i in range(n):
    for j in range(i + 1, n):
        if random.random() < p:      # independent coin flip per pair
            er_degree[i] += 1
            er_degree[j] += 1

# ---------------------------------------------------------------
# Barabasi-Albert from scratch (preferential attachment)
# Start with a small clique of m0 = m fully connected vertices.
# Each new vertex attaches to m existing distinct vertices chosen
# with probability proportional to their current degree. We use a
# "target list" that repeats each vertex once per edge (endpoint
# multiset), so drawing uniformly from it gives degree-proportional
# selection.
# ---------------------------------------------------------------
m0 = m
ba_degree = [0] * n

# Initial clique among the first m0 vertices
targets = []          # multiset of endpoints for degree-proportional sampling
for i in range(m0):
    for j in range(i + 1, m0):
        ba_degree[i] += 1
        ba_degree[j] += 1
        targets.append(i)
        targets.append(j)

# Grow the network one new vertex at a time
for new in range(m0, n):
    chosen = set()
    # pick m distinct existing vertices, proportional to degree
    while len(chosen) < m:
        candidate = random.choice(targets)   # degree-proportional draw
        chosen.add(candidate)
    for t in chosen:
        ba_degree[new] += 1
        ba_degree[t] += 1
        targets.append(new)   # add both endpoints to the multiset
        targets.append(t)

# ---------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------
def mean(x):
    return sum(x) / len(x)

def std(x):
    mu = mean(x)
    return math.sqrt(sum((v - mu) ** 2 for v in x) / len(x))

er_mean, er_std = mean(er_degree), std(er_degree)
ba_mean, ba_std = mean(ba_degree), std(ba_degree)

print(f"Erdos-Renyi   mean degree: {er_mean:.4f}")
print(f"Erdos-Renyi   std  degree: {er_std:.4f}")
print(f"Erdos-Renyi   theoretical mean (n-1)*p: {(n - 1) * p:.4f}")
print(f"Erdos-Renyi   coefficient of variation (std/mean): {er_std / er_mean:.4f}")
print(f"Erdos-Renyi   min degree: {min(er_degree)}")
print(f"Erdos-Renyi   max degree: {max(er_degree)}")
print(f"Barabasi-Albert mean degree: {ba_mean:.4f}")
print(f"Barabasi-Albert std  degree: {ba_std:.4f}")
print(f"Barabasi-Albert theoretical mean 2*m: {2 * m:.4f}")
print(f"Barabasi-Albert coefficient of variation (std/mean): {ba_std / ba_mean:.4f}")
print(f"Barabasi-Albert min degree: {min(ba_degree)}")
print(f"Barabasi-Albert max degree: {max(ba_degree)}")

# ---------------------------------------------------------------
# Power-law fit for BA: estimate exponent gamma from a straight-line
# fit to log10(k) vs log10(P(k)) (least squares), which quantifies
# the scale-free behaviour we expect (theory: gamma ~ 3).
# ---------------------------------------------------------------
ba_counts = Counter(ba_degree)
total_ba = len(ba_degree)
ks = sorted(k for k in ba_counts if k > 0)
logk = [math.log10(k) for k in ks]
logp = [math.log10(ba_counts[k] / total_ba) for k in ks]

# simple least-squares slope on log-log points
lk_mean, lp_mean = mean(logk), mean(logp)
num = sum((a - lk_mean) * (b - lp_mean) for a, b in zip(logk, logp))
den = sum((a - lk_mean) ** 2 for a in logk)
slope = num / den
intercept = lp_mean - slope * lk_mean
print(f"Barabasi-Albert log-log fit slope: {slope:.4f}")
print(f"Barabasi-Albert estimated power-law exponent gamma = -slope: {-slope:.4f}")

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Erdos-Renyi: histogram (degrees cluster tightly around the mean)
axes[0].hist(er_degree, bins=range(min(er_degree), max(er_degree) + 2),
             color="steelblue", edgecolor="black", align="left")
axes[0].axvline(er_mean, color="red", linestyle="--",
                label=f"mean = {er_mean:.2f}")
axes[0].set_title("Erdos-Renyi degree distribution (histogram)")
axes[0].set_xlabel("degree k")
axes[0].set_ylabel("count")
axes[0].legend()

# Barabasi-Albert: log-log scatter (straight line => scale-free power law)
axes[1].scatter(ks, [ba_counts[k] / total_ba for k in ks],
                color="darkorange", edgecolor="black")
fit_x = [min(ks), max(ks)]
fit_y = [10 ** (intercept + slope * math.log10(x)) for x in fit_x]
axes[1].plot(fit_x, fit_y, "r--",
             label=f"fit slope = {slope:.2f}")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_title("Barabasi-Albert degree distribution (log-log)")
axes[1].set_xlabel("degree k")
axes[1].set_ylabel("P(k)")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.2.1_s3.png")

# ---------------------------------------------------------------
# Explanation of why the check confirms the result
# ---------------------------------------------------------------
print("Check explanation: a tightly clustered ER histogram (small std/mean) reflects "
      "its Poisson-like homogeneity, while a straight BA line on log-log axes means "
      "P(k) ~ k^(-gamma), the hallmark of a scale-free power-law with hubs like real biological networks.")
