import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random
from collections import Counter
import math

# ---------------- Parameters ----------------
n = 1000            # number of vertices
p = 0.006           # Erdos-Renyi edge probability
m = 3               # Barabasi-Albert attachments per new vertex
seed = 1
random.seed(seed)   # reproducibility

# ---------------- Erdos-Renyi G(n, p) from scratch ----------------
# Consider every unordered pair (i, j); connect them independently with prob p.
er_degrees = [0] * n
for i in range(n):
    for j in range(i + 1, n):
        if random.random() < p:   # independent Bernoulli(p) trial per pair
            er_degrees[i] += 1
            er_degrees[j] += 1

# ---------------- Barabasi-Albert (preferential attachment) from scratch ----------------
# Start with a small seed clique of m vertices, then add vertices one at a time.
# Each new vertex connects to m existing vertices chosen with probability
# proportional to their current degree. We implement this with a "repeated node"
# list: each vertex appears once per edge endpoint, so uniform sampling from the
# list gives probability proportional to degree.
ba_degrees = [0] * n
repeated_nodes = []   # multiset of endpoints; frequency == degree

# initial seed: connect the first m vertices in a small clique
for i in range(m):
    for j in range(i + 1, m):
        ba_degrees[i] += 1
        ba_degrees[j] += 1
        repeated_nodes.extend([i, j])

# add the remaining vertices
for new in range(m, n):
    targets = set()
    # pick m distinct targets, each proportional to current degree
    while len(targets) < m:
        if repeated_nodes:
            candidate = random.choice(repeated_nodes)  # prob ~ degree
        else:
            candidate = random.randrange(new)          # fallback: uniform
        if candidate != new:
            targets.add(candidate)
    # wire up the new vertex to the chosen targets
    for t in targets:
        ba_degrees[new] += 1
        ba_degrees[t] += 1
        repeated_nodes.extend([new, t])

# ---------------- Summary statistics ----------------
er_mean = sum(er_degrees) / n
er_var = sum((d - er_mean) ** 2 for d in er_degrees) / n
er_std = math.sqrt(er_var)

ba_mean = sum(ba_degrees) / n
ba_var = sum((d - ba_mean) ** 2 for d in ba_degrees) / n
ba_std = math.sqrt(ba_var)

print(f"Erdos-Renyi  n={n}, p={p}")
print(f"ER expected mean degree (n-1)*p: {(n - 1) * p:.4f}")
print(f"ER observed mean degree: {er_mean:.4f}")
print(f"ER degree std dev: {er_std:.4f}")
print(f"ER coefficient of variation (std/mean): {er_std / er_mean:.4f}")
print(f"ER min degree: {min(er_degrees)}")
print(f"ER max degree: {max(er_degrees)}")

print(f"Barabasi-Albert  n={n}, m={m}")
print(f"BA expected mean degree ~2*m: {2 * m:.4f}")
print(f"BA observed mean degree: {ba_mean:.4f}")
print(f"BA degree std dev: {ba_std:.4f}")
print(f"BA coefficient of variation (std/mean): {ba_std / ba_mean:.4f}")
print(f"BA min degree: {min(ba_degrees)}")
print(f"BA max degree: {max(ba_degrees)}")

# ---------------- Check: fit a power law to the BA log-log tail ----------------
# Build the degree distribution P(k) and fit a line in log-log space.
ba_counts = Counter(ba_degrees)
ks = sorted(k for k in ba_counts if k > 0)
logk = [math.log10(k) for k in ks]
logp = [math.log10(ba_counts[k] / n) for k in ks]

# simple least-squares line: log P(k) = a + b*log k  (b is the negative exponent)
N = len(logk)
sx = sum(logk); sy = sum(logp)
sxx = sum(x * x for x in logk); sxy = sum(x * y for x, y in zip(logk, logp))
slope = (N * sxy - sx * sy) / (N * sxx - sx * sx)
intercept = (sy - slope * sx) / N
# coefficient of determination R^2 of the log-log fit
ybar = sy / N
ss_tot = sum((y - ybar) ** 2 for y in logp)
ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(logk, logp))
r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

print(f"BA log-log fit slope: {slope:.4f}")
print(f"BA implied power-law exponent gamma (= -slope): {-slope:.4f}")
print(f"BA log-log fit R^2: {r2:.4f}")

# Explanation of the check:
print("Explanation: A tight, symmetric ER histogram near its mean signals a "
      "Poisson-like distribution with a characteristic scale, whereas a "
      "straight BA line on log-log axes (high R^2) means P(k) ~ k^-gamma, "
      "confirming the scale-free power law that better matches biological networks.")

# ---------------- Plots ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# ER: histogram, expected to cluster tightly around the mean
ax1.hist(er_degrees, bins=range(min(er_degrees), max(er_degrees) + 2),
         color="steelblue", edgecolor="black", align="left")
ax1.axvline(er_mean, color="red", linestyle="--", label=f"mean={er_mean:.2f}")
ax1.set_title("Erdos-Renyi degree distribution (histogram)")
ax1.set_xlabel("degree k")
ax1.set_ylabel("count")
ax1.legend()

# BA: log-log scatter of P(k), expected to follow a straight line
ax2.scatter(ks, [ba_counts[k] / n for k in ks], color="darkorange",
            edgecolor="black", zorder=3, label="P(k)")
fit_p = [10 ** (intercept + slope * x) for x in logk]
ax2.plot(ks, fit_p, "r--", label=f"fit slope={slope:.2f}")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_title("Barabasi-Albert degree distribution (log-log)")
ax2.set_xlabel("degree k")
ax2.set_ylabel("P(k)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.2.1_s2.png")
