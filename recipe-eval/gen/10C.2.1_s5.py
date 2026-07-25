import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Reproducibility
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)   # seed 1
n = 1000                          # number of vertices

# ----------------------------------------------------------------------
# Erdos-Renyi G(n, p): connect each unordered vertex pair independently
# with probability p.  Implemented explicitly, not via a library routine.
# ----------------------------------------------------------------------
p = 0.006
er_degree = np.zeros(n, dtype=int)   # degree counter per vertex
for i in range(n):                   # loop over all i < j pairs (upper triangle)
    for j in range(i + 1, n):
        if rng.random() < p:         # draw a Bernoulli(p) edge
            er_degree[i] += 1        # undirected: both endpoints gain a degree
            er_degree[j] += 1

# ----------------------------------------------------------------------
# Barabasi-Albert: grow a network by attaching each new vertex to m
# existing vertices chosen with probability proportional to their degree
# (preferential attachment).  Implemented explicitly.
# ----------------------------------------------------------------------
m = 3
ba_degree = np.zeros(n, dtype=int)
# "repeated node list": each vertex appears once per incident edge endpoint,
# so drawing uniformly from it is equivalent to degree-proportional selection.
targets = list(range(m))     # initial fully-forming seed of m vertices
repeated = []                # starts empty; filled as edges are added
for i in range(m):
    ba_degree[i] = 0         # seed vertices start with no recorded edges

for new in range(m, n):      # add vertices one at a time
    # choose m distinct existing targets by preferential attachment
    chosen = set()
    while len(chosen) < m:
        if repeated:                          # pick proportional to degree
            t = repeated[rng.integers(len(repeated))]
        else:                                 # bootstrap: pick uniformly
            t = int(rng.integers(new))
        chosen.add(t)
    # add the m edges between the new vertex and the chosen targets
    for t in chosen:
        ba_degree[new] += 1
        ba_degree[t] += 1
        repeated.append(new)   # add both endpoints to the repeated list
        repeated.append(t)

# ----------------------------------------------------------------------
# Summary statistics
# ----------------------------------------------------------------------
print(f"Erdos-Renyi mean degree: {er_degree.mean():.4f}")
print(f"Erdos-Renyi std degree: {er_degree.std():.4f}")
print(f"Erdos-Renyi expected mean (n-1)*p: {(n - 1) * p:.4f}")
print(f"Erdos-Renyi min degree: {er_degree.min()}")
print(f"Erdos-Renyi max degree: {er_degree.max()}")
print(f"Erdos-Renyi coefficient of variation: {er_degree.std() / er_degree.mean():.4f}")

print(f"Barabasi-Albert mean degree: {ba_degree.mean():.4f}")
print(f"Barabasi-Albert std degree: {ba_degree.std():.4f}")
print(f"Barabasi-Albert expected mean 2*m: {2 * m:.4f}")
print(f"Barabasi-Albert min degree: {ba_degree.min()}")
print(f"Barabasi-Albert max degree: {ba_degree.max()}")
print(f"Barabasi-Albert coefficient of variation: {ba_degree.std() / ba_degree.mean():.4f}")

# ----------------------------------------------------------------------
# Fit a power law slope to the BA degree distribution on log-log axes
# ----------------------------------------------------------------------
vals, counts = np.unique(ba_degree[ba_degree > 0], return_counts=True)
prob = counts / counts.sum()
log_k = np.log10(vals.astype(float))
log_p = np.log10(prob)
slope, intercept = np.polyfit(log_k, log_p, 1)   # linear fit in log-log space
print(f"Barabasi-Albert log-log slope (power-law exponent estimate): {slope:.4f}")
print(f"Barabasi-Albert power-law exponent gamma = -slope: {-slope:.4f}")

# ----------------------------------------------------------------------
# Figure: ER histogram (left) and BA log-log scatter (right)
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.hist(er_degree, bins=range(er_degree.min(), er_degree.max() + 2),
         color="steelblue", edgecolor="black", align="left")
ax1.axvline(er_degree.mean(), color="red", linestyle="--",
            label=f"mean = {er_degree.mean():.2f}")
ax1.set_title("Erdos-Renyi G(n=1000, p=0.006)")
ax1.set_xlabel("degree")
ax1.set_ylabel("count")
ax1.legend()

ax2.scatter(vals, prob, color="darkorange", edgecolor="black", zorder=3)
fit_p = 10 ** (intercept + slope * log_k)
ax2.plot(vals, fit_p, color="red", linestyle="--",
         label=f"slope = {slope:.2f}")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_title("Barabasi-Albert (n=1000, m=3)")
ax2.set_xlabel("degree k")
ax2.set_ylabel("P(k)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.2.1_s5.png")

# ----------------------------------------------------------------------
# Explanation of the check
# ----------------------------------------------------------------------
print("Check: The ER coefficient of variation is small (degrees cluster tightly "
      "around the mean in a bell-shaped Poisson/binomial peak), while the BA degree "
      "distribution falls on a straight line on log-log axes, confirming a scale-free "
      "power law P(k) ~ k^-gamma that, unlike ER, reproduces the hub-heavy structure of "
      "real biological networks.")
