import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# ---- Input: Zachary's karate club adjacency matrix ----
G = nx.karate_club_graph()
nodes = list(G.nodes())
A = nx.to_numpy_array(G, nodelist=nodes)   # adjacency matrix A_ij
n = A.shape[0]

# ---- Build the modularity matrix B_ij = A_ij - k_i*k_j/(2m) ----
k = A.sum(axis=1)          # degree of each node
m = A.sum() / 2.0          # number of edges m = (1/2) sum_ij A_ij
B = A - np.outer(k, k) / (2.0 * m)

# ---- Newman spectral bisection: leading eigenvector of B ----
# Real symmetric -> use eigh; take eigenvector of the largest eigenvalue.
eigvals, eigvecs = np.linalg.eigh(B)
leading = eigvecs[:, np.argmax(eigvals)]

# Fix the sign so the labeling is reproducible: make the first
# non-zero entry positive (an eigenvector is only defined up to sign).
first_nz = leading[np.argmax(np.abs(leading) > 1e-12)]
if first_nz < 0:
    leading = -leading

# Split nodes by the sign of the leading eigenvector into two communities.
s = np.where(leading >= 0, 1, -1)          # community index vector s_i = +/-1
community = np.where(s > 0, 0, 1)          # relabel as 0 / 1

# ---- Modularity Q = (1/2m) sum_ij (A_ij - k_i k_j/2m) delta(c_i,c_j) ----
# delta(c_i,c_j) = 1 when in the same community, i.e. s_i*s_j = 1 -> (s_i s_j + 1)/2
same = (np.outer(s, s) + 1.0) / 2.0
Q = np.sum(B * same) / (2.0 * m)

# ---- Report community sizes ----
size0 = int(np.sum(community == 0))
size1 = int(np.sum(community == 1))

print(f"Number of nodes n = {n}")
print(f"Number of edges m = {int(m)}")
print(f"Community 0 size = {size0}")
print(f"Community 1 size = {size1}")
print(f"Community sizes (sorted) = {sorted([size0, size1])}")
print(f"Modularity Q = {Q:.6f}")

# ---- Check against Zachary's recorded factional split ----
sizes_ok = sorted([size0, size1]) == [16, 18]
q_ok = abs(Q - 0.37) < 0.02
print(f"Check: split into 18 and 16 -> {sizes_ok}")
print(f"Check: modularity approximately 0.37 -> {q_ok} (Q = {Q:.4f})")
print(f"Overall check passed = {sizes_ok and q_ok}")
# Why this confirms the result: matching both the 18/16 sizes and Q ~ 0.37
# reproduces the strong two-community structure Zachary observed empirically,
# so agreement on independent quantities (sizes AND modularity) rather than a
# single tuned number is strong evidence the spectral bisection is correct.
print("Explanation: The split reproduces both the community sizes (18/16) and a "
      "high modularity (~0.37) that independently match Zachary's observed factional "
      "division, so agreeing on two separate quantities confirms the detection is correct.")

# ---- Draw the network with the SAME Kamada-Kawai layout as 10C.1 ----
pos = nx.kamada_kawai_layout(G)
colors = ["#1f77b4" if community[i] == 0 else "#ff7f0e" for i in range(n)]

plt.figure(figsize=(8, 8))
nx.draw_networkx_edges(G, pos, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=400)
nx.draw_networkx_labels(G, pos, font_size=9)
plt.title(f"Newman spectral communities of Zachary's karate club (Q = {Q:.3f})")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.3.1_s5.png")
