import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# ---- Input: Zachary's karate club adjacency matrix ----
G = nx.karate_club_graph()
nodes = list(G.nodes())
A = nx.to_numpy_array(G, nodelist=nodes)  # adjacency matrix A_ij
n = A.shape[0]

# ---- Build the modularity matrix B ----
# degrees k_i and total edge endpoints 2m
k = A.sum(axis=1)          # degree of each node
two_m = k.sum()            # 2m = sum of degrees
m = two_m / 2.0            # number of edges
# B_ij = A_ij - k_i*k_j/2m  (the modularity matrix)
B = A - np.outer(k, k) / two_m

# ---- Newman spectral bisection: leading eigenvector of B ----
eigvals, eigvecs = np.linalg.eigh(B)     # B is symmetric
leading = eigvecs[:, np.argmax(eigvals)] # eigenvector of largest eigenvalue

# Fix the sign so the labeling is reproducible:
# force the first nonzero (largest-magnitude) component to be positive.
pivot = leading[np.argmax(np.abs(leading))]
if pivot < 0:
    leading = -leading

# Split by the sign of the leading eigenvector components.
# s_i = +1 or -1 assigns node i to one of two communities.
s = np.where(leading >= 0, 1, -1)
community = np.where(s == 1, 0, 1)   # community labels 0 / 1

# ---- Compute modularity Q for this split ----
# Q = (1/2m) * sum_ij (A_ij - k_i k_j/2m) * delta(c_i,c_j)
delta = (community[:, None] == community[None, :]).astype(float)
Q = np.sum(B * delta) / two_m

# Community sizes
size0 = int(np.sum(community == 0))
size1 = int(np.sum(community == 1))

# ---- Report numerical results ----
print(f"Number of nodes: {n}")
print(f"Number of edges m: {m:.0f}")
print(f"Community 0 size: {size0}")
print(f"Community 1 size: {size1}")
print(f"Modularity Q: {Q:.6f}")

# ---- Separate check ----
sizes = sorted([size0, size1], reverse=True)
check_sizes = (sizes == [18, 16])
check_Q = abs(Q - 0.37) < 0.02
print(f"Check split sizes are 18 and 16: {check_sizes}")
print(f"Check modularity ~ 0.37: {check_Q}")
print(f"Overall check passed: {check_sizes and check_Q}")
# This check confirms the result because reproducing Zachary's real-world 18/16
# factional split with a high modularity (~0.37) shows the spectral partition
# recovers the true community structure rather than an arbitrary cut.

# ---- Draw the network with the same Kamada-Kawai layout as 10C.1 ----
pos = nx.kamada_kawai_layout(G)
colors = ["tab:blue" if c == 0 else "tab:orange" for c in community]

plt.figure(figsize=(8, 8))
nx.draw_networkx_edges(G, pos, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=400)
nx.draw_networkx_labels(G, pos, font_size=9)
plt.title(f"Newman spectral communities (Q = {Q:.3f}, sizes {size0}/{size1})")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.3.1_s4.png")
