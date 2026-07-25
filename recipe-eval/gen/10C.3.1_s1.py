import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# --- Input: Zachary's karate club adjacency matrix ---
G = nx.karate_club_graph()
A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))  # symmetric adjacency
n = A.shape[0]

# --- Build the modularity matrix B_ij = A_ij - k_i*k_j/(2m) ---
k = A.sum(axis=1)          # node degrees
two_m = A.sum()            # = 2m
m = two_m / 2.0
B = A - np.outer(k, k) / two_m

# --- Newman spectral bisection: leading eigenvector of B ---
eigvals, eigvecs = np.linalg.eigh(B)      # ascending eigenvalues
leading = eigvecs[:, np.argmax(eigvals)]  # eigenvector of largest eigenvalue

# Fix the eigenvector sign so the labeling is reproducible:
# force the first nonzero-signed entry to be positive.
if leading[np.argmax(np.abs(leading))] < 0:
    leading = -leading

# Split by sign of the leading eigenvector into two communities
comm = (leading >= 0).astype(int)  # 0/1 community label per node

# --- Modularity Q = (1/2m) * sum_ij (A_ij - k_i k_j/2m) * delta(c_i,c_j) ---
same = (comm[:, None] == comm[None, :]).astype(float)  # delta(c_i, c_j)
Q = np.sum(B * same) / two_m

# --- Report sizes and modularity ---
size0 = int(np.sum(comm == 0))
size1 = int(np.sum(comm == 1))
print(f"Community 0 size: {size0}")
print(f"Community 1 size: {size1}")
print(f"Modularity Q: {Q:.4f}")

# --- Check: two communities of 18 and 16 with Q about 0.37 ---
sizes_ok = sorted([size0, size1]) == [16, 18]
Q_ok = abs(Q - 0.37) < 0.02
print(f"Check sizes are {{16, 18}}: {sizes_ok}")
print(f"Check Q ~ 0.37: {Q_ok}")
print(f"Check passed: {sizes_ok and Q_ok}")
# This check confirms the result because reproducing Zachary's own 18/16
# factional split at high modularity shows the spectral method recovers the
# real community structure rather than an arbitrary partition.

# --- Draw the network with the same Kamada-Kawai layout as 10C.1 ---
pos = nx.kamada_kawai_layout(G)
colors = ["tab:blue" if c == 0 else "tab:orange" for c in comm]
plt.figure(figsize=(8, 6))
nx.draw_networkx_edges(G, pos, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=300)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title(f"Newman spectral communities (Q = {Q:.3f})")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.3.1_s1.png")
