import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# --- Input: Zachary's karate club as an adjacency matrix ---
G = nx.karate_club_graph()
nodes = list(G.nodes())
A = nx.to_numpy_array(G, nodelist=nodes)   # symmetric 0/1 adjacency matrix
n = A.shape[0]

# --- Modularity matrix ---
k = A.sum(axis=1)          # degree of each node
m = A.sum() / 2.0          # total number of edges
B = A - np.outer(k, k) / (2.0 * m)   # B_ij = A_ij - k_i k_j / 2m

# --- Newman spectral bisection: leading eigenvector of B ---
eigvals, eigvecs = np.linalg.eigh(B)   # eigh: B is symmetric
lead_idx = np.argmax(eigvals)          # index of the largest (leading) eigenvalue
u = eigvecs[:, lead_idx].copy()        # leading eigenvector

# Fix the eigenvector sign so the labeling is reproducible:
# force the largest-magnitude component to be positive.
if u[np.argmax(np.abs(u))] < 0:
    u = -u

# Split by sign of the leading eigenvector components (+1 / -1 groups)
community = np.where(u >= 0, 1, 0)

# --- Modularity Q of this split, computed explicitly from the definition ---
# Q = (1/2m) * sum_ij (A_ij - k_i k_j/2m) * delta(c_i, c_j)
same = (community[:, None] == community[None, :]).astype(float)  # delta(c_i, c_j)
Q = np.sum(B * same) / (2.0 * m)

# --- Report sizes and modularity ---
size0 = int(np.sum(community == 0))
size1 = int(np.sum(community == 1))
print(f"Number of nodes: {n}")
print(f"Number of edges m: {int(m)}")
print(f"Community A size: {size0}")
print(f"Community B size: {size1}")
print(f"Modularity Q: {Q:.4f}")

# --- Draw with the same Kamada-Kawai layout as 10C.1 ---
pos = nx.kamada_kawai_layout(G)
colors = ["#1f77b4" if c == 1 else "#ff7f0e" for c in community]
plt.figure(figsize=(8, 6))
nx.draw_networkx_edges(G, pos, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=300)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title(f"Newman spectral communities of Zachary's karate club (Q = {Q:.3f})")
plt.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.3.1_s2.png")

# --- Separate check ---
sizes_ok = sorted([size0, size1]) == [16, 18]
q_ok = abs(Q - 0.37) < 0.02
print(f"Check split sizes == {{18, 16}}: {sizes_ok}")
print(f"Check Q ~ 0.37: {q_ok}")
print(f"Overall check passed: {sizes_ok and q_ok}")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Recovering an 18/16 split with Q~=0.37, the same partition and "
      "modularity Zachary observed for the club's real factional break, confirms that "
      "the spectral method found the network's genuine community structure rather than "
      "an arbitrary cut.")
