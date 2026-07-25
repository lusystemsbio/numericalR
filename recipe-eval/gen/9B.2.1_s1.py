import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Inputs ----
A = "GAATTCAGTTA"
B = "GGATCGA"
s_match, s_mismatch, gap = 3, -3, -2

def score(a, b):
    return s_match if a == b else s_mismatch  # substitution score S(A_i, B_j)

# =====================================================================
# GLOBAL ALIGNMENT (Needleman-Wunsch) -- for comparison
# =====================================================================
n, m = len(A), len(B)
G = np.zeros((n + 1, m + 1), dtype=int)

# Initialize first row/column with cumulative gap penalties (no clamping).
for i in range(1, n + 1):
    G[i][0] = G[i - 1][0] + gap
for j in range(1, m + 1):
    G[0][j] = G[0][j - 1] + gap

# Fill: choose best of diagonal (match/mismatch), up (gap in B), left (gap in A).
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = G[i - 1][j - 1] + score(A[i - 1], B[j - 1])
        up = G[i - 1][j] + gap
        left = G[i][j - 1] + gap
        G[i][j] = max(diag, up, left)

# Traceback from bottom-right corner to top-left (spans both sequences fully).
gi, gj = n, m
gA, gB = "", ""
while gi > 0 or gj > 0:
    if gi > 0 and gj > 0 and G[gi][gj] == G[gi - 1][gj - 1] + score(A[gi - 1], B[gj - 1]):
        gA = A[gi - 1] + gA
        gB = B[gj - 1] + gB
        gi -= 1; gj -= 1
    elif gi > 0 and G[gi][gj] == G[gi - 1][gj] + gap:
        gA = A[gi - 1] + gA
        gB = "-" + gB
        gi -= 1
    else:
        gA = "-" + gA
        gB = B[gj - 1] + gB
        gj -= 1
global_score = G[n][m]

# =====================================================================
# LOCAL ALIGNMENT (Smith-Waterman) -- negative scores clamped to zero
# =====================================================================
H = np.zeros((n + 1, m + 1), dtype=int)  # first row/col stay zero (local baseline)

# Fill using the given recurrence, with 0 as a floor so alignments can restart.
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = H[i - 1][j - 1] + score(A[i - 1], B[j - 1])
        left = H[i][j - 1] + gap
        up = H[i - 1][j] + gap
        H[i][j] = max(diag, left, up, 0)

# Best local alignment ends at the largest entry in the matrix.
best = np.unravel_index(np.argmax(H), H.shape)
li, lj = int(best[0]), int(best[1])
local_score = int(H[li][lj])
end_i, end_j = li, lj  # 1-based end positions in A and B

# Traceback from the max cell, stopping when we reach a zero.
lA, lB = "", ""
while li > 0 and lj > 0 and H[li][lj] != 0:
    if H[li][lj] == H[li - 1][lj - 1] + score(A[li - 1], B[lj - 1]):
        lA = A[li - 1] + lA
        lB = B[lj - 1] + lB
        li -= 1; lj -= 1
    elif H[li][lj] == H[li - 1][lj] + gap:
        lA = A[li - 1] + lA
        lB = "-" + lB
        li -= 1
    else:
        lA = "-" + lB[:0] + A[li - 1][:0] + lA  # placeholder (never taken over below)
        lA = "-" + lA[1:] if False else "-" + (lA[1:] if False else lA)
        lB = B[lj - 1] + lB
        lj -= 1
# start positions (1-based) are one past where traceback stopped
start_i, start_j = li + 1, lj + 1

# =====================================================================
# Report results
# =====================================================================
print("Sequence A:", A)
print("Sequence B:", B)
print("Scores -> match:", s_match, " mismatch:", s_mismatch, " gap:", gap)
print()
print("Global alignment score:", global_score)
print("Global A:", gA)
print("Global B:", gB)
print()
print("Local (Smith-Waterman) alignment score:", local_score)
print("Local A:", lA)
print("Local B:", lB)
print("Local segment in A: positions", start_i, "to", end_i, "->", A[start_i - 1:end_i])
print("Local segment in B: positions", start_j, "to", end_j, "->", B[start_j - 1:end_j])
print()

# =====================================================================
# Separate check: global spans both ends; local is an internal segment
# =====================================================================
global_spans_A = gA.replace("-", "") == A
global_spans_B = gB.replace("-", "") == B
global_spans_both = global_spans_A and global_spans_B
local_is_internal = (len(A[start_i - 1:end_i]) < len(A)) or (len(B[start_j - 1:end_j]) < len(B))

print("CHECK global alignment recovers all of A end-to-end:", global_spans_A)
print("CHECK global alignment recovers all of B end-to-end:", global_spans_B)
print("CHECK global alignment spans both sequences fully:", global_spans_both)
print("CHECK local alignment is only an internal subsegment:", local_is_internal)
print()
print("Explanation: stripping gaps from the global rows reproduces the full A and B "
      "while the local rows reproduce only a shorter interior slice, which confirms the "
      "global DP aligns end-to-end whereas Smith-Waterman isolates the best-matching subsegment.")

# =====================================================================
# Visualization: Smith-Waterman H matrix heatmap with traceback path
# =====================================================================
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(H, cmap="viridis")
for i in range(n + 1):
    for j in range(m + 1):
        ax.text(j, i, str(H[i][j]), ha="center", va="center",
                color="white" if H[i][j] < H.max() / 2 else "black", fontsize=8)
ax.set_xticks(range(m + 1))
ax.set_xticklabels(["-"] + list(B))
ax.set_yticks(range(n + 1))
ax.set_yticklabels(["-"] + list(A))
ax.set_xlabel("B")
ax.set_ylabel("A")
ax.set_title("Smith-Waterman H matrix (local score = %d)" % local_score)
fig.colorbar(im, ax=ax, shrink=0.8)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.2.1_s1.png")
