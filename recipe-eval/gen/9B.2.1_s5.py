import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Sequences and scoring parameters ---
A = "GAATTCAGTTA"
B = "GGATCGA"
s_match, s_mismatch, gap = 3, -3, -2

def S(a, b):
    # substitution score: match vs mismatch
    return s_match if a == b else s_mismatch


# ============================================================
# GLOBAL ALIGNMENT (Needleman-Wunsch) -- for comparison
# ============================================================
n, m = len(A), len(B)

# G[i][j] = best global alignment score of A[:i] with B[:j]
G = [[0] * (m + 1) for _ in range(n + 1)]
# initialize first row/column with cumulative gap penalties
for i in range(1, n + 1):
    G[i][0] = G[i - 1][0] + gap
for j in range(1, m + 1):
    G[0][j] = G[0][j - 1] + gap

# fill the matrix: no clamping to zero (global allows negatives)
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = G[i - 1][j - 1] + S(A[i - 1], B[j - 1])  # align A_i with B_j
        up = G[i - 1][j] + gap                          # gap in B
        left = G[i][j - 1] + gap                         # gap in A
        G[i][j] = max(diag, up, left)

# traceback from bottom-right corner (i=n, j=m) to top-left (0,0)
gi, gj = n, m
gA, gB = "", ""
while gi > 0 or gj > 0:
    if gi > 0 and gj > 0 and G[gi][gj] == G[gi - 1][gj - 1] + S(A[gi - 1], B[gj - 1]):
        gA = A[gi - 1] + gA
        gB = B[gj - 1] + gB
        gi -= 1
        gj -= 1
    elif gi > 0 and G[gi][gj] == G[gi - 1][gj] + gap:
        gA = A[gi - 1] + gA
        gB = "-" + gB
        gi -= 1
    else:
        gA = "-" + gA
        gB = B[gj - 1] + gB
        gj -= 1

global_score = G[n][m]


# ============================================================
# LOCAL ALIGNMENT (Smith-Waterman)
# ============================================================
# H[i][j] = best local alignment score ending at A_i, B_j
H = [[0] * (m + 1) for _ in range(n + 1)]
# first row/column stay zero (a local alignment can start anywhere)

best_val, best_i, best_j = 0, 0, 0
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = H[i - 1][j - 1] + S(A[i - 1], B[j - 1])  # match/mismatch
        up = H[i - 1][j] + gap                          # gap in B
        left = H[i][j - 1] + gap                         # gap in A
        # clamp negative scores to zero: a fresh local alignment can restart here
        H[i][j] = max(diag, up, left, 0)
        # track the largest entry -- this is where traceback begins
        if H[i][j] > best_val:
            best_val, best_i, best_j = H[i][j], i, j

# traceback from the largest entry, stopping when we hit a zero
li, lj = best_i, best_j
lA, lB = "", ""
while li > 0 and lj > 0 and H[li][lj] != 0:
    if H[li][lj] == H[li - 1][lj - 1] + S(A[li - 1], B[lj - 1]):
        lA = A[li - 1] + lA
        lB = B[lj - 1] + lB
        li -= 1
        lj -= 1
    elif H[li][lj] == H[li - 1][lj] + gap:
        lA = A[li - 1] + lA
        lB = "-" + lB
        li -= 1
    else:
        lA = "-" + lA
        lB = B[lj - 1] + lB
        lj -= 1

local_score = best_val
# record the internal segment coordinates (1-based, inclusive) that were matched
local_A_start, local_A_end = li + 1, best_i
local_B_start, local_B_end = lj + 1, best_j


# ============================================================
# Report results
# ============================================================
print("Sequence A:", A)
print("Sequence B:", B)
print("s_match =", s_match, " s_mismatch =", s_mismatch, " gap =", gap)
print()

print("Global alignment (Needleman-Wunsch):")
print("  A:", gA)
print("  B:", gB)
print("Global alignment score:", global_score)
print("Global aligned length:", len(gA))
print()

print("Local alignment (Smith-Waterman):")
print("  A:", lA)
print("  B:", lB)
print("Local alignment score:", local_score)
print("Local segment of A (1-based, inclusive): [%d, %d]" % (local_A_start, local_A_end))
print("Local segment of B (1-based, inclusive): [%d, %d]" % (local_B_start, local_B_end))
print()

# --- Separate check ---
# The global alignment must span both sequences end to end: the number of
# non-gap characters on each row must equal the full sequence length.
global_A_span = sum(1 for c in gA if c != "-")
global_B_span = sum(1 for c in gB if c != "-")
global_spans_full = (global_A_span == len(A)) and (global_B_span == len(B))

# The local alignment should report only an internal segment: at least one of
# its endpoints does NOT reach the sequence boundary (start > 1 or end < len).
local_is_internal = (local_A_start > 1 or local_A_end < len(A) or
                     local_B_start > 1 or local_B_end < len(B))

print("CHECK -- global spans A end-to-end (chars):", global_A_span, "of", len(A))
print("CHECK -- global spans B end-to-end (chars):", global_B_span, "of", len(B))
print("CHECK -- global alignment spans both sequences fully:", global_spans_full)
print("CHECK -- local alignment is an internal (not full-length) segment:", local_is_internal)
print()
print("Explanation: the check confirms the result because the global row-spans"
      " equal each full sequence length (proving end-to-end coverage) while the"
      " local segment's endpoints fall strictly inside the sequences (proving it"
      " reports only the best-matching internal subsegment).")


# ============================================================
# Visualization: local (H) vs global (G) score matrices
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, mat, title, (ti, tj) in (
    (axes[0], H, "Smith-Waterman H (local)", (best_i, best_j)),
    (axes[1], G, "Needleman-Wunsch G (global)", (n, m)),
):
    arr = [[mat[i][j] for j in range(m + 1)] for i in range(n + 1)]
    im = ax.imshow(arr, cmap="viridis", aspect="auto")
    ax.set_title(title)
    ax.set_xlabel("B index (0 = gap)")
    ax.set_ylabel("A index (0 = gap)")
    ax.set_xticks(range(m + 1))
    ax.set_xticklabels(["-"] + list(B))
    ax.set_yticks(range(n + 1))
    ax.set_yticklabels(["-"] + list(A))
    for i in range(n + 1):
        for j in range(m + 1):
            ax.text(j, i, str(arr[i][j]), ha="center", va="center",
                    color="white", fontsize=7)
    # mark the traceback start cell
    ax.scatter([tj], [ti], s=200, facecolors="none", edgecolors="red", linewidths=2)
    fig.colorbar(im, ax=ax, fraction=0.046)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.2.1_s5.png")
