import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Sequences to align ----
A = "GAATTCAGTTA"
B = "GGATCGA"

def needleman_wunsch(A, B, s_match, s_mismatch, gap):
    """Global (end-to-end) alignment via the Needleman-Wunsch DP."""
    n, m = len(A), len(B)
    # F is the (n+1) x (m+1) score matrix
    F = np.zeros((n + 1, m + 1), dtype=int)
    # Initialize first column/row: aligning a prefix against all gaps
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + gap
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + gap
    # Substitution score S(a, b)
    def S(a, b):
        return s_match if a == b else s_mismatch
    # Fill the matrix using F_{i,j} = max(diag+S, left+gap, up+gap)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1][j - 1] + S(A[i - 1], B[j - 1])  # match/mismatch
            left = F[i][j - 1] + gap                         # gap in A
            up = F[i - 1][j] + gap                           # gap in B
            F[i][j] = max(diag, left, up)
    # ---- Traceback from the last cell (n, m) ----
    ai, bj = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + S(A[i - 1], B[j - 1]):
            ai.append(A[i - 1]); bj.append(B[j - 1]); i -= 1; j -= 1   # diagonal
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            ai.append(A[i - 1]); bj.append('-'); i -= 1                # gap in B
        else:
            ai.append('-'); bj.append(B[j - 1]); j -= 1               # gap in A
    ai.reverse(); bj.reverse()
    aln_A = "".join(ai)
    aln_B = "".join(bj)
    return F, F[n][m], aln_A, aln_B

def count_matches(aln_A, aln_B):
    """Count aligned columns that are identical (non-gap matches)."""
    return sum(1 for x, y in zip(aln_A, aln_B) if x == y and x != '-')

def count_gaps(aln_A, aln_B):
    """Count gap columns in either sequence."""
    return sum(1 for x, y in zip(aln_A, aln_B) if x == '-' or y == '-')

# ---- Scheme 1: match-only (reward matches, no penalties) ----
F1, score1, a1, b1 = needleman_wunsch(A, B, s_match=1, s_mismatch=0, gap=0)
matches1 = count_matches(a1, b1)
gaps1 = count_gaps(a1, b1)

# ---- Scheme 2: match=3, mismatch=-3, gap=-2 ----
F2, score2, a2, b2 = needleman_wunsch(A, B, s_match=3, s_mismatch=-3, gap=-2)
matches2 = count_matches(a2, b2)
gaps2 = count_gaps(a2, b2)

# ---- Report ----
print("Sequence A:", A)
print("Sequence B:", B)
print()
print("=== Scheme 1: match-only (s_match=1, s_mismatch=0, gap=0) ===")
print("Optimal score:", score1)
print("Alignment A:", a1)
print("Alignment B:", b1)
print("Number of matches:", matches1)
print("Number of gaps:", gaps1)
print()
print("=== Scheme 2: (s_match=3, s_mismatch=-3, gap=-2) ===")
print("Optimal score:", score2)
print("Alignment A:", a2)
print("Alignment B:", b2)
print("Number of matches:", matches2)
print("Number of gaps:", gaps2)
print()

# ---- Separate check ----
# Under match-only scoring the score equals the match count, so the DP
# maximizes matches with no regard for how many gaps that costs.
print("=== Check ===")
print("Scheme 1 score equals its match count (score==matches):", score1 == matches1)
print("Scheme 1 achieves >= matches than Scheme 2:", matches1 >= matches2)
print("Scheme 1 gaps:", gaps1, " Scheme 2 gaps:", gaps2)
print("Alignments differ between schemes:", (a1, b1) != (a2, b2))
print("Explanation: because the match-only score IS the match count, its optimum "
      "must maximize matches irrespective of gaps, whereas the penalized scheme "
      "trades some matches away to avoid gap cost, yielding a different alignment.")

# ---- Visualization of the two score matrices ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, F, title, seq in [
    (axes[0], F1, "Scheme 1 (match-only)", (A, B)),
    (axes[1], F2, "Scheme 2 (3,-3,-2)", (A, B)),
]:
    im = ax.imshow(F, cmap="viridis")
    ax.set_title(f"{title}\nscore={F[len(A)][len(B)]}")
    ax.set_xticks(range(len(B) + 1))
    ax.set_xticklabels(['-'] + list(B))
    ax.set_yticks(range(len(A) + 1))
    ax.set_yticklabels(['-'] + list(A))
    ax.set_xlabel("B"); ax.set_ylabel("A")
    for i in range(F.shape[0]):
        for j in range(F.shape[1]):
            ax.text(j, i, str(F[i][j]), ha="center", va="center",
                    color="white", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.1.1_s1.png")
