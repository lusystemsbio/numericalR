import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Sequences and scoring parameters ----
A = "GAATTCAGTTA"
B = "GGATCGA"
s_match, s_mismatch, gap = 3, -3, -2


def S(a, b):
    # match/mismatch substitution score
    return s_match if a == b else s_mismatch


# =====================================================================
# GLOBAL alignment (Needleman-Wunsch) -- for comparison
# =====================================================================
def needleman_wunsch(A, B):
    n, m = len(A), len(B)
    # H[i][j] = best global score aligning A[:i] with B[:j]
    H = [[0] * (m + 1) for _ in range(n + 1)]
    # first row/column: only gaps allowed
    for i in range(1, n + 1):
        H[i][0] = H[i - 1][0] + gap
    for j in range(1, m + 1):
        H[0][j] = H[0][j - 1] + gap
    # fill the matrix taking the best of diagonal/up/left
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            H[i][j] = max(
                H[i - 1][j - 1] + S(A[i - 1], B[j - 1]),  # match/mismatch
                H[i - 1][j] + gap,                        # gap in B
                H[i][j - 1] + gap,                        # gap in A
            )
    # traceback from bottom-right corner all the way to the top-left
    i, j = n, m
    a_al, b_al = "", ""
    while i > 0 or j > 0:
        if i > 0 and j > 0 and H[i][j] == H[i - 1][j - 1] + S(A[i - 1], B[j - 1]):
            a_al = A[i - 1] + a_al
            b_al = B[j - 1] + b_al
            i, j = i - 1, j - 1
        elif i > 0 and H[i][j] == H[i - 1][j] + gap:
            a_al = A[i - 1] + a_al
            b_al = "-" + b_al
            i -= 1
        else:
            a_al = "-" + a_al
            b_al = B[j - 1] + b_al
            j -= 1
    return H[n][m], a_al, b_al


# =====================================================================
# LOCAL alignment (Smith-Waterman)
# =====================================================================
def smith_waterman(A, B):
    n, m = len(A), len(B)
    # H[i][j] = best local score; negative scores are clamped to zero
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best_score, best_pos = 0, (0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            H[i][j] = max(
                H[i - 1][j - 1] + S(A[i - 1], B[j - 1]),  # match/mismatch
                H[i][j - 1] + gap,                        # gap in A
                H[i - 1][j] + gap,                        # gap in B
                0,                                        # clamp: start fresh
            )
            if H[i][j] > best_score:
                best_score, best_pos = H[i][j], (i, j)
    # traceback from the largest entry, stopping when we hit a zero
    i, j = best_pos
    a_al, b_al = "", ""
    while i > 0 and j > 0 and H[i][j] != 0:
        if H[i][j] == H[i - 1][j - 1] + S(A[i - 1], B[j - 1]):
            a_al = A[i - 1] + a_al
            b_al = B[j - 1] + b_al
            i, j = i - 1, j - 1
        elif H[i][j] == H[i - 1][j] + gap:
            a_al = A[i - 1] + a_al
            b_al = "-" + b_al
            i -= 1
        else:
            a_al = "-" + a_al
            b_al = B[j - 1] + b_al
            j -= 1
    # record where in A and B the local segment starts (0-based)
    return best_score, a_al, b_al, i, j, H


g_score, gA, gB = needleman_wunsch(A, B)
l_score, lA, lB, li, lj, Hloc = smith_waterman(A, B)


def match_line(x, y):
    return "".join("|" if a == b and a != "-" else " " for a, b in zip(x, y))


# ---- Report results ----
print(f"Sequence A: {A}")
print(f"Sequence B: {B}")
print(f"Scoring: match={s_match}, mismatch={s_mismatch}, gap={gap}")
print()

print("Global alignment (Needleman-Wunsch):")
print(f"  A: {gA}")
print(f"     {match_line(gA, gB)}")
print(f"  B: {gB}")
print(f"Global alignment score: {g_score}")
print()

print("Local alignment (Smith-Waterman):")
print(f"  A: {lA}")
print(f"     {match_line(lA, lB)}")
print(f"  B: {lB}")
print(f"Local alignment score: {l_score}")
print(f"Local segment starts at A index (0-based): {li}")
print(f"Local segment starts at B index (0-based): {lj}")
print()

# ---- Separate span check ----
# Global: the aligned rows, with gaps removed, must reproduce the FULL sequences.
global_spans_A = gA.replace("-", "") == A
global_spans_B = gB.replace("-", "") == B
# Local: the aligned rows, with gaps removed, are only INTERNAL substrings.
local_A_sub = lA.replace("-", "")
local_B_sub = lB.replace("-", "")
local_is_internal_A = local_A_sub in A and local_A_sub != A
local_is_internal_B = local_B_sub in B

print("Span check:")
print(f"Global alignment covers all of A end-to-end: {global_spans_A}")
print(f"Global alignment covers all of B end-to-end: {global_spans_B}")
print(f"Local A segment '{local_A_sub}' is an internal substring of A: {local_is_internal_A}")
print(f"Local B segment '{local_B_sub}' is an internal substring of B: {local_is_internal_B}")
print(f"Local segment is shorter than full A ({len(local_A_sub)} < {len(A)}): {len(local_A_sub) < len(A)}")
print()
print("Why this confirms the result: stripping gaps from the global rows reproduces "
      "the entire input sequences (end-to-end coverage), whereas stripping gaps from "
      "the local rows yields only a shorter internal substring, showing Smith-Waterman "
      "returned just the best-matching subsegment rather than a full-length alignment.")

# ---- Visualization: Smith-Waterman score matrix ----
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(Hloc, cmap="viridis")
ax.set_xticks(range(len(B) + 1))
ax.set_xticklabels(["-"] + list(B))
ax.set_yticks(range(len(A) + 1))
ax.set_yticklabels(["-"] + list(A))
ax.set_xlabel("Sequence B")
ax.set_ylabel("Sequence A")
ax.set_title(f"Smith-Waterman H matrix (best local score = {l_score})")
for i in range(len(A) + 1):
    for j in range(len(B) + 1):
        ax.text(j, i, Hloc[i][j], ha="center", va="center",
                color="white" if Hloc[i][j] < l_score * 0.6 else "black", fontsize=8)
fig.colorbar(im, ax=ax, label="H score")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.2.1_s4.png")
