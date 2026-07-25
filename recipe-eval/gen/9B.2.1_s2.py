import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Sequences and scoring parameters
# ---------------------------------------------------------------
A = "GAATTCAGTTA"
B = "GGATCGA"
s_match, s_mismatch, gap = 3, -3, -2

def S(a, b):
    # substitution score: match vs. mismatch
    return s_match if a == b else s_mismatch

# ---------------------------------------------------------------
# Smith-Waterman LOCAL alignment
# H[i][j] = max(diagonal + S, left + gap, up + gap, 0)  (clamp to 0)
# ---------------------------------------------------------------
def smith_waterman(A, B):
    n, m = len(A), len(B)
    H = [[0] * (m + 1) for _ in range(n + 1)]   # score matrix, first row/col = 0
    P = [[0] * (m + 1) for _ in range(n + 1)]   # pointer: 0=stop,1=diag,2=up,3=left
    best_val, best_pos = 0, (0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + S(A[i - 1], B[j - 1])
            up   = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            # take the largest option, but never below zero (local: fresh restart)
            H[i][j] = max(diag, up, left, 0)
            if H[i][j] == 0:
                P[i][j] = 0
            elif H[i][j] == diag:
                P[i][j] = 1
            elif H[i][j] == up:
                P[i][j] = 2
            else:
                P[i][j] = 3
            # track the global maximum entry (start of traceback)
            if H[i][j] > best_val:
                best_val, best_pos = H[i][j], (i, j)
    # traceback from the largest entry, stop when we hit a zero
    i, j = best_pos
    a_al, b_al = [], []
    while i > 0 and j > 0 and H[i][j] != 0:
        if P[i][j] == 1:      # diagonal -> aligned pair
            a_al.append(A[i - 1]); b_al.append(B[j - 1]); i -= 1; j -= 1
        elif P[i][j] == 2:    # up -> gap in B
            a_al.append(A[i - 1]); b_al.append("-"); i -= 1
        else:                 # left -> gap in A
            a_al.append("-"); b_al.append(B[j - 1]); j -= 1
    a_al.reverse(); b_al.reverse()
    return H, best_val, "".join(a_al), "".join(b_al)

# ---------------------------------------------------------------
# Needleman-Wunsch GLOBAL alignment (for comparison; same S and gap)
# G[i][j] = max(diag + S, up + gap, left + gap)  -- NO clamping to zero
# ---------------------------------------------------------------
def needleman_wunsch(A, B):
    n, m = len(A), len(B)
    G = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]  # 1=diag,2=up,3=left
    for i in range(1, n + 1):
        G[i][0] = i * gap; P[i][0] = 2          # leading gaps down first column
    for j in range(1, m + 1):
        G[0][j] = j * gap; P[0][j] = 3          # leading gaps across first row
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = G[i - 1][j - 1] + S(A[i - 1], B[j - 1])
            up   = G[i - 1][j] + gap
            left = G[i][j - 1] + gap
            G[i][j] = max(diag, up, left)
            P[i][j] = 1 if G[i][j] == diag else (2 if G[i][j] == up else 3)
    # traceback from bottom-right corner all the way to (0,0): spans both ends
    i, j = n, m
    a_al, b_al = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and P[i][j] == 1:
            a_al.append(A[i - 1]); b_al.append(B[j - 1]); i -= 1; j -= 1
        elif i > 0 and P[i][j] == 2:
            a_al.append(A[i - 1]); b_al.append("-"); i -= 1
        else:
            a_al.append("-"); b_al.append(B[j - 1]); j -= 1
    a_al.reverse(); b_al.reverse()
    return G, G[n][m], "".join(a_al), "".join(b_al)

# ---------------------------------------------------------------
# Run both algorithms
# ---------------------------------------------------------------
H, local_score, la_A, la_B = smith_waterman(A, B)
G, global_score, ga_A, ga_B = needleman_wunsch(A, B)

# match/middle line helper
def midline(x, y):
    return "".join("|" if (a == b and a != "-") else " " for a, b in zip(x, y))

# ---------------------------------------------------------------
# Report results
# ---------------------------------------------------------------
print("Sequence A:", A)
print("Sequence B:", B)
print("Scores: match =", s_match, " mismatch =", s_mismatch, " gap =", gap)
print()
print("LOCAL alignment score:", local_score)
print("Local alignment A:", la_A)
print("Local alignment  :", midline(la_A, la_B))
print("Local alignment B:", la_B)
print()
print("GLOBAL alignment score:", global_score)
print("Global alignment A:", ga_A)
print("Global alignment  :", midline(ga_A, ga_B))
print("Global alignment B:", ga_B)
print()

# ---------------------------------------------------------------
# Separate check: global spans both sequences end to end; local is internal
# ---------------------------------------------------------------
global_A_full = ga_A.replace("-", "") == A
global_B_full = ga_B.replace("-", "") == B
local_A_sub = la_A.replace("-", "") in A
local_B_sub = la_B.replace("-", "") in B
local_is_internal = (len(la_A.replace("-", "")) < len(A)) or (len(la_B.replace("-", "")) < len(B))

print("CHECK global covers all of A end to end:", global_A_full)
print("CHECK global covers all of B end to end:", global_B_full)
print("CHECK local A residues form a subsegment of A:", local_A_sub)
print("CHECK local B residues form a subsegment of B:", local_B_sub)
print("CHECK local alignment is a strictly internal segment (shorter than full):", local_is_internal)
print("Local aligned length (columns):", len(la_A))
print("Global aligned length (columns):", len(ga_A))

# ---------------------------------------------------------------
# Figure: side-by-side text panels of the two alignments
# ---------------------------------------------------------------
out_path = "/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.2.1_s2.png"
os.makedirs(os.path.dirname(out_path), exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
local_text = ("LOCAL (Smith-Waterman)\nscore = %d\n\n%s\n%s\n%s"
              % (local_score, la_A, midline(la_A, la_B), la_B))
global_text = ("GLOBAL (Needleman-Wunsch)\nscore = %d\n\n%s\n%s\n%s"
               % (global_score, ga_A, midline(ga_A, ga_B), ga_B))
for ax, title, text in zip(axes, ["Local", "Global"], [local_text, global_text]):
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.02, 0.5, text, family="monospace", fontsize=13, va="center")
fig.suptitle("Local vs. Global alignment of A and B")
plt.tight_layout()
plt.savefig(out_path)

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result
# ---------------------------------------------------------------
print()
print("Explanation: The check confirms the result because the global alignment's "
      "gap-removed rows exactly reproduce all of A and all of B (end-to-end coverage), "
      "whereas the local alignment's rows reproduce only contiguous internal subsegments "
      "of each sequence, which is precisely the distinction between global and local alignment.")
