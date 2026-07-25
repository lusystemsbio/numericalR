import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Sequences and scoring scheme
# ---------------------------------------------------------------
A = "GAATTCAGTTA"
B = "GGATCGA"
s_match, s_mismatch, gap = 3, -3, -2

def S(a, b):
    # substitution score: match vs mismatch
    return s_match if a == b else s_mismatch

# ===============================================================
# SMITH-WATERMAN  (local alignment, negatives clamped to zero)
# ===============================================================
n, m = len(A), len(B)
# H[i][j] = best local alignment score ending at A[i-1], B[j-1]
H = [[0] * (m + 1) for _ in range(n + 1)]
# pointer for traceback: 'd'=diagonal, 'u'=up (gap in B), 'l'=left (gap in A), '0'=stop
ptr = [[None] * (m + 1) for _ in range(n + 1)]

best_val, best_i, best_j = 0, 0, 0
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = H[i-1][j-1] + S(A[i-1], B[j-1])  # align A_i with B_j
        up   = H[i-1][j]   + gap                 # gap in B (skip A_i)
        left = H[i][j-1]   + gap                 # gap in A (skip B_j)
        cell = max(diag, up, left, 0)            # clamp negatives to 0
        H[i][j] = cell
        # record which move produced the cell (0 means a fresh start)
        if cell == 0:
            ptr[i][j] = '0'
        elif cell == diag:
            ptr[i][j] = 'd'
        elif cell == up:
            ptr[i][j] = 'u'
        else:
            ptr[i][j] = 'l'
        # track the global maximum entry (start point of traceback)
        if cell > best_val:
            best_val, best_i, best_j = cell, i, j

# traceback from the largest entry, stopping when we reach a zero cell
la_A, la_B = [], []
i, j = best_i, best_j
while i > 0 and j > 0 and H[i][j] != 0:
    move = ptr[i][j]
    if move == 'd':
        la_A.append(A[i-1]); la_B.append(B[j-1]); i -= 1; j -= 1
    elif move == 'u':
        la_A.append(A[i-1]); la_B.append('-');    i -= 1
    elif move == 'l':
        la_A.append('-');    la_B.append(B[j-1]); j -= 1
    else:
        break
local_A = "".join(reversed(la_A))
local_B = "".join(reversed(la_B))
local_score = best_val
# 1-based inclusive coordinates of the matched internal segment
local_A_start, local_A_end = i + 1, best_i
local_B_start, local_B_end = j + 1, best_j

# ===============================================================
# NEEDLEMAN-WUNSCH  (global alignment, for comparison; no clamping)
# ===============================================================
G = [[0] * (m + 1) for _ in range(n + 1)]
gptr = [[None] * (m + 1) for _ in range(n + 1)]
# initialize first row/column with accumulated gap penalties
for i in range(1, n + 1):
    G[i][0] = G[i-1][0] + gap; gptr[i][0] = 'u'
for j in range(1, m + 1):
    G[0][j] = G[0][j-1] + gap; gptr[0][j] = 'l'
for i in range(1, n + 1):
    for j in range(1, m + 1):
        diag = G[i-1][j-1] + S(A[i-1], B[j-1])
        up   = G[i-1][j]   + gap
        left = G[i][j-1]   + gap
        cell = max(diag, up, left)   # global: no zero clamp
        G[i][j] = cell
        if cell == diag:
            gptr[i][j] = 'd'
        elif cell == up:
            gptr[i][j] = 'u'
        else:
            gptr[i][j] = 'l'

# traceback from the bottom-right corner back to the origin (spans both ends)
ga_A, ga_B = [], []
i, j = n, m
while i > 0 or j > 0:
    move = gptr[i][j]
    if move == 'd':
        ga_A.append(A[i-1]); ga_B.append(B[j-1]); i -= 1; j -= 1
    elif move == 'u':
        ga_A.append(A[i-1]); ga_B.append('-');    i -= 1
    else:
        ga_A.append('-');    ga_B.append(B[j-1]); j -= 1
global_A = "".join(reversed(ga_A))
global_B = "".join(reversed(ga_B))
global_score = G[n][m]

# ---------------------------------------------------------------
# alignment "middle" markers for readability
# ---------------------------------------------------------------
def mid(a, b):
    return "".join('|' if x == y and x != '-' else ' ' for x, y in zip(a, b))

# ---------------------------------------------------------------
# Reported results
# ---------------------------------------------------------------
print("Sequence A: " + A)
print("Sequence B: " + B)
print("Scoring: match=%d, mismatch=%d, gap=%d" % (s_match, s_mismatch, gap))
print("")
print("=== Local alignment (Smith-Waterman) ===")
print("Local score: %d" % local_score)
print("A segment [%d..%d]: %s" % (local_A_start, local_A_end, local_A))
print("                    " + mid(local_A, local_B))
print("B segment [%d..%d]: %s" % (local_B_start, local_B_end, local_B))
print("")
print("=== Global alignment (Needleman-Wunsch) ===")
print("Global score: %d" % global_score)
print("A: " + global_A)
print("   " + mid(global_A, global_B))
print("B: " + global_B)
print("")

# ---------------------------------------------------------------
# Separate check: coverage of each sequence
# ---------------------------------------------------------------
global_A_letters = len(global_A.replace('-', ''))
global_B_letters = len(global_B.replace('-', ''))
local_A_letters  = len(local_A.replace('-', ''))
local_B_letters  = len(local_B.replace('-', ''))
print("=== Coverage check ===")
print("Global covers A end-to-end: %s (%d of %d letters)" %
      (global_A_letters == len(A), global_A_letters, len(A)))
print("Global covers B end-to-end: %s (%d of %d letters)" %
      (global_B_letters == len(B), global_B_letters, len(B)))
print("Local A span is internal (not full): %s (%d of %d letters, positions %d..%d)" %
      (local_A_letters < len(A), local_A_letters, len(A), local_A_start, local_A_end))
print("Local B span is internal (not full): %s (%d of %d letters, positions %d..%d)" %
      (local_B_letters < len(B), local_B_letters, len(B), local_B_start, local_B_end))
print("")
print("Check meaning: the global alignment consuming every letter of both A and B "
      "while the local alignment consumes only a strict interior sub-range confirms that "
      "each algorithm did its intended job (end-to-end vs. best internal segment).")

# ---------------------------------------------------------------
# Figure: Smith-Waterman scoring matrix H with the traceback path
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
Hmat = [[H[i][j] for j in range(m + 1)] for i in range(n + 1)]
im = ax.imshow(Hmat, cmap="Blues")
ax.set_xticks(range(m + 1)); ax.set_yticks(range(n + 1))
ax.set_xticklabels(['-'] + list(B)); ax.set_yticklabels(['-'] + list(A))
ax.set_xlabel("B"); ax.set_ylabel("A")
ax.set_title("Smith-Waterman H matrix (local score = %d)" % local_score)
for i in range(n + 1):
    for j in range(m + 1):
        ax.text(j, i, str(H[i][j]), ha="center", va="center", fontsize=8)
# overlay the traceback path
ti, tj = best_i, best_j
path_i, path_j = [], []
while ti > 0 and tj > 0 and H[ti][tj] != 0:
    path_i.append(ti); path_j.append(tj)
    mv = ptr[ti][tj]
    if mv == 'd': ti -= 1; tj -= 1
    elif mv == 'u': ti -= 1
    elif mv == 'l': tj -= 1
    else: break
ax.plot(path_j, path_i, color="red", marker="o", markersize=6, linewidth=2)
fig.colorbar(im, ax=ax, shrink=0.8)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.2.1_s3.png")
