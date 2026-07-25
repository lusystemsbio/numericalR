import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Needleman-Wunsch global alignment, implemented explicitly.
# ----------------------------------------------------------------------

def needleman_wunsch(A, B, s_match, s_mismatch, g):
    """Return (aligned_A, aligned_B, score, n_matches, n_gaps)."""
    n, m = len(A), len(B)

    # S(a,b): substitution score for aligning residue a with residue b
    def S(a, b):
        return s_match if a == b else s_mismatch

    # F is the (n+1) x (m+1) score matrix; F[i][j] scores prefixes A[:i], B[:j]
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize first row/column: aligning a prefix against all gaps
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + g
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + g

    # Fill the matrix using the recurrence
    #   F[i][j] = max(diag + S, left + g, up + g)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1][j - 1] + S(A[i - 1], B[j - 1])   # match / mismatch
            left = F[i][j - 1] + g                            # gap in A
            up   = F[i - 1][j] + g                            # gap in B
            F[i][j] = max(diag, left, up)

    # Traceback from the bottom-right cell to reconstruct the alignment
    ai, bi = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + S(A[i - 1], B[j - 1]):
            ai.append(A[i - 1]); bi.append(B[j - 1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] + g:
            ai.append(A[i - 1]); bi.append('-'); i -= 1     # gap in B
        else:
            ai.append('-'); bi.append(B[j - 1]); j -= 1     # gap in A

    aligned_A = ''.join(reversed(ai))
    aligned_B = ''.join(reversed(bi))

    # Count matches and gaps in the produced alignment
    n_matches = sum(1 for x, y in zip(aligned_A, aligned_B) if x == y and x != '-')
    n_gaps = aligned_A.count('-') + aligned_B.count('-')

    return aligned_A, aligned_B, F[n][m], n_matches, n_gaps


A = "GAATTCAGTTA"
B = "GGATCGA"

schemes = [
    ("match-only (s_match=1, s_mismatch=0, gap=0)", 1, 0, 0),
    ("penalized (s_match=3, s_mismatch=-3, gap=-2)", 3, -3, -2),
]

results = []
for name, sm, smm, g in schemes:
    aA, aB, score, nm, ng = needleman_wunsch(A, B, sm, smm, g)
    results.append((name, aA, aB, score, nm, ng))
    print(f"Scheme: {name}")
    print(f"  Aligned A:      {aA}")
    print(f"  Aligned B:      {aB}")
    print(f"  Optimal score:  {score}")
    print(f"  Matches:        {nm}")
    print(f"  Gaps:           {ng}")

# Separate check: match-only should maximize the number of matches.
match_only_matches = results[0][4]
penalized_matches = results[1][4]
print(f"Match-only match count:  {match_only_matches}")
print(f"Penalized match count:   {penalized_matches}")
print(f"Match-only maximizes matches (>= penalized): {match_only_matches >= penalized_matches}")
print(f"Alignments differ between schemes: {results[0][1:3] != results[1][1:3]}")

# Explanation (one sentence):
# Rewarding only matches makes the DP objective identical to counting matches,
# so its optimum is the maximum-match alignment; adding mismatch and gap
# penalties changes the objective, yielding a different optimal alignment that
# trades some matches against gap cost -- confirming the two schemes optimize
# different things.
print("Explanation: match-only scoring makes the DP objective equal to the "
      "match count (so it is maximized regardless of gaps), whereas adding "
      "mismatch/gap penalties changes the objective and thus the optimal "
      "alignment, confirming matches are balanced against gap cost.")

# ----------------------------------------------------------------------
# Visualization of the two alignments.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(len(results), 1, figsize=(10, 4))
for ax, (name, aA, aB, score, nm, ng) in zip(axes, results):
    ax.axis("off")
    ax.set_title(f"{name}\nscore={score}, matches={nm}, gaps={ng}", fontsize=9)
    for k, (x, y) in enumerate(zip(aA, aB)):
        color = "green" if x == y and x != '-' else ("red" if x != '-' and y != '-' else "gray")
        ax.text(k, 1, x, ha="center", va="center", fontsize=12, color=color, family="monospace")
        ax.text(k, 0, y, ha="center", va="center", fontsize=12, color=color, family="monospace")
    ax.set_xlim(-1, max(len(aA), len(aB)))
    ax.set_ylim(-0.5, 1.5)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.1.1_s5.png")
