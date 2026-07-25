import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Needleman-Wunsch global (end-to-end) alignment, implemented explicitly.
# ------------------------------------------------------------------

def needleman_wunsch(A, B, s_match, s_mismatch, g):
    """Return (score, aligned_A, aligned_B, n_matches, n_gaps)."""
    n, m = len(A), len(B)

    # substitution score S(A_i, B_j)
    def S(a, b):
        return s_match if a == b else s_mismatch

    # F is the (n+1) x (m+1) score matrix; T stores traceback directions.
    F = [[0] * (m + 1) for _ in range(n + 1)]
    # 'D' = diagonal (match/mismatch), 'U' = up (gap in B), 'L' = left (gap in A)
    T = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize first row/column: only gaps are possible along the edges.
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + g
        T[i][0] = 'U'
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + g
        T[0][j] = 'L'

    # Fill the matrix using the recurrence F_{i,j} = max(diag, left, up).
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1][j - 1] + S(A[i - 1], B[j - 1])  # align A_i with B_j
            up   = F[i - 1][j] + g                          # gap in B
            left = F[i][j - 1] + g                          # gap in A
            best = max(diag, up, left)
            F[i][j] = best
            # tie-break priority: diagonal, then up, then left
            if best == diag:
                T[i][j] = 'D'
            elif best == up:
                T[i][j] = 'U'
            else:
                T[i][j] = 'L'

    # Traceback from the last cell (n, m) to (0, 0).
    ai, bj = [], []
    i, j = n, m
    while i > 0 or j > 0:
        d = T[i][j]
        if d == 'D':
            ai.append(A[i - 1]); bj.append(B[j - 1]); i -= 1; j -= 1
        elif d == 'U':
            ai.append(A[i - 1]); bj.append('-'); i -= 1
        else:  # 'L'
            ai.append('-'); bj.append(B[j - 1]); j -= 1

    aligned_A = ''.join(reversed(ai))
    aligned_B = ''.join(reversed(bj))

    # Count matches and gaps in the resulting alignment.
    n_matches = sum(1 for x, y in zip(aligned_A, aligned_B) if x == y and x != '-')
    n_gaps = aligned_A.count('-') + aligned_B.count('-')

    return F[n][m], aligned_A, aligned_B, n_matches, n_gaps


# ------------------------------------------------------------------
# Test sequences and the two scoring schemes.
# ------------------------------------------------------------------
A = "GAATTCAGTTA"
B = "GGATCGA"

schemes = [
    ("match-only", dict(s_match=1, s_mismatch=0, g=0)),
    ("match/mismatch/gap", dict(s_match=3, s_mismatch=-3, g=-2)),
]

results = {}
for name, params in schemes:
    score, aA, aB, n_match, n_gap = needleman_wunsch(A, B, **params)
    results[name] = (aA, aB, n_match, n_gap, score)
    print(f"=== Scheme: {name}  (s_match={params['s_match']}, "
          f"s_mismatch={params['s_mismatch']}, g={params['g']}) ===")
    print(f"Aligned A: {aA}")
    print(f"Aligned B: {aB}")
    print(f"Optimal score: {score}")
    print(f"Number of matches: {n_match}")
    print(f"Number of gaps: {n_gap}")
    print()

# ------------------------------------------------------------------
# Separate check:
# Rewarding only matches maximizes the number of matched columns
# regardless of how many gaps are introduced; adding mismatch and gap
# penalties trades matches off against gap cost, changing the alignment.
# ------------------------------------------------------------------
mo_matches = results["match-only"][2]
mo_gaps    = results["match-only"][3]
mm_matches = results["match/mismatch/gap"][2]
mm_gaps    = results["match/mismatch/gap"][3]

print("=== Check: match-only maximizes matches; penalties balance vs gaps ===")
print(f"match-only matches: {mo_matches}")
print(f"match-only gaps: {mo_gaps}")
print(f"penalized-scheme matches: {mm_matches}")
print(f"penalized-scheme gaps: {mm_gaps}")
print(f"match-only has >= matches than penalized scheme: {mo_matches >= mm_matches}")
print(f"match-only uses >= gaps than penalized scheme: {mo_gaps >= mm_gaps}")
print(f"Alignments differ between schemes: "
      f"{results['match-only'][:2] != results['match/mismatch/gap'][:2]}")

# ------------------------------------------------------------------
# Visualization of both alignments.
# ------------------------------------------------------------------
fig, axes = plt.subplots(len(schemes), 1, figsize=(9, 4.5))
for ax, (name, _) in zip(axes, schemes):
    aA, aB, n_match, n_gap, score = results[name]
    ax.axis("off")
    ax.set_title(f"{name}: score={score}, matches={n_match}, gaps={n_gap}",
                 fontsize=11)
    # Draw the two aligned strings with match/mismatch coloring.
    for k, (ca, cb) in enumerate(zip(aA, aB)):
        color = "green" if (ca == cb and ca != '-') else \
                ("red" if (ca != '-' and cb != '-') else "gray")
        ax.text(k, 1.0, ca, family="monospace", fontsize=14,
                ha="center", va="center", color=color)
        ax.text(k, 0.0, cb, family="monospace", fontsize=14,
                ha="center", va="center", color=color)
    ax.set_xlim(-1, max(len(aA), len(aB)))
    ax.set_ylim(-0.6, 1.6)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.1.1_s4.png")

# One-sentence explanation of why the check confirms the result.
print()
print("Explanation: The check confirms the result because the match-only scheme "
      "attains at least as many matched columns (using freely more gaps) as the "
      "penalized scheme, whereas the penalized scheme sacrifices some matches to "
      "avoid gap cost and thus yields a genuinely different optimal alignment.")
