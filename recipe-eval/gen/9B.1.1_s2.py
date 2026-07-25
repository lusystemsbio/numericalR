import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def needleman_wunsch(A, B, s_match, s_mismatch, gap):
    """Global alignment via the Needleman-Wunsch dynamic program."""
    n, m = len(A), len(B)

    # S(A_i, B_j): substitution score for a pair of residues
    def S(a, b):
        return s_match if a == b else s_mismatch

    # F is the (n+1) x (m+1) score matrix; ptr records the chosen move for traceback
    F = np.zeros((n + 1, m + 1))
    ptr = np.empty((n + 1, m + 1), dtype=object)

    # Initialize first column and row: aligning a prefix against all-gaps
    for i in range(1, n + 1):
        F[i, 0] = F[i - 1, 0] + gap
        ptr[i, 0] = "up"      # gap in B
    for j in range(1, m + 1):
        F[0, j] = F[0, j - 1] + gap
        ptr[0, j] = "left"    # gap in A

    # Fill the matrix using the recurrence F_{i,j} = max(diag+S, left+g, up+g)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1, j - 1] + S(A[i - 1], B[j - 1])  # match/mismatch
            left = F[i, j - 1] + gap                        # gap in A
            up = F[i - 1, j] + gap                          # gap in B
            best = max(diag, left, up)
            F[i, j] = best
            if best == diag:
                ptr[i, j] = "diag"
            elif best == left:
                ptr[i, j] = "left"
            else:
                ptr[i, j] = "up"

    # Traceback from the bottom-right cell back to the origin
    align_A, align_B = [], []
    i, j = n, m
    while i > 0 or j > 0:
        move = ptr[i, j]
        if move == "diag":
            align_A.append(A[i - 1]); align_B.append(B[j - 1]); i -= 1; j -= 1
        elif move == "left":
            align_A.append("-"); align_B.append(B[j - 1]); j -= 1
        else:  # up
            align_A.append(A[i - 1]); align_B.append("-"); i -= 1

    return "".join(reversed(align_A)), "".join(reversed(align_B)), F[n, m], F


def count_matches(a, b):
    """Number of aligned positions where residues are identical (not gaps)."""
    return sum(1 for x, y in zip(a, b) if x == y and x != "-")


A = "GAATTCAGTTA"
B = "GGATCGA"

schemes = [
    ("match-only (match=1, mismatch=0, gap=0)", (1, 0, 0)),
    ("penalized (match=3, mismatch=-3, gap=-2)", (3, -3, -2)),
]

results = {}
for name, (sm, smm, g) in schemes:
    a, b, score, F = needleman_wunsch(A, B, sm, smm, g)
    matches = count_matches(a, b)
    gaps = a.count("-") + b.count("-")
    results[name] = (a, b, score, matches, gaps)
    print(f"=== Scheme: {name} ===")
    print(f"Aligned A: {a}")
    print(f"Aligned B: {b}")
    print(f"Optimal score: {score}")
    print(f"Number of matches: {matches}")
    print(f"Number of gap symbols: {gaps}")
    print()

# --- Separate check ---
# Under match-only scoring the total score equals the number of matches, so the
# optimum IS the maximum achievable match count; under the penalized scheme the
# alignment differs because gaps and mismatches now carry a cost.
mo_name = "match-only (match=1, mismatch=0, gap=0)"
pen_name = "penalized (match=3, mismatch=-3, gap=-2)"

mo_a, mo_b, mo_score, mo_matches, mo_gaps = results[mo_name]
pen_a, pen_b, pen_score, pen_matches, pen_gaps = results[pen_name]

print("=== Check ===")
print(f"Match-only score equals its match count: {mo_score} == {mo_matches} -> {mo_score == mo_matches}")
print(f"Match-only match count ({mo_matches}) >= penalized match count ({pen_matches}): {mo_matches >= pen_matches}")
print(f"Match-only gaps ({mo_gaps}) >= penalized gaps ({pen_gaps}): {mo_gaps >= pen_gaps}")
print(f"Alignments differ between schemes: {(mo_a, mo_b) != (pen_a, pen_b)}")
print("Explanation: because the match-only objective is literally the match count, its optimum "
      "must be the maximum possible matches (indifferent to gaps), and the penalized scheme's "
      "different alignment with fewer gaps confirms that gap/mismatch costs trade matches for a tighter alignment.")

# --- Figure ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, name in zip(axes, [mo_name, pen_name]):
    a, b, score, matches, gaps = results[name]
    ax.axis("off")
    match_line = "".join("|" if (x == y and x != "-") else " " for x, y in zip(a, b))
    text = f"A: {a}\n   {match_line}\nB: {b}"
    ax.text(0.5, 0.6, text, family="monospace", fontsize=14, ha="center", va="center")
    ax.set_title(f"{name}\nscore={score}, matches={matches}, gaps={gaps}", fontsize=10)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9B.1.1_s2.png")
print("\nSaved figure to 9B.1.1_s2.png")
