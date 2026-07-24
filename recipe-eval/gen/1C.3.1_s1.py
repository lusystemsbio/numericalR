import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Babylonian square-root iteration (no built-in sqrt used).
def babylonian_sqrt(a, x0, epsilon=1e-6):
    # Start from the positive guess x0 and keep the successive iterates.
    iterates = [x0]
    x = x0
    while True:
        # Babylonian update: average of x and a/x.
        x_next = (x + a / x) / 2.0
        iterates.append(x_next)
        # Stop once x_next^2 is within epsilon of a.
        if abs(x_next * x_next - a) < epsilon:
            break
        x = x_next
    return iterates

a = 2.0
epsilon = 1e-6
guesses = [1, 100, 200]

results = {}  # store (iterates, final estimate, iteration count) per guess

for g in guesses:
    iterates = babylonian_sqrt(a, float(g), epsilon)
    final = iterates[-1]
    n_iters = len(iterates) - 1  # number of update steps performed
    results[g] = (iterates, final, n_iters)

    print(f"=== Starting guess x0 = {g} ===")
    for i, xi in enumerate(iterates):
        print(f"  iterate[{i}] = {xi:.15f}")
    print(f"  final estimate       = {final:.15f}")
    print(f"  iterations to converge = {n_iters}")
    print(f"  final^2              = {final * final:.15f}")
    print(f"  |final^2 - a|        = {abs(final * final - a):.3e}")
    print()

# --- Separate check --------------------------------------------------------
print("=== Convergence check ===")
target = 1.41421  # approximate expected value of sqrt(2)
for g in guesses:
    iterates, final, n_iters = results[g]
    close_to_target = abs(final - target) < 1e-3        # converges to ~1.41421
    squared_in_tol = abs(final * final - a) < epsilon    # final^2 within tol of 2
    print(f"guess {g:>4}: final={final:.15f}, "
          f"|final-1.41421|<1e-3 -> {close_to_target}, "
          f"|final^2-2|<eps -> {squared_in_tol}, "
          f"iters={n_iters}")

far = results[200][2]
near = results[1][2]
print(f"far guess (200) used {far} iterations vs {near} for guess 1; "
      f"extra iterations = {far - near}")

# Explanation of why the check confirms the result:
print("Explanation: Since squaring the final estimate reproduces a to within "
      "epsilon, the estimate must be the positive root of x^2=a, so it is a "
      "genuine square root regardless of the (possibly far) starting guess.")

# --- Plot: iterate value vs iteration index for each starting guess --------
plt.figure(figsize=(8, 5))
for g in guesses:
    iterates, final, n_iters = results[g]
    plt.plot(range(len(iterates)), iterates, marker="o", label=f"x0 = {g}")
plt.axhline(target, color="k", linestyle="--", linewidth=1, label="~1.41421")
plt.xlabel("iteration index")
plt.ylabel("iterate value")
plt.title("Babylonian iteration converging to sqrt(2)")
plt.yscale("log")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1C.3.1_s1.png")
