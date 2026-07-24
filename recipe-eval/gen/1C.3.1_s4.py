import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Babylonian (Heron's) iteration for square root of a > 0.
# Model: sqrt(a) is the positive x solving x^2 = a.
# Iteration: x_next = (x + a/x) / 2, repeated until |x_next^2 - a| < epsilon.

def babylonian_sqrt(a, x0, epsilon):
    x = float(x0)              # current estimate, starts at the guess
    iterates = [x]             # record successive iterates, including the start
    while True:
        x_next = (x + a / x) / 2.0   # averaging step: pull estimate toward sqrt(a)
        iterates.append(x_next)      # store the new iterate
        if abs(x_next * x_next - a) < epsilon:  # stop when x_next^2 is within tol of a
            break
        x = x_next                   # otherwise iterate again
    return x_next, iterates

a = 2.0
epsilon = 1e-6
guesses = [1, 100, 200]

results = {}
for g in guesses:
    final, iterates = babylonian_sqrt(a, g, epsilon)
    results[g] = (final, iterates)
    print(f"--- Starting guess x0 = {g} ---")
    for k, xi in enumerate(iterates):
        print(f"guess {g}: iterate[{k}] = {xi:.10f}")
    print(f"guess {g}: final estimate    = {final:.10f}")
    print(f"guess {g}: number of steps   = {len(iterates) - 1}")
    print(f"guess {g}: final^2           = {final * final:.10f}")
    print(f"guess {g}: |final^2 - a|     = {abs(final * final - a):.3e}")
    print()

# Separate check: every guess should converge to ~1.41421, its square within tol of 2.
print("=== Verification check ===")
reference = 1.41421
for g in guesses:
    final, iterates = results[g]
    close_to_ref = abs(final - reference) < 1e-3
    sq_within_tol = abs(final * final - a) < epsilon
    print(f"guess {g}: final = {final:.5f}, close to {reference}? {close_to_ref}; "
          f"final^2 within tol of {a}? {sq_within_tol}; steps = {len(iterates) - 1}")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: The check confirms the result because all starting guesses reach "
      "the same value whose square is within epsilon of a, which is exactly the defining "
      "condition x^2 = a for the positive square root.")

# Plot the convergence of the iterates for each starting guess.
plt.figure(figsize=(8, 5))
for g in guesses:
    _, iterates = results[g]
    plt.plot(range(len(iterates)), iterates, marker="o", label=f"x0 = {g}")
plt.axhline(reference, color="gray", linestyle="--", label=f"sqrt(2) ≈ {reference}")
plt.yscale("log")
plt.xlabel("iteration")
plt.ylabel("iterate value (log scale)")
plt.title("Babylonian iteration converging to sqrt(2)")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1C.3.1_s4.png")
