import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Babylonian (Heron's) iteration for sqrt(a), a > 0.
# Model: sqrt(a) is the positive root of x^2 = a.
# Update rule: x_next = (x + a/x) / 2, iterate until |x_next^2 - a| < epsilon.

def babylonian_sqrt(a, x0, epsilon):
    x = float(x0)              # current estimate, starting from the guess
    iterates = [x]             # record successive iterates for inspection
    while True:
        x_next = (x + a / x) / 2.0        # one Babylonian averaging step
        iterates.append(x_next)           # store the new iterate
        if abs(x_next * x_next - a) < epsilon:   # stop when x_next^2 is close to a
            break
        x = x_next                        # otherwise continue iterating
    return x_next, iterates

# Test parameters
a = 2.0
guesses = [1.0, 100.0, 200.0]
epsilon = 1e-6

results = {}   # store final estimate and iterate count per guess

print(f"Target: sqrt({a}) with tolerance epsilon = {epsilon}")
print("=" * 60)

for x0 in guesses:
    final, iterates = babylonian_sqrt(a, x0, epsilon)
    results[x0] = (final, iterates)
    n_steps = len(iterates) - 1   # number of iteration steps taken
    print(f"\nStarting guess x0 = {x0}")
    for k, xk in enumerate(iterates):
        print(f"  iterate[{k}] = {xk:.10f}")
    print(f"  Final estimate       = {final:.10f}")
    print(f"  Number of iterations = {n_steps}")
    print(f"  Final estimate^2     = {final * final:.10f}")

# ---- Separate check ----
print("\n" + "=" * 60)
print("CHECK:")
reference = 1.41421   # expected value of sqrt(2) to 5 decimals
for x0 in guesses:
    final, iterates = results[x0]
    n_steps = len(iterates) - 1
    close_to_ref = abs(final - reference) < 1e-4          # converges to ~1.41421
    sq_within_tol = abs(final * final - a) < epsilon       # square within tolerance of a
    print(f"  x0 = {x0:>6}: estimate = {final:.5f}, "
          f"|estimate-1.41421| < 1e-4 = {close_to_ref}, "
          f"|estimate^2 - 2| < eps = {sq_within_tol}, "
          f"steps = {n_steps}")

steps_list = [len(results[x0][1]) - 1 for x0 in guesses]
print(f"\n  Iterations per guess (1, 100, 200): {steps_list}")
print("  A far starting guess only needs a few more iterations to converge.")

# Explanation of why the check confirms the result:
print("\n  Explanation: because x^2 = a defines sqrt(a), showing every guess "
      "lands on the same value ~1.41421 whose square is within epsilon of 2 "
      "confirms we found the true positive root regardless of the starting point.")

# ---- Plot the iterates for each starting guess ----
plt.figure(figsize=(8, 5))
for x0 in guesses:
    _, iterates = results[x0]
    plt.plot(range(len(iterates)), iterates, marker="o", label=f"x0 = {x0}")
plt.axhline(reference, color="black", linestyle="--", linewidth=1, label="sqrt(2) ~ 1.41421")
plt.xlabel("Iteration k")
plt.ylabel("Iterate value x_k")
plt.title("Babylonian iteration converging to sqrt(2)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1C.3.1_s3.png")
