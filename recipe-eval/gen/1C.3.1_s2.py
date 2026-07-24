import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Babylonian (Newton) iteration for sqrt(a): the positive root of x^2 = a.
def babylonian_sqrt(a, x0, epsilon=1e-6):
    x = float(x0)          # start from the given positive guess
    iterates = [x]         # record the successive iterates
    # repeat until the squared estimate is within epsilon of a
    while abs(x * x - a) >= epsilon:
        x = (x + a / x) / 2.0   # averaging step: pull x toward sqrt(a)
        iterates.append(x)
    return iterates

a = 2.0
epsilon = 1e-6
guesses = [1, 100, 200]

results = {}
for g in guesses:
    iters = babylonian_sqrt(a, g, epsilon)
    results[g] = iters
    print(f"--- Starting guess x0 = {g} ---")
    for k, xk in enumerate(iters):
        print(f"iterate[{k}] = {xk:.15f}")
    final = iters[-1]
    print(f"final estimate (x0={g}) = {final:.15f}")
    print(f"iterations needed (x0={g}) = {len(iters) - 1}")
    print(f"final^2 (x0={g}) = {final * final:.15f}")
    print(f"|final^2 - a| (x0={g}) = {abs(final * final - a):.3e}")

# --- Separate check ---
reference = 1.41421
print("--- Check ---")
all_close_to_ref = True
all_within_tol = True
for g in guesses:
    final = results[g][-1]
    close_ref = abs(final - reference) < 1e-3
    within_tol = abs(final * final - a) < epsilon
    all_close_to_ref = all_close_to_ref and close_ref
    all_within_tol = all_within_tol and within_tol
    print(f"guess {g}: final={final:.5f}, close to 1.41421? {close_ref}, "
          f"square within tol of 2? {within_tol}")

print(f"all guesses converge to ~1.41421 : {all_close_to_ref}")
print(f"all final squares within tolerance of 2 : {all_within_tol}")

iters_needed = {g: len(results[g]) - 1 for g in guesses}
print(f"iterations per guess : {iters_needed}")
far_needs_a_few_more = iters_needed[200] > iters_needed[1]
print(f"far guess (200) needs a few more iterations than near guess (1) : "
      f"{far_needs_a_few_more}")
print(f"extra iterations for far guess (200 vs 1) : "
      f"{iters_needed[200] - iters_needed[1]}")

# Explanation of why the check confirms the result:
print("Explanation: The check confirms the result because a value whose square "
      "is within epsilon of 2 must be within tolerance of the true positive "
      "root sqrt(2)=1.41421..., and seeing every starting guess reach it "
      "(the far one only needing a few extra steps) shows convergence is "
      "correct and independent of the initial guess.")

# --- Plot the iterates for each starting guess ---
plt.figure(figsize=(8, 5))
for g in guesses:
    iters = results[g]
    plt.plot(range(len(iters)), iters, marker="o", label=f"x0 = {g}")
plt.axhline(reference, color="k", linestyle="--", linewidth=1,
            label="sqrt(2) ≈ 1.41421")
plt.xlabel("iteration k")
plt.ylabel("iterate x_k")
plt.title("Babylonian iteration converging to sqrt(2)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1C.3.1_s2.png")
