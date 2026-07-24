import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1C.3.1_s5.png"


def babylonian_sqrt(a, x0, epsilon=1e-6, max_iter=1000):
    # Babylonian iteration: repeatedly average x with a/x.
    # Store every iterate so we can inspect convergence.
    x = float(x0)          # start from the positive guess
    iterates = [x]         # record the initial guess
    for _ in range(max_iter):
        x_next = (x + a / x) / 2.0   # Babylonian update step
        iterates.append(x_next)      # keep the new iterate
        if abs(x_next * x_next - a) < epsilon:  # stop when x^2 is close to a
            break
        x = x_next                   # advance for the next step
    return iterates


a = 2.0
epsilon = 1e-6
guesses = [1, 100, 200]

results = {}  # store iterates per starting guess

for x0 in guesses:
    iterates = babylonian_sqrt(a, x0, epsilon)
    results[x0] = iterates
    final = iterates[-1]
    print(f"--- Starting guess x0 = {x0} ---")
    for i, xi in enumerate(iterates):
        print(f"  iterate[{i}] = {xi:.10f}")
    print(f"Final estimate (x0={x0}): {final:.10f}")
    print(f"Iterations needed (x0={x0}): {len(iterates) - 1}")
    print()

# --- Separate check: verify convergence, accuracy, and iteration cost ---
print("=== Verification check ===")
all_near_root = True
all_within_tol = True
for x0 in guesses:
    final = results[x0][-1]
    residual = abs(final * final - a)
    near_root = abs(final - 1.41421) < 1e-3
    within_tol = residual < epsilon
    all_near_root = all_near_root and near_root
    all_within_tol = all_within_tol and within_tol
    print(f"x0={x0}: final={final:.6f}, final^2={final*final:.10f}, "
          f"|final^2 - a|={residual:.3e}, converges~1.41421={near_root}, "
          f"squared-within-tol={within_tol}")

iters = {x0: len(results[x0]) - 1 for x0 in guesses}
print(f"Iteration counts by guess: {iters}")
print(f"All guesses converge to ~1.41421: {all_near_root}")
print(f"All squared results within tolerance of 2: {all_within_tol}")
print(f"Far guess (200) needs only a few more iterations than near guess (1): "
      f"{iters[200] - iters[1]} extra")

# Explanation (one sentence):
print("Explanation: The check confirms the result because the stopping "
      "criterion guarantees |x^2 - a| < epsilon, so x^2 is provably within "
      "tolerance of 2 and x therefore approximates the true positive root "
      "sqrt(2) ~ 1.41421 regardless of the starting guess.")

# --- Plot the convergence of iterates for each starting guess ---
plt.figure(figsize=(8, 5))
for x0 in guesses:
    it = results[x0]
    plt.plot(range(len(it)), it, marker="o", label=f"start = {x0}")
plt.axhline(1.4142135624, color="gray", linestyle="--", label="sqrt(2)")
plt.xlabel("iteration index")
plt.ylabel("iterate value")
plt.yscale("log")
plt.title("Babylonian iteration for sqrt(2) from different guesses")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig(OUT)
print(f"Figure saved to: {OUT}")
