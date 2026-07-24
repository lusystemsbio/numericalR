import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The list of numbers to classify
numbers = [4, 7, 10, 13]

# Loop over each number and label it explicitly
labels = {}  # store results for the separate check
for n in numbers:
    # Use the modulo operator: remainder when dividing by 2
    remainder = n % 2
    if remainder == 0:
        # No remainder means the number is divisible by 2 -> even
        label = "even"
    else:
        # A remainder of 1 means it is not divisible by 2 -> odd
        label = "odd"
    labels[n] = label
    # Print the label for this number
    print(f"{n}: {label}")

# --- Separate check ---
# Confirm the expected classifications explicitly
print(f"Check 4 is even: {labels[4] == 'even'}")
print(f"Check 10 is even: {labels[10] == 'even'}")
print(f"Check 7 is odd: {labels[7] == 'odd'}")
print(f"Check 13 is odd: {labels[13] == 'odd'}")

# Overall confirmation
all_correct = (labels[4] == "even" and labels[10] == "even"
               and labels[7] == "odd" and labels[13] == "odd")
print(f"All checks pass: {all_correct}")

# Explanation (one sentence):
# This check confirms the result because a number's parity is defined by
# whether it is divisible by 2, so matching each printed label against the
# known-correct even/odd values for 4, 10, 7, and 13 verifies the modulo
# logic produced the right answer.
print("Explanation: The check confirms the result because comparing each "
      "printed label to the known parity of 4, 10 (even) and 7, 13 (odd) "
      "directly verifies that the modulo test classified every number correctly.")

# Simple visualization of the results
colors = ["tab:blue" if labels[n] == "even" else "tab:orange" for n in numbers]
plt.bar([str(n) for n in numbers], [n % 2 for n in numbers], color=colors)
plt.xlabel("Number")
plt.ylabel("n % 2 (0 = even, 1 = odd)")
plt.title("Even/Odd classification via modulo")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.3.1_s2.png")
