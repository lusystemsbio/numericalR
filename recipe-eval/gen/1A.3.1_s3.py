# Label each number in a list as "even" or "odd" using explicit modulo testing.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The list to test
numbers = [4, 7, 10, 13]

# Store labels so we can do a separate check afterward
labels = {}

# Loop over each number and determine its label explicitly
for n in numbers:
    # The modulo operator gives the remainder after division by 2.
    # A remainder of 0 means the number is exactly divisible by 2 -> even.
    if n % 2 == 0:
        label = "even"
    else:
        # Any non-zero remainder (which is 1 for integers) means odd.
        label = "odd"
    labels[n] = label
    print(f"{n}: {label}")

# Separate check: confirm the expected labels
print(f"Check 4 is even: {labels[4] == 'even'}")
print(f"Check 10 is even: {labels[10] == 'even'}")
print(f"Check 7 is odd: {labels[7] == 'odd'}")
print(f"Check 13 is odd: {labels[13] == 'odd'}")

# Overall pass/fail of the check
all_correct = (labels[4] == "even" and labels[10] == "even"
               and labels[7] == "odd" and labels[13] == "odd")
print(f"All labels match expected values: {all_correct}")

# Explanation (one sentence):
# This check confirms the result because 4 and 10 have remainder 0 when divided
# by 2 (so they must be even) while 7 and 13 have remainder 1 (so they must be
# odd), matching the definition of even and odd exactly.
explanation = ("This check confirms the result because it verifies each number's "
               "computed label against its known parity, so matching all four "
               "means the modulo test correctly distinguishes even from odd.")
print(f"Explanation: {explanation}")

# Simple visualization of the labels
fig, ax = plt.subplots()
colors = ["tab:blue" if labels[n] == "even" else "tab:orange" for n in numbers]
ax.bar([str(n) for n in numbers], [n % 2 for n in numbers], color=colors)
ax.set_xlabel("Number")
ax.set_ylabel("n % 2  (0 = even, 1 = odd)")
ax.set_title("Even/Odd via modulo operator")
for i, n in enumerate(numbers):
    ax.text(i, (n % 2) + 0.02, labels[n], ha="center")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.3.1_s3.png")
