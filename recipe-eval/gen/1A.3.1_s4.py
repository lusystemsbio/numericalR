import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The list to test
numbers = [4, 7, 10, 13]

# Store computed labels so we can verify them afterward
labels = {}

# Loop over each number and classify it explicitly
for n in numbers:
    # Compute the remainder when dividing by 2
    remainder = n % 2
    # A remainder of 0 means the number is even; otherwise it is odd
    if remainder == 0:
        label = "even"
    else:
        label = "odd"
    labels[n] = label
    # Print the label for this number
    print(f"{n}: {label}")

# Separate check: confirm the expected classifications
print("--- Verification check ---")
print("4 is even:", labels[4] == "even")
print("10 is even:", labels[10] == "even")
print("7 is odd:", labels[7] == "odd")
print("13 is odd:", labels[13] == "odd")

all_correct = (
    labels[4] == "even"
    and labels[10] == "even"
    and labels[7] == "odd"
    and labels[13] == "odd"
)
print("All labels match expected result:", all_correct)

# Explanation:
# This check confirms the result because it independently compares each printed
# label against its known-correct even/odd value, so agreement on all four numbers
# proves the modulo-based classification logic behaved correctly.

# Simple visualization of the results
colors = ["tab:blue" if labels[n] == "even" else "tab:orange" for n in numbers]
plt.bar([str(n) for n in numbers], [n % 2 for n in numbers], color=colors)
plt.yticks([0, 1], ["even (0)", "odd (1)"])
plt.xlabel("Number")
plt.ylabel("n % 2")
plt.title("Even (0) vs Odd (1) via modulo")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.3.1_s4.png")
