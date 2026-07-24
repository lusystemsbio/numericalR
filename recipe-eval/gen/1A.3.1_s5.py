import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The list of numbers to classify
numbers = [4, 7, 10, 13]

# Store labels so we can both print them and run a separate check afterward
labels = {}

# Loop over each number and decide even vs odd explicitly
for n in numbers:
    remainder = n % 2          # remainder when divided by 2
    if remainder == 0:         # no remainder means the number is even
        label = "even"
    else:                      # any remainder means the number is odd
        label = "odd"
    labels[n] = label
    print(f"{n} is {label}")   # print the label for this number

# --- Separate check ---------------------------------------------------------
# Confirm the expected classifications hold.
check_even = (labels[4] == "even") and (labels[10] == "even")
check_odd = (labels[7] == "odd") and (labels[13] == "odd")
print(f"Check 4 and 10 are even: {check_even}")
print(f"Check 7 and 13 are odd: {check_odd}")
print(f"All checks passed: {check_even and check_odd}")

# Explanation (one sentence):
# This check confirms the result because a number is even exactly when it is
# divisible by 2 (remainder 0) and odd otherwise, so verifying that 4 and 10
# yield "even" while 7 and 13 yield "odd" matches the known parity of these values.
explanation = ("The check confirms the result because parity is defined by "
               "divisibility by 2, and 4 and 10 (remainder 0) are indeed even "
               "while 7 and 13 (remainder 1) are indeed odd.")
print(explanation)

# Produce a simple figure visualizing the classification
fig, ax = plt.subplots()
colors = ["tab:blue" if labels[n] == "even" else "tab:orange" for n in numbers]
ax.bar([str(n) for n in numbers], [n % 2 for n in numbers], color=colors)
ax.set_xlabel("Number")
ax.set_ylabel("Remainder mod 2 (0 = even, 1 = odd)")
ax.set_title("Even/Odd classification")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.3.1_s5.png")
