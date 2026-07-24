import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The list to test
numbers = [4, 7, 10, 13]

# Store labels so we can run the separate check afterward
labels = {}

# Loop over each number in the list
for n in numbers:
    # Use the modulo operator: remainder when dividing by 2
    # If the remainder is 0, the number is divisible by 2 -> even
    if n % 2 == 0:
        label = "even"
    else:
        # A non-zero remainder means it is not divisible by 2 -> odd
        label = "odd"
    labels[n] = label
    # Print the number with its computed label
    print(f"{n}: {label}")

# Separate check: confirm the expected labels explicitly
print("Check 4 is even:", labels[4] == "even")
print("Check 10 is even:", labels[10] == "even")
print("Check 7 is odd:", labels[7] == "odd")
print("Check 13 is odd:", labels[13] == "odd")

all_correct = (labels[4] == "even" and labels[10] == "even"
               and labels[7] == "odd" and labels[13] == "odd")
print("All labels correct:", all_correct)

# Explanation:
# This check confirms the result because 4 and 10 are known to be divisible by 2
# (even) while 7 and 13 are not (odd), so matching these known values verifies
# the modulo logic assigned the correct label to every number.

# Simple bar chart visualizing the even/odd labeling (0 = odd, 1 = even)
codes = [1 if labels[n] == "even" else 0 for n in numbers]
colors = ["tab:blue" if c == 1 else "tab:orange" for c in codes]
plt.bar([str(n) for n in numbers], codes, color=colors)
plt.yticks([0, 1], ["odd", "even"])
plt.xlabel("Number")
plt.ylabel("Label")
plt.title("Even/Odd labeling of [4, 7, 10, 13]")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.3.1_s1.png")
