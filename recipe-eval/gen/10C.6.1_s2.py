import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import product

# ------------------------------------------------------------------
# Repressilator update rule on n = 3 genes (ring of three repressions)
#   X(t+1) = NOT Z(t)
#   Y(t+1) = NOT X(t)
#   Z(t+1) = NOT Y(t)
# State encoded as a tuple (x, y, z) of 0/1 values.
# ------------------------------------------------------------------
def update(state):
    x, y, z = state
    return (1 - z, 1 - x, 1 - y)  # NOT is 1 - value for Boolean 0/1

# ------------------------------------------------------------------
# Boolean attractor enumeration (from 10C.4):
# Because the state space is finite (2^n states) and the dynamics are
# deterministic, every trajectory must eventually revisit a state and
# fall into a cycle (the attractor). We enumerate all states, follow
# each trajectory until it repeats, and record the cycle it lands in.
# ------------------------------------------------------------------
n = 3
all_states = list(product([0, 1], repeat=n))   # all eight states

def fmt(s):
    return "".join(str(b) for b in s)          # tuple -> "010" string

attractors = []          # list of cycles, each a list of states
seen_in_attractor = set()

for start in all_states:
    # Walk the trajectory, remembering the order states were first seen.
    path = []
    pos = {}             # state -> index along this path
    s = start
    while s not in pos:
        pos[s] = len(path)
        path.append(s)
        s = update(s)
    # s is the first repeated state: the cycle runs from pos[s] to the end.
    cycle = path[pos[s]:]
    # Canonicalize the cycle (rotate to its lexicographically smallest
    # start) so we can deduplicate cycles found from different starts.
    k = min(range(len(cycle)), key=lambda i: cycle[i])
    canon = tuple(cycle[k:] + cycle[:k])
    if canon not in seen_in_attractor:
        seen_in_attractor.add(canon)
        attractors.append(list(canon))

# ------------------------------------------------------------------
# Report attractors as arrow chains.
# ------------------------------------------------------------------
print("Number of attractors found:", len(attractors))
fixed_points = []
cyclic = []
for cyc in attractors:
    # Arrow chain closes back on the first state to show it is a cycle.
    chain = " -> ".join(fmt(s) for s in cyc) + " -> " + fmt(cyc[0])
    length = len(cyc)
    print(f"Attractor (period {length}): {chain}")
    if length == 1:
        fixed_points.append(cyc)
    else:
        cyclic.append(cyc)

# ------------------------------------------------------------------
# Separate check: no fixed points, exactly two cyclic attractors,
# one of period 2 (000 -> 111 -> 000) and one of period 6.
# ------------------------------------------------------------------
print("Number of fixed points:", len(fixed_points))
print("Number of cyclic attractors:", len(cyclic))
periods = sorted(len(c) for c in cyclic)
print("Cyclic attractor periods:", periods)

check_no_fixed = (len(fixed_points) == 0)
check_two_cycles = (len(cyclic) == 2)
check_periods = (periods == [2, 6])
print("Check - no fixed points:", check_no_fixed)
print("Check - exactly two cyclic attractors:", check_two_cycles)
print("Check - periods are [2, 6]:", check_periods)
print("All checks pass:", check_no_fixed and check_two_cycles and check_periods)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Finding zero fixed points and only cyclic attractors "
      "(a period-2 and a period-6 cycle) confirms the repressilator has no "
      "stable steady state and instead settles into sustained oscillation, "
      "the discrete analog of the continuous limit cycle.")

# ------------------------------------------------------------------
# Visualization: draw each attractor as a ring of its states.
# ------------------------------------------------------------------
import math
fig, axes = plt.subplots(1, len(attractors), figsize=(5 * len(attractors), 5))
if len(attractors) == 1:
    axes = [axes]
for ax, cyc in zip(axes, attractors):
    m = len(cyc)
    xs = [math.cos(2 * math.pi * i / m) for i in range(m)]
    ys = [math.sin(2 * math.pi * i / m) for i in range(m)]
    for i in range(m):
        j = (i + 1) % m
        ax.annotate("", xy=(xs[j], ys[j]), xytext=(xs[i], ys[i]),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.scatter(xs, ys, s=600, c="lightblue", zorder=3, edgecolors="k")
    for i, s in enumerate(cyc):
        ax.text(xs[i], ys[i], fmt(s), ha="center", va="center", zorder=4)
    ax.set_title(f"Attractor (period {m})")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.6.1_s2.png")
