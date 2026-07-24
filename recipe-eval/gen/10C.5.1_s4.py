import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Boolean toggle switch on n = 2 genes (X = gene 0, Y = gene 1):
#   X(t+1) = NOT Y(t)
#   Y(t+1) = NOT X(t)
# We enumerate attractors over all 2^n states (10C.4 style).
# ---------------------------------------------------------------

n = 2  # number of genes

# State representation: a tuple (X, Y) of 0/1 values.
def update(state):
    X, Y = state
    # Apply the toggle-switch update rule explicitly.
    newX = 1 - Y   # X(t+1) = NOT Y(t)
    newY = 1 - X   # Y(t+1) = NOT X(t)
    return (newX, newY)

# Build the list of all 2^n states.
all_states = [(X, Y) for X in (0, 1) for Y in (0, 1)]

def fmt(state):
    # Format a state as a bit string like "01".
    return "".join(str(b) for b in state)

# ---------------------------------------------------------------
# Attractor enumeration: from each start state, follow the
# deterministic trajectory until we revisit a state. The set of
# states in the repeating cycle is an attractor.
# ---------------------------------------------------------------
attractors = []          # list of cycles (each a list of states)
attractor_sets = []      # matching set of the states in each cycle

for start in all_states:
    seen = []            # trajectory in visitation order
    seen_index = {}      # state -> position in trajectory
    s = start
    while s not in seen_index:
        seen_index[s] = len(seen)
        seen.append(s)
        s = update(s)    # step the dynamics forward one tick
    # s is the first repeated state -> the cycle starts there.
    cycle_start = seen_index[s]
    cycle = seen[cycle_start:]
    cycle_set = frozenset(cycle)
    # Record each distinct attractor only once.
    if cycle_set not in attractor_sets:
        attractor_sets.append(cycle_set)
        attractors.append(cycle)

# ---------------------------------------------------------------
# Print the attractors as arrow chains.
# ---------------------------------------------------------------
print("Toggle-switch attractors (arrow chains):")
for cyc in attractors:
    # For a cycle, show it looping back to its first state.
    chain = " -> ".join(fmt(st) for st in cyc + [cyc[0]])
    kind = "fixed point" if len(cyc) == 1 else f"oscillation (period {len(cyc)})"
    print(f"  {chain}   [{kind}]")

# ---------------------------------------------------------------
# Separate check: confirm the expected attractors are present.
# ---------------------------------------------------------------
found_fixed = sorted(fmt(next(iter(a))) for a in attractor_sets if len(a) == 1)
found_cycles = sorted(
    tuple(sorted(fmt(s) for s in a)) for a in attractor_sets if len(a) > 1
)

expected_fixed = ["01", "10"]
expected_cycle = ("00", "11")

check_fixed = (found_fixed == expected_fixed)
check_cycle = (found_cycles == [expected_cycle])

print()
print("Number of attractors found:", len(attractors))
print("Fixed points found:", found_fixed)
print("Expected fixed points:", expected_fixed)
print("Fixed-point check passed:", check_fixed)
print("Oscillations found (as state sets):", [list(c) for c in found_cycles])
print("Expected oscillation:", list(expected_cycle), "i.e. 00 -> 11 -> 00")
print("Oscillation check passed:", check_cycle)
print("Overall check passed:", check_fixed and check_cycle)

# ---------------------------------------------------------------
# Visualize the full state-transition graph and save the figure.
# ---------------------------------------------------------------
import math
pos = {}
for i, st in enumerate(all_states):
    ang = 2 * math.pi * i / len(all_states) + math.pi / 4
    pos[st] = (math.cos(ang), math.sin(ang))

fig, ax = plt.subplots(figsize=(5, 5))
for st in all_states:
    nxt = update(st)
    x0, y0 = pos[st]
    x1, y1 = pos[nxt]
    if st == nxt:
        ax.annotate("", xy=(x0 + 0.12, y0 + 0.12), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="green", lw=2))
    else:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.5,
                                    shrinkA=15, shrinkB=15))
for st in all_states:
    x, y = pos[st]
    is_fp = update(st) == st
    ax.scatter([x], [y], s=1400,
               color=("lightgreen" if is_fp else "lightyellow"),
               edgecolor="black", zorder=3)
    ax.text(x, y, fmt(st), ha="center", va="center", fontsize=14, zorder=4)

ax.set_title("Boolean toggle switch: state transitions\n(01,10 fixed; 00<->11 oscillate)")
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.6, 1.6)
ax.axis("off")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10C.5.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print()
print("Explanation: The check confirms the result because it verifies the "
      "enumeration recovers exactly the two mutually-exclusive stable on-states "
      "(01 and 10) that make the switch bistable plus the single 00<->11 "
      "oscillation, which together account for every one of the four states.")
