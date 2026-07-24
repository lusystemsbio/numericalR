import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Boolean toggle switch on n = 2 genes -------------------------------
# State = (x, y) with x, y in {0, 1}. Update rule (mutual repression):
#   x(t+1) = NOT y(t)
#   y(t+1) = NOT x(t)
def update(state):
    x, y = state
    # each gene turns ON only if its repressor is OFF
    return (1 - y, 1 - x)

# ---- Enumerate all 2^n = 4 states ---------------------------------------
n = 2
states = [(x, y) for x in (0, 1) for y in (0, 1)]

def fmt(s):
    # print a state (x, y) as the 2-bit string "xy"
    return f"{s[0]}{s[1]}"

# ---- Attractor enumeration (method from 10C.4) --------------------------
# For each start state, follow the deterministic trajectory until a state
# repeats. The set of states in the repeating loop is an attractor
# (a cycle of length 1 = fixed point, length > 1 = oscillation).
attractors = []          # list of attractor cycles (each a list of states)
seen_in_attractor = set()  # states already assigned to a found attractor

for start in states:
    trajectory = []      # ordered list of visited states this run
    seen_index = {}      # state -> position in trajectory
    s = start
    # walk forward until we revisit a state on THIS trajectory
    while s not in seen_index:
        seen_index[s] = len(trajectory)
        trajectory.append(s)
        s = update(s)
    # the cycle is the tail of the trajectory from the first repeat onward
    cycle = trajectory[seen_index[s]:]
    # canonicalize the cycle so we don't record the same loop twice
    key = frozenset(cycle)
    if key not in seen_in_attractor:
        seen_in_attractor.add(key)
        attractors.append(cycle)

# ---- Print attractors as arrow chains -----------------------------------
print("Toggle-switch attractors (arrow chains):")
fixed_points = []
oscillations = []
for cyc in attractors:
    if len(cyc) == 1:
        # fixed point: show as state -> state
        chain = f"{fmt(cyc[0])} -> {fmt(cyc[0])}"
        fixed_points.append(cyc[0])
    else:
        # oscillation: show the full loop returning to its start
        chain = " -> ".join(fmt(st) for st in cyc) + f" -> {fmt(cyc[0])}"
        oscillations.append(cyc)
    print("  " + chain)

# ---- Separate check ------------------------------------------------------
print()
print("CHECK:")
fp_set = set(fixed_points)
print("Number of fixed points:", len(fixed_points))
print("Fixed points found:", sorted(fmt(s) for s in fixed_points))
print("Contains fixed point 01:", (0, 1) in fp_set)
print("Contains fixed point 10:", (1, 0) in fp_set)

# check the two-state oscillation 00 -> 11 -> 00
osc_ok = any(frozenset(cyc) == frozenset([(0, 0), (1, 1)]) and len(cyc) == 2
             for cyc in oscillations)
print("Number of oscillations:", len(oscillations))
print("Contains oscillation {00, 11}:", osc_ok)

check_passed = (len(fixed_points) == 2 and (0, 1) in fp_set and (1, 0) in fp_set
                and osc_ok)
print("Check passed:", check_passed)

# Explanation (one sentence):
print()
print("Why this confirms the result: the enumeration recovers exactly the two "
      "mutually exclusive on-states 01 and 10 as stable fixed points plus the "
      "00<->11 oscillation, which are precisely the attractors a bistable "
      "toggle switch must have, so a correct method reproducing them confirms it.")

# ---- Visualization: state-transition graph -------------------------------
pos = {(0, 0): (0, 0), (0, 1): (0, 1), (1, 0): (1, 0), (1, 1): (1, 1)}
fig, ax = plt.subplots(figsize=(5, 5))
for s in states:
    t = update(s)
    x0, y0 = pos[s]
    x1, y1 = pos[t]
    color = "red" if s in fp_set else "blue"
    if s == t:
        ax.annotate("", xy=(x0 + 0.06, y0 + 0.06), xytext=(x0 + 0.18, y0 + 0.18),
                    arrowprops=dict(arrowstyle="->", color=color))
    else:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color=color, shrinkA=12, shrinkB=12))
for s in states:
    x0, y0 = pos[s]
    ax.plot(x0, y0, "o", ms=28, color="lightgray", zorder=3)
    ax.text(x0, y0, fmt(s), ha="center", va="center", zorder=4, fontsize=12)
ax.set_xlim(-0.4, 1.4)
ax.set_ylim(-0.4, 1.4)
ax.set_title("Boolean toggle switch: state-transition graph")
ax.axis("off")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10C.5.1_s1.png")
