import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Boolean toggle switch model on n=2 genes ---
# State is (X, Y). Update rule:
#   X(t+1) = NOT Y(t)
#   Y(t+1) = NOT X(t)
n = 2

def update(state):
    # state is a tuple (X, Y) of 0/1 ints; apply the toggle-switch rule
    X, Y = state
    Xn = 1 - Y   # NOT Y
    Yn = 1 - X   # NOT X
    return (Xn, Yn)

def fmt(state):
    # render a state tuple as a bit string, e.g. (0,1) -> "01"
    return "".join(str(b) for b in state)

# --- Enumerate all 2^n = 4 states (from 10C.4 attractor enumeration) ---
all_states = [(x, y) for x in (0, 1) for y in (0, 1)]

# --- Attractor enumeration ---
# For each start state, follow the (deterministic) trajectory until we revisit
# a state we've already seen on THIS trajectory. The cycle from the first
# repeat onward is an attractor. Deduplicate cycles by their canonical form.
def find_attractor_from(start):
    trajectory = []          # ordered list of visited states
    seen_index = {}          # state -> position in trajectory
    s = start
    while s not in seen_index:
        seen_index[s] = len(trajectory)
        trajectory.append(s)
        s = update(s)
    # s is the first repeated state; the cycle is trajectory[idx:]
    idx = seen_index[s]
    cycle = trajectory[idx:]
    return cycle

def canonical(cycle):
    # rotate the cycle so its lexicographically smallest state is first,
    # giving a unique key so the same attractor isn't counted twice
    reps = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(reps)

attractors = {}
for start in all_states:
    cyc = find_attractor_from(start)
    attractors[canonical(cyc)] = cyc

# --- Print attractors as arrow chains ---
fixed_points = []
oscillations = []
print("Toggle-switch attractors (as arrow chains):")
for key, cyc in attractors.items():
    # a cycle of length 1 is a fixed point; longer is an oscillation
    chain_states = cyc + [cyc[0]]  # close the loop for display
    chain = " -> ".join(fmt(st) for st in chain_states)
    if len(cyc) == 1:
        fixed_points.append(fmt(cyc[0]))
        print("  fixed point:", chain)
    else:
        oscillations.append([fmt(st) for st in cyc])
        print("  oscillation:", chain)

# --- Separate check: confirm expected attractors ---
found_fps = sorted(fixed_points)
print("Number of fixed points found:", len(found_fps))
print("Fixed points found:", found_fps)
check_fps = (found_fps == ["01", "10"])
print("Fixed points are exactly {01, 10}:", check_fps)

# check the oscillation 00 -> 11 -> 00
osc_ok = False
for osc in oscillations:
    if canonical([tuple(int(c) for c in s) for s in osc]) == \
       canonical([(0, 0), (1, 1)]):
        osc_ok = True
print("Number of oscillations found:", len(oscillations))
print("Oscillation 00 <-> 11 present:", osc_ok)
print("Overall check passed:", check_fps and osc_ok and len(oscillations) == 1)

# One-sentence explanation:
# WHY THIS CHECK CONFIRMS THE RESULT: recovering exactly the two fixed points
# 01 and 10 plus the single 00<->11 oscillation accounts for all 4 states with
# the biologically expected bistable-switch behavior, so the enumeration is correct.
print("Explanation: The check confirms the result because recovering exactly the "
      "two mutually exclusive fixed points 01 and 10 plus the single 00->11->00 "
      "oscillation accounts for all four states with the expected bistable-switch "
      "dynamics, verifying the enumeration is correct.")

# --- Visualization: state-transition graph ---
import math
coords = {s: (math.cos(2*math.pi*i/4), math.sin(2*math.pi*i/4))
          for i, s in enumerate(all_states)}
fig, ax = plt.subplots(figsize=(6, 6))
for s in all_states:
    t = update(s)
    x0, y0 = coords[s]
    x1, y1 = coords[t]
    if s == t:
        ax.annotate("", xy=(x0, y0+0.12), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="tab:red"))
    else:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="tab:blue",
                                    shrinkA=15, shrinkB=15))
for s, (x, y) in coords.items():
    color = "tab:green" if fmt(s) in found_fps else "lightgray"
    ax.scatter([x], [y], s=1200, color=color, zorder=3, edgecolors="black")
    ax.text(x, y, fmt(s), ha="center", va="center", fontsize=14, zorder=4)
ax.set_title("Boolean toggle switch: state-transition graph")
ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6)
ax.axis("off")
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10C.5.1_s2.png")
