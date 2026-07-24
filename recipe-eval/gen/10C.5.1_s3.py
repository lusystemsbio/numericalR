import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Boolean toggle switch on n = 2 genes ---
# State is (X, Y). Update rule: X(t+1) = NOT Y(t), Y(t+1) = NOT X(t).
def update(state):
    x, y = state
    return (1 - y, 1 - x)   # NOT Y, NOT X

n = 2
# Enumerate all 2^n = 4 states as tuples.
states = [(x, y) for x in (0, 1) for y in (0, 1)]

def fmt(s):
    # Format a state tuple as a bit string, e.g. (0,1) -> "01".
    return "".join(str(b) for b in s)

# --- Boolean attractor enumeration (10C.4) ---
# From each start state, follow the deterministic trajectory until a state
# repeats. The repeated state marks entry into the attractor cycle; the
# states between its first and second occurrence form the attractor.
attractors = []          # list of attractors, each a list of states (the cycle)
seen_cycles = set()      # canonical (frozenset) form to avoid duplicates

for start in states:
    trajectory = []
    s = start
    # Walk forward, recording states, until we hit one already in the path.
    while s not in trajectory:
        trajectory.append(s)
        s = update(s)
    # s is the first repeated state: extract the cycle from its first occurrence.
    idx = trajectory.index(s)
    cycle = trajectory[idx:]
    key = frozenset(cycle)   # cycles are the same regardless of start phase
    if key not in seen_cycles:
        seen_cycles.add(key)
        attractors.append(cycle)

# --- Print attractors as arrow chains (closing the loop back to start) ---
print("Toggle-switch attractors:")
fixed_points = []
oscillations = []
for cyc in attractors:
    chain = " -> ".join(fmt(s) for s in cyc + [cyc[0]])
    print(chain)
    if len(cyc) == 1:
        fixed_points.append(fmt(cyc[0]))
    else:
        oscillations.append([fmt(s) for s in cyc])

# --- Separate check: two fixed points 01 and 10, one oscillation 00 <-> 11 ---
print()
print("Number of fixed points:", len(fixed_points))
print("Fixed points found:", sorted(fixed_points))
print("Number of oscillations:", len(oscillations))
for osc in oscillations:
    print("Oscillation cycle:", " -> ".join(osc + [osc[0]]))

check_fixed = sorted(fixed_points) == ["01", "10"]
# The oscillation set {00, 11} identifies the 2-cycle regardless of phase.
check_osc = any(set(osc) == {"00", "11"} for osc in oscillations)
print("Check fixed points == {01, 10}:", check_fixed)
print("Check oscillation == {00 <-> 11}:", check_osc)
print("All checks pass:", check_fixed and check_osc)

# One-sentence explanation of why the check confirms the result:
print("Explanation: because enumeration over all four states must partition "
      "them into exactly these attractors, recovering the two mutually "
      "exclusive on-states (01, 10) as the bistable fixed points and the "
      "00<->11 2-cycle confirms the switch's known behavior.")

# --- Visualization: state transition graph on the four states ---
pos = {(0, 0): (0, 0), (0, 1): (1, 1), (1, 0): (1, -1), (1, 1): (2, 0)}
fig, ax = plt.subplots(figsize=(6, 4))
for s in states:
    t = update(s)
    x0, y0 = pos[s]
    x1, y1 = pos[t]
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color="steelblue", lw=1.5))
for s in states:
    x0, y0 = pos[s]
    is_fp = fmt(s) in fixed_points
    ax.scatter([x0], [y0], s=900,
               color="lightgreen" if is_fp else "lightsalmon",
               edgecolors="black", zorder=3)
    ax.text(x0, y0, fmt(s), ha="center", va="center", fontweight="bold", zorder=4)
ax.set_title("Boolean toggle switch: state transition graph")
ax.axis("off")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10C.5.1_s3.png")
