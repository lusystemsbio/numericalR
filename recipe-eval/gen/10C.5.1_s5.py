import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Boolean toggle switch on n = 2 genes (X, Y).
# Update rule:  X(t+1) = NOT Y(t),  Y(t+1) = NOT X(t).
# We enumerate attractors explicitly (method of 10C.4):
#   1. build the full state-transition map over all 2^n states,
#   2. iterate the map from each state until a state repeats,
#   3. the repeating cycle is the attractor.
# ---------------------------------------------------------------

n = 2  # number of genes

# --- state helpers: a state is a tuple (X, Y) of 0/1 bits ---
def to_str(state):
    # render a state tuple as a bit string, e.g. (0,1) -> "01"
    return "".join(str(b) for b in state)

# --- the explicit toggle-switch update rule ---
def update(state):
    X, Y = state
    newX = 1 - Y   # NOT Y
    newY = 1 - X   # NOT X
    return (newX, newY)

# --- enumerate all 2^n states ---
all_states = [(x, y) for x in (0, 1) for y in (0, 1)]

# --- build the transition map: state -> next state ---
transition = {s: update(s) for s in all_states}

print("Transition map (all four states):")
for s in all_states:
    print(f"  {to_str(s)} -> {to_str(transition[s])}")

# ---------------------------------------------------------------
# Attractor enumeration: from each state, follow transitions,
# recording the visited order. When we revisit a state, the
# portion of the trajectory from that state onward is the cycle.
# ---------------------------------------------------------------
attractors = []          # list of cycles (each a list of state tuples)
seen_attractor_sets = [] # frozensets used to dedupe cycles

for start in all_states:
    trajectory = []          # ordered visited states
    positions = {}           # state -> index in trajectory
    cur = start
    # walk until we hit a state we've already visited on this walk
    while cur not in positions:
        positions[cur] = len(trajectory)
        trajectory.append(cur)
        cur = transition[cur]
    # the cycle starts where we re-entered the trajectory
    cycle_start = positions[cur]
    cycle = trajectory[cycle_start:]
    # dedupe: same cycle reached from different starts is one attractor
    key = frozenset(cycle)
    if key not in seen_attractor_sets:
        seen_attractor_sets.append(key)
        attractors.append(cycle)

# ---------------------------------------------------------------
# Report the attractors as arrow chains.
# A fixed point loops back to itself; an oscillation lists its cycle.
# ---------------------------------------------------------------
fixed_points = []
oscillations = []

print("\nToggle-switch attractors (as arrow chains):")
for cycle in attractors:
    # build the arrow chain, closing the loop back to the first state
    chain = " -> ".join(to_str(s) for s in cycle) + " -> " + to_str(cycle[0])
    if len(cycle) == 1:
        kind = "fixed point"
        fixed_points.append(to_str(cycle[0]))
    else:
        kind = f"oscillation (period {len(cycle)})"
        oscillations.append([to_str(s) for s in cycle])
    print(f"  [{kind}] {chain}")

# ---------------------------------------------------------------
# Separate check: confirm two fixed points 01 and 10, plus the
# oscillation 00 -> 11 -> 00.
# ---------------------------------------------------------------
print("\nCheck results:")
print(f"Number of fixed points found: {len(fixed_points)}")
print(f"Fixed points found: {sorted(fixed_points)}")
check_fp = (sorted(fixed_points) == ["01", "10"])
print(f"Fixed points are exactly {{01, 10}}: {check_fp}")

print(f"Number of oscillations found: {len(oscillations)}")
osc_sets = [set(o) for o in oscillations]
check_osc = ({"00", "11"} in osc_sets)
print(f"Oscillation {{00, 11}} (00 -> 11 -> 00) present: {check_osc}")

print(f"Overall check passed: {check_fp and check_osc}")

# ---------------------------------------------------------------
# Visualization: draw the state-transition graph with attractors.
# ---------------------------------------------------------------
import math

fig, ax = plt.subplots(figsize=(6, 6))
# place the four states on a circle
coords = {}
for i, s in enumerate(all_states):
    theta = math.pi / 2 - 2 * math.pi * i / len(all_states)
    coords[s] = (math.cos(theta), math.sin(theta))

# draw transition arrows
for s in all_states:
    x0, y0 = coords[s]
    x1, y1 = coords[transition[s]]
    if s == transition[s]:
        # self-loop for a fixed point
        ax.annotate("", xy=(x0 * 1.18, y0 * 1.18), xytext=(x0 * 1.02, y0 * 1.02),
                    arrowprops=dict(arrowstyle="->", color="green", lw=2,
                                    connectionstyle="arc3,rad=1.5"))
    else:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="royalblue", lw=2,
                                    shrinkA=15, shrinkB=15,
                                    connectionstyle="arc3,rad=0.15"))

# draw state nodes
for s in all_states:
    x, y = coords[s]
    is_fp = (s == transition[s])
    ax.scatter([x], [y], s=1400,
               color=("mediumseagreen" if is_fp else "lightsteelblue"),
               edgecolors="black", zorder=3)
    ax.text(x, y, to_str(s), ha="center", va="center",
            fontsize=14, fontweight="bold", zorder=4)

ax.set_title("Boolean toggle switch: state-transition graph\n"
             "green = fixed points (01, 10), blue = oscillation (00<->11)")
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.6, 1.6)
ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/10C.5.1_s5.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("\nExplanation:")
print("The check confirms the result because recovering exactly the two "
      "mutually exclusive fixed points 01 and 10 plus the 00<->11 oscillation "
      "matches the known bistable-plus-antiphase behavior of the toggle switch, "
      "so the enumeration correctly captured every attractor of the 2^2 state space.")
