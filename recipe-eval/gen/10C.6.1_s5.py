import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Boolean repressilator update rule on n = 3 genes ----
# State is a tuple (X, Y, Z) of 0/1 values.
# Ring of three repressions:
#   X(t+1) = NOT Z(t)
#   Y(t+1) = NOT X(t)
#   Z(t+1) = NOT Y(t)
def update(state):
    X, Y, Z = state
    return (1 - Z, 1 - X, 1 - Y)


# ---- Boolean attractor enumeration (from 10C.4) over all 2^n states ----
# Because the state space is finite, iterating the deterministic update from
# any state must eventually enter a cycle (the attractor). We follow the
# trajectory from each of the 8 states, detect the first repeated state,
# and extract the cycle it settled into.
n = 3
all_states = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]


def trajectory_cycle(start):
    # Walk forward, remembering the order states were first seen.
    seen = {}            # state -> position in the walk
    order = []           # states in the order visited
    s = start
    while s not in seen:
        seen[s] = len(order)
        order.append(s)
        s = update(s)
    # s is the first state that repeats; the cycle is the tail from its
    # first occurrence to the end of the recorded order.
    cycle = order[seen[s]:]
    return cycle


# Collect the distinct attractors (cycles) reachable from every state.
attractors = []                      # list of cycles (each a list of states)
canonical_seen = set()               # canonical keys to dedupe cycles


def canonical(cycle):
    # A cycle is the same regardless of starting point; use the rotation
    # that begins at its lexicographically smallest state as an identity key.
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(rotations)


for start in all_states:
    cyc = trajectory_cycle(start)
    key = canonical(cyc)
    if key not in canonical_seen:
        canonical_seen.add(key)
        attractors.append(list(key))


# ---- Format helpers ----
def bits(state):
    return "".join(str(b) for b in state)


def chain(cycle):
    # Print an attractor as an arrow chain that returns to its start.
    return " -> ".join(bits(s) for s in cycle) + " -> " + bits(cycle[0])


# ---- Report the attractors as arrow chains ----
print("Boolean repressilator attractors (n = 3, over all 8 states):")
for i, cyc in enumerate(attractors, 1):
    print("Attractor %d (period %d): %s" % (i, len(cyc), chain(cyc)))

# ---- Separate check: no fixed points, exactly two cyclic attractors ----
fixed_points = [s for s in all_states if update(s) == s]
periods = sorted(len(c) for c in attractors)

print("")
print("CHECK - number of fixed points found:", len(fixed_points))
print("CHECK - number of attractors found:", len(attractors))
print("CHECK - attractor periods (sorted):", periods)
print("CHECK - fixed points (period-1 attractors):", fixed_points if fixed_points else "none")

# Confirm the expected structure explicitly.
has_no_fixed_points = (len(fixed_points) == 0)
has_two_cycles = (len(attractors) == 2)
expected_periods = (periods == [2, 6])
print("CHECK - no fixed points:", has_no_fixed_points)
print("CHECK - exactly two cyclic attractors:", has_two_cycles)
print("CHECK - periods are {2, 6}:", expected_periods)
print("CHECK - all conditions satisfied:", has_no_fixed_points and has_two_cycles and expected_periods)

# Why this check confirms the result, in one sentence:
print("")
print("Why: finding zero period-1 states rules out any steady state, so the")
print("only long-term behaviors are the two cycles - a period-2 (000<->111)")
print("and a period-6 orbit, the discrete analog of the continuous limit cycle.")

# ---- Visualize the two attractors as directed cycles ----
import math
fig, axes = plt.subplots(1, len(attractors), figsize=(6 * len(attractors), 6))
if len(attractors) == 1:
    axes = [axes]

for ax, cyc in zip(axes, attractors):
    m = len(cyc)
    # Place states evenly around a circle.
    pts = [(math.cos(2 * math.pi * k / m + math.pi / 2),
            math.sin(2 * math.pi * k / m + math.pi / 2)) for k in range(m)]
    for k in range(m):
        x0, y0 = pts[k]
        x1, y1 = pts[(k + 1) % m]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color="tab:blue",
                                    shrinkA=18, shrinkB=18, lw=1.5))
    for (px, py), s in zip(pts, cyc):
        ax.scatter([px], [py], s=900, c="white", edgecolors="black", zorder=3)
        ax.text(px, py, bits(s), ha="center", va="center",
                fontsize=12, zorder=4)
    ax.set_title("Period-%d attractor" % m)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")

fig.suptitle("Boolean repressilator attractors", fontsize=14)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.6.1_s5.png")
print("")
print("Figure saved.")
