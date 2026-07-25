import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Boolean repressilator: a ring of three repressions on n = 3 genes.
#   X(t+1) = NOT Z(t)
#   Y(t+1) = NOT X(t)
#   Z(t+1) = NOT Y(t)
# State is a tuple (x, y, z) of bits. We enumerate all 2^n states,
# follow each deterministic trajectory until it re-enters a visited
# state, and report the cycle it settles into (the attractor).
# ----------------------------------------------------------------------

n = 3


def update(state):
    """One synchronous Boolean update of the repressilator."""
    x, y, z = state
    nx = 1 - z          # X(t+1) = NOT Z
    ny = 1 - x          # Y(t+1) = NOT X
    nz = 1 - y          # Z(t+1) = NOT Y
    return (nx, ny, nz)


def fmt(state):
    """Format a state tuple as a bit string like '010'."""
    return "".join(str(b) for b in state)


# ----------------------------------------------------------------------
# Attractor enumeration (method from 10C.4), done explicitly.
#
# For each of the 2^n starting states we walk the trajectory, recording
# the order in which states are first seen. The first time we hit a
# state we have already seen, everything from that repeat onward forms
# the attractor cycle. We canonicalize each cycle (rotate so its
# lexicographically smallest state comes first) so that different phases
# of the same cycle are recognized as one attractor.
# ----------------------------------------------------------------------

# Build the list of all 2^n states.
all_states = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]

attractors = {}  # canonical-cycle (tuple of states) -> cycle list

for start in all_states:
    seen_order = []          # states in the order first visited
    seen_index = {}          # state -> position in seen_order
    s = start
    while s not in seen_index:
        seen_index[s] = len(seen_order)
        seen_order.append(s)
        s = update(s)
    # s is the first repeated state: the cycle runs from it to the end.
    cycle = seen_order[seen_index[s]:]
    # Canonicalize by rotating so the smallest state leads.
    k = min(range(len(cycle)), key=lambda i: cycle[i])
    canon = tuple(cycle[k:] + cycle[:k])
    attractors[canon] = list(canon)

# ----------------------------------------------------------------------
# Report the attractors as arrow chains.
# ----------------------------------------------------------------------
print("Repressilator attractors (as arrow chains):")
fixed_points = []
cyclic = []
for canon, cycle in sorted(attractors.items()):
    chain = " -> ".join(fmt(st) for st in cycle) + " -> " + fmt(cycle[0])
    print("  " + chain)
    if len(cycle) == 1:
        fixed_points.append(cycle)
    else:
        cyclic.append(cycle)

# ----------------------------------------------------------------------
# Separate check: no fixed points, exactly two cyclic attractors,
# one of length 2 (000 <-> 111) and one of length 6.
# ----------------------------------------------------------------------
print()
print("Number of attractors total:", len(attractors))
print("Number of fixed points:", len(fixed_points))
print("Number of cyclic attractors:", len(cyclic))
cycle_lengths = sorted(len(c) for c in cyclic)
print("Cyclic attractor lengths:", cycle_lengths)

check = (len(fixed_points) == 0
         and len(cyclic) == 2
         and cycle_lengths == [2, 6])
print("Check passed (no fixed points, one 2-cycle + one 6-cycle):", check)

# One-sentence explanation of why the check confirms the result.
print()
print("Why the check confirms the result: the absence of any fixed point")
print("together with a single short 2-cycle (000<->111) and a single")
print("6-state cycle shows the ring of odd-length repression admits no")
print("steady state and instead sustains oscillation, the discrete analog")
print("of the continuous repressilator's limit cycle.")

# ----------------------------------------------------------------------
# Visualize the attractors as directed cycles.
# ----------------------------------------------------------------------
import math

fig, axes = plt.subplots(1, len(attractors), figsize=(5 * len(attractors), 5))
if len(attractors) == 1:
    axes = [axes]

for ax, (canon, cycle) in zip(axes, sorted(attractors.items())):
    m = len(cycle)
    # Place states evenly around a circle.
    angles = [math.pi / 2 - 2 * math.pi * i / m for i in range(m)]
    xs = [math.cos(a) for a in angles]
    ys = [math.sin(a) for a in angles]
    for i in range(m):
        j = (i + 1) % m
        ax.annotate("", xy=(xs[j], ys[j]), xytext=(xs[i], ys[i]),
                    arrowprops=dict(arrowstyle="->", color="steelblue", lw=2,
                                    shrinkA=15, shrinkB=15))
    for i in range(m):
        ax.text(xs[i], ys[i], fmt(cycle[i]), ha="center", va="center",
                fontsize=13, fontweight="bold",
                bbox=dict(boxstyle="circle", fc="lightyellow", ec="black"))
    ax.set_title("%d-state attractor" % m)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")

fig.suptitle("Boolean repressilator attractors (n = 3)", fontsize=14)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.6.1_s4.png")
