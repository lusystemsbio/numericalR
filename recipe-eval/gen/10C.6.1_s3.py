import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import product

# --- Boolean repressilator update rule on n = 3 genes -------------------
# State = (x, y, z), each 0/1.  Ring of three repressions:
#   X(t+1) = NOT Z(t)
#   Y(t+1) = NOT X(t)
#   Z(t+1) = NOT Y(t)
def update(state):
    x, y, z = state
    return (1 - z, 1 - x, 1 - y)

# --- Boolean attractor enumeration (method from 10C.4) ------------------
# Every state is eventually periodic in a finite deterministic system.
# From each of the 2^n states, follow the trajectory until a state repeats;
# the repeating portion is the attractor (cycle).  We canonicalize each
# cycle (rotate so its lexicographically smallest state is first) so that
# each distinct attractor is recorded exactly once.
n = 3
all_states = list(product([0, 1], repeat=n))   # all eight states

def trajectory_to_cycle(start):
    seen = {}          # state -> position in path
    path = []
    s = start
    # walk forward until we revisit a state
    while s not in seen:
        seen[s] = len(path)
        path.append(s)
        s = update(s)
    # s is the first repeated state; the cycle is path[seen[s]:]
    return path[seen[s]:]

def canonical(cycle):
    # rotate the cycle so it begins at its smallest state -> unique key
    i = min(range(len(cycle)), key=lambda k: cycle[k])
    return tuple(cycle[i:] + cycle[:i])

attractors = {}
for st in all_states:
    cyc = canonical(trajectory_to_cycle(st))
    attractors[cyc] = len(cyc)   # store cycle with its length (period)

# --- Format attractors as arrow chains ----------------------------------
def bits(s):
    return "".join(str(b) for b in s)

def arrow_chain(cycle):
    # close the loop by returning to the first state
    chain = list(cycle) + [cycle[0]]
    return " -> ".join(bits(s) for s in chain)

print("Repressilator attractors (arrow chains):")
attractor_list = sorted(attractors.keys(), key=lambda c: (len(c), c))
for cyc in attractor_list:
    print("period %d: %s" % (len(cyc), arrow_chain(cyc)))

# --- Separate check: no fixed points, exactly two cyclic attractors -----
fixed_points = [c for c in attractor_list if len(c) == 1]
cyclic = [c for c in attractor_list if len(c) > 1]

print()
print("Number of fixed points (period-1 attractors):", len(fixed_points))
print("Number of cyclic attractors (period > 1):", len(cyclic))
print("Total distinct attractors:", len(attractor_list))
for cyc in cyclic:
    print("Cyclic attractor period:", len(cyc))

# verify the two expected cycles are present
two_cycle = canonical(trajectory_to_cycle((0, 0, 0)))   # 000 <-> 111
six_cycle_seed = canonical(trajectory_to_cycle((1, 0, 0)))
print("Found the 000<->111 two-cycle:", two_cycle in attractors and len(two_cycle) == 2)
print("Found a six-state cycle:", six_cycle_seed in attractors and len(six_cycle_seed) == 6)

# Why the check confirms the result: with no fixed points every one of the
# eight states must flow into one of the two cycles, and 2 + 6 = 8 accounts
# for all states exactly once, so these two cycles are the complete and
# only attractors of the system.
print()
print("Explanation: with zero fixed points the two cycles of length 2 and 6"
      " partition all 2+6=8 states, so they are exactly the system's attractors.")

# --- Draw the attractors as ring diagrams -------------------------------
import math
fig, axes = plt.subplots(1, len(attractor_list), figsize=(5 * len(attractor_list), 5))
if len(attractor_list) == 1:
    axes = [axes]
for ax, cyc in zip(axes, attractor_list):
    m = len(cyc)
    angles = [2 * math.pi * k / m + math.pi / 2 for k in range(m)]
    xs = [math.cos(a) for a in angles]
    ys = [math.sin(a) for a in angles]
    for k in range(m):
        x0, y0 = xs[k], ys[k]
        x1, y1 = xs[(k + 1) % m], ys[(k + 1) % m]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="steelblue", lw=2))
        ax.text(x0 * 1.15, y0 * 1.15, bits(cyc[k]),
                ha="center", va="center", fontsize=13, fontweight="bold")
    ax.set_title("period %d cycle" % m)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
fig.suptitle("Boolean repressilator attractors", fontsize=15)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.6.1_s3.png")
