import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Boolean repressilator update rule on n = 3 genes ----
# State is a tuple (x, y, z) of bits.
# Ring of three repressions:
#   X(t+1) = NOT Z(t)
#   Y(t+1) = NOT X(t)
#   Z(t+1) = NOT Y(t)
def update(state):
    x, y, z = state
    return (1 - z, 1 - x, 1 - y)

# ---- Enumerate all 2^n = 8 states (attractor enumeration from 10C.4) ----
n = 3
all_states = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]

# ---- Find the attractor reachable from each start state ----
# Because the state space is finite and the dynamics deterministic, iterating
# from any start must eventually revisit a state; the set of states in that
# repeating loop is the attractor (a cycle; a fixed point is a cycle of length 1).
def attractor_from(start):
    seen = []          # ordered list of visited states
    s = start
    while s not in seen:
        seen.append(s)
        s = update(s)
    # s is the first repeated state -> loop starts there
    loop_start = seen.index(s)
    return seen[loop_start:]   # the cyclic attractor, in order

# ---- Collect the distinct attractors (canonicalize each cycle) ----
def canonical(cycle):
    # rotate so the lexicographically smallest state is first -> unique key
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(rotations)

attractors = {}   # canonical key -> cycle (as first found)
for start in all_states:
    cyc = attractor_from(start)
    key = canonical(cyc)
    if key not in attractors:
        attractors[key] = cyc

# ---- Format helpers ----
def bits(state):
    return "".join(str(b) for b in state)

def arrow_chain(cycle):
    # close the loop back to the first state to show it is cyclic
    return " -> ".join(bits(s) for s in cycle) + " -> " + bits(cycle[0])

# ---- Report attractors as arrow chains ----
print("Number of states enumerated:", len(all_states))
print("Number of distinct attractors:", len(attractors))
print()

fixed_points = []
cyclic = []
for i, cyc in enumerate(attractors.values(), start=1):
    print(f"Attractor {i} (length {len(cyc)}): {arrow_chain(cyc)}")
    if len(cyc) == 1:
        fixed_points.append(cyc)
    else:
        cyclic.append(cyc)

# ---- Separate check: no fixed points, exactly two cyclic attractors ----
print()
print("Number of fixed points found:", len(fixed_points))
print("Number of cyclic attractors found:", len(cyclic))
cycle_lengths = sorted(len(c) for c in cyclic)
print("Cyclic attractor lengths (sorted):", cycle_lengths)
print("Check: no fixed points, only cyclic attractors:", len(fixed_points) == 0)
print("Check: exactly two cyclic attractors of lengths [2, 6]:", cycle_lengths == [2, 6])
# Why this confirms the result: a fixed point would require a state equal to its
# own update, but every gene is the negation of another, so no state is stable;
# finding only the length-2 cycle 000<->111 and one length-6 cycle exhausts all
# 8 states and matches the discrete analog of the continuous limit cycle.
print()
print("Explanation: Because each gene's next value is the NOT of another gene,")
print("no state can equal its own successor, so there are no fixed points; the")
print("length-2 and length-6 cycles together cover all 8 states, confirming the")
print("six-state cycle is the discrete analog of the continuous limit cycle.")

# ---- Visualization: draw each attractor as a labeled ring/chain ----
import math
fig, axes = plt.subplots(1, len(attractors), figsize=(5 * len(attractors), 5))
if len(attractors) == 1:
    axes = [axes]
for ax, (i, cyc) in zip(axes, enumerate(attractors.values(), start=1)):
    m = len(cyc)
    angles = [2 * math.pi * k / m - math.pi / 2 for k in range(m)]
    xs = [math.cos(a) for a in angles]
    ys = [math.sin(a) for a in angles]
    for k in range(m):
        x0, y0 = xs[k], ys[k]
        x1, y1 = xs[(k + 1) % m], ys[(k + 1) % m]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="steelblue", lw=2))
    for k in range(m):
        ax.plot(xs[k], ys[k], "o", ms=28, color="white",
                markeredgecolor="black", zorder=3)
        ax.text(xs[k], ys[k], bits(cyc[k]), ha="center", va="center",
                fontsize=12, zorder=4)
    ax.set_title(f"Attractor {i} (length {m})")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")

fig.suptitle("Boolean repressilator attractors")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.6.1_s1.png")
