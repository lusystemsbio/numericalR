import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import product

# ---------------------------------------------------------------
# Represent a network state as a tuple of 0/1 of length n.
# An "update" is a function: state (tuple) -> next state (tuple).
# Synchronous updating: all genes update at once from the SAME x(t).
# ---------------------------------------------------------------


def boolean_attractors(update, n):
    """Generic driver: run every one of the 2^n start states to its
    first repeat, and collect the distinct attractors (fixed points
    or cycles) that the dynamics fall into."""

    attractors = []          # list of attractors, each a list of states (the cycle)
    seen_states = set()      # states already known to lead to a known attractor

    # Enumerate ALL 2^n initial conditions exhaustively.
    for start in product((0, 1), repeat=n):

        if start in seen_states:
            # This start already flows into an attractor we found; skip.
            continue

        # Iterate x(t+1) = f(x(t)) and record the trajectory order,
        # stopping as soon as we revisit a state (the "first repeat").
        trajectory = []              # ordered list of visited states
        position = {}                # state -> index in trajectory
        x = start
        while x not in position:
            position[x] = len(trajectory)
            trajectory.append(x)
            x = tuple(update(x))     # one synchronous step

        # The repeated state 'x' marks where the trajectory closes.
        # Everything from that index onward is the attractor (cycle);
        # the earlier states are the transient tail leading into it.
        loop_start = position[x]
        cycle = trajectory[loop_start:]
        canonical = min(_rotations(cycle))   # canonical form for de-dup

        # Mark all trajectory states as accounted for (they all flow here).
        seen_states.update(trajectory)

        # Register the attractor only once.
        if canonical not in {min(_rotations(a)) for a in attractors}:
            attractors.append(cycle)

    return attractors


def _rotations(cycle):
    """All rotations of a cycle, so cyclic attractors compare equal
    regardless of which state we happened to enter them on."""
    return [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]


# ---------------------------------------------------------------
# The specific update rule applied in 10C.5 and 10C.6.
# A 3-gene circuit:
#   gene0 next = gene1 AND gene2
#   gene1 next = NOT gene0
#   gene2 next = gene1 OR gene2
# ---------------------------------------------------------------


def update_rule(x):
    g0, g1, g2 = x
    return (
        g1 & g2,          # f_0
        1 - g0,           # f_1 (NOT)
        g1 | g2,          # f_2
    )


n = 3
attractors = boolean_attractors(update_rule, n)

# ---- Report the enumerated attractors ----
print("Number of genes n = %d" % n)
print("Number of start states explored = %d" % (2 ** n))
print("Number of attractors found = %d" % len(attractors))
for k, cyc in enumerate(attractors):
    kind = "fixed-point" if len(cyc) == 1 else ("cycle (period %d)" % len(cyc))
    states = " -> ".join("".join(map(str, s)) for s in cyc)
    print("Attractor %d [%s]: %s" % (k + 1, kind, states))

# ---------------------------------------------------------------
# Separate check: confirm the routine returns EVERY attractor.
# Independently, from each of the 2^n starts, follow the dynamics to
# its first repeat and record which attractor (canonical cycle) it
# reaches. Then verify: (a) every start maps to one of the returned
# attractors, and (b) every returned attractor is actually reached.
# ---------------------------------------------------------------


def _reach_attractor(start):
    position = {}
    traj = []
    x = start
    while x not in position:
        position[x] = len(traj)
        traj.append(x)
        x = tuple(update_rule(x))
    return min(_rotations(traj[position[x]:]))


returned_canon = {min(_rotations(a)) for a in attractors}
reached_canon = set()
all_starts_covered = True
for start in product((0, 1), repeat=n):
    c = _reach_attractor(start)
    reached_canon.add(c)
    if c not in returned_canon:
        all_starts_covered = False

check_all_reached = reached_canon == returned_canon

print("Check - every start state lands in a returned attractor:", all_starts_covered)
print("Check - returned attractors exactly equal reached attractors:", check_all_reached)
print("Check PASSED:", all_starts_covered and check_all_reached)
# One sentence: because we independently trace all 2^n states to their
# first repeat, this check confirms the result by showing the routine's
# attractor set is exactly the set of terminal cycles reachable from
# every possible initial condition -- none missing and none spurious.

# ---------------------------------------------------------------
# Figure: state-transition graph coloring attractor states.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
import math

states = list(product((0, 1), repeat=n))
pos = {s: (math.cos(2 * math.pi * i / len(states)),
           math.sin(2 * math.pi * i / len(states)))
       for i, s in enumerate(states)}
attractor_states = set()
for a in attractors:
    attractor_states.update(a)

for s in states:
    nxt = tuple(update_rule(s))
    x0, y0 = pos[s]
    x1, y1 = pos[nxt]
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color="0.6", lw=1))

for s in states:
    x0, y0 = pos[s]
    color = "crimson" if s in attractor_states else "lightsteelblue"
    ax.plot(x0, y0, "o", ms=26, color=color, zorder=3)
    ax.text(x0, y0, "".join(map(str, s)), ha="center", va="center",
            zorder=4, fontsize=9)

ax.set_title("Boolean network state-transition graph\n(red = attractor states)")
ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.4.1_s1.png")
print("Figure saved.")
