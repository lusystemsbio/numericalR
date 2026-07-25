import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Boolean network attractor enumeration (synchronous updating).
#
# State: an n-bit tuple x = (x_0, ..., x_{n-1}), each x_i in {0,1}.
# Dynamics: x(t+1) = update(x(t)), applied to ALL genes at once.
# Because the state space is finite (2^n states), iterating from any
# start MUST eventually repeat a state; the repeated states form the
# attractor (a fixed point if length 1, a cycle if length > 1).
# -------------------------------------------------------------------

def boolean_attractors(update, n):
    """Enumerate all attractors of an n-gene Boolean network exhaustively.

    'update' maps an n-tuple of 0/1 to the next n-tuple (synchronous).
    We run every one of the 2^n start states forward to its first repeat.
    Returns a list of attractors, each a list of states (a canonical cycle).
    """
    attractors = []            # collected distinct attractors
    seen_attractor_states = set()  # states already known to lie on some attractor

    for start in range(2 ** n):
        # decode integer 'start' into an n-bit state tuple (bit i -> gene i)
        state = tuple((start >> i) & 1 for i in range(n))

        # follow the trajectory, recording the order states are first visited
        path = []              # ordered list of states on this trajectory
        index_of = {}          # state -> position in path (for detecting repeat)
        while state not in index_of:
            index_of[state] = len(path)
            path.append(state)
            state = update(state)   # synchronous one-step update

        # 'state' is the first repeated state; the attractor is the tail
        # of the path from that state's first occurrence to the end (the loop)
        cycle = path[index_of[state]:]

        # canonicalize the cycle (rotate so its smallest state is first) so
        # the same attractor found from different starts is recognized as one
        m = cycle.index(min(cycle))
        canonical = tuple(cycle[m:] + cycle[:m])

        if canonical[0] not in seen_attractor_states or \
           all(s not in seen_attractor_states for s in canonical):
            # only add if we have not already recorded this attractor
            if not any(set(canonical) == set(a) for a in attractors):
                attractors.append(list(canonical))
                seen_attractor_states.update(canonical)

    return attractors


# -------------------------------------------------------------------
# The update rule applied in 10C.5 and 10C.6: a 3-gene circuit.
#   gene 0 (A): A' = B AND C
#   gene 1 (B): B' = NOT A
#   gene 2 (C): C' = B OR C
# -------------------------------------------------------------------
def update(x):
    A, B, C = x
    A_next = 1 if (B and C) else 0
    B_next = 1 if (not A) else 0
    C_next = 1 if (B or C) else 0
    return (A_next, B_next, C_next)


n = 3

# ---- Main result: enumerate the attractors ----
attractors = boolean_attractors(update, n)

print("Number of genes n:", n)
print("Total start states scanned (2^n):", 2 ** n)
print("Number of attractors found:", len(attractors))

fixed_points = []
cyclic = []
for k, att in enumerate(attractors):
    kind = "fixed-point" if len(att) == 1 else ("cycle (period %d)" % len(att))
    print("Attractor %d: %s : %s" % (k, kind, att))
    if len(att) == 1:
        fixed_points.append(att[0])
    else:
        cyclic.append(att)

print("Number of fixed-point attractors:", len(fixed_points))
print("Number of cyclic attractors:", len(cyclic))

# -------------------------------------------------------------------
# Separate check: confirm the routine returns EVERY attractor by
# verifying that every one of the 2^n states flows into exactly one
# of the enumerated attractors (i.e. the attractors partition the
# reachable long-run behavior of the whole state space).
# -------------------------------------------------------------------
attractor_state_set = set()
for att in attractors:
    attractor_state_set.update(att)

all_states_covered = True
for start in range(2 ** n):
    state = tuple((start >> i) & 1 for i in range(n))
    # iterate enough steps to certainly reach an attractor (<= 2^n steps)
    for _ in range(2 ** n):
        state = update(state)
    landed_in_attractor = state in attractor_state_set
    if not landed_in_attractor:
        all_states_covered = False

print("Every start state lands in an enumerated attractor:", all_states_covered)
print("Total states lying on some attractor:", len(attractor_state_set))

# One-sentence explanation of why this check confirms the result:
# Because the state space is finite every trajectory must end on a cycle,
# so if every one of the 2^n starts flows into one of the collected sets,
# those sets are exactly all the fixed-point and cyclic attractors.
print("Check rationale: since every finite-state trajectory must terminate "
      "on a cycle, all 2^n starts landing in the collected sets proves those "
      "sets are the complete list of fixed-point and cyclic attractors.")

# -------------------------------------------------------------------
# Visualization: the state-transition graph, attractor states shaded.
# -------------------------------------------------------------------
import math

fig, ax = plt.subplots(figsize=(7, 7))
positions = {}
for s in range(2 ** n):
    ang = 2 * math.pi * s / (2 ** n)
    positions[s] = (math.cos(ang), math.sin(ang))

def to_int(state):
    return sum(b << i for i, b in enumerate(state))

# draw transition edges
for s in range(2 ** n):
    state = tuple((s >> i) & 1 for i in range(n))
    t = to_int(update(state))
    x0, y0 = positions[s]
    x1, y1 = positions[t]
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color="gray", alpha=0.6))

# draw nodes
for s in range(2 ** n):
    state = tuple((s >> i) & 1 for i in range(n))
    x0, y0 = positions[s]
    on_attr = state in attractor_state_set
    ax.scatter([x0], [y0], s=600,
               color=("tomato" if on_attr else "lightsteelblue"),
               edgecolor="black", zorder=3)
    ax.text(x0, y0, "".join(map(str, state)), ha="center", va="center",
            fontsize=9, zorder=4)

ax.set_title("Boolean network state-transition graph\n(red = attractor states)")
ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.4.1_s2.png")
print("Figure saved to 10C.4.1_s2.png")
