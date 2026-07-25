import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Boolean network dynamics: synchronous updating.
# A state is a tuple of 0/1 values, one per gene.  The update rule maps
# the whole state at time t to the whole state at time t+1 via
#   x_i(t+1) = f_i(x(t)).
# Iterating from any start eventually repeats (finite state space), and
# the repeated portion is the attractor (fixed point if length 1,
# cyclic if length > 1).
# ----------------------------------------------------------------------


def int_to_state(k, n):
    # decode integer k in [0, 2^n) into an n-bit state tuple (gene 0 = MSB)
    return tuple((k >> (n - 1 - i)) & 1 for i in range(n))


def boolean_attractors(update, n):
    # Generic driver: run every one of the 2^n start states to its first
    # repeat, collect the distinct attractors (as canonical cycles).
    attractors = []           # list of canonical cycles (tuples of states)
    seen_attractor_sets = []  # frozensets of states, to dedupe attractors

    for k in range(2 ** n):               # loop over ALL initial states
        state = int_to_state(k, n)
        trajectory = []                   # ordered states visited this run
        index_of = {}                     # state -> position in trajectory

        # iterate synchronously until we hit a state we've already seen
        while state not in index_of:
            index_of[state] = len(trajectory)
            trajectory.append(state)
            state = update(state)         # apply x(t+1) = f(x(t))

        # 'state' is the first repeat; the cycle is everything from its
        # first occurrence onward -- that closed loop is the attractor
        start = index_of[state]
        cycle = trajectory[start:]

        # canonicalize the cycle so the same attractor found from
        # different starts / different rotations is recorded once
        cyc_set = frozenset(cycle)
        if cyc_set not in seen_attractor_sets:
            seen_attractor_sets.append(cyc_set)
            # rotate so the lexicographically smallest state leads
            m = cycle.index(min(cycle))
            canonical = tuple(cycle[m:] + cycle[:m])
            attractors.append(canonical)

    return attractors


# ----------------------------------------------------------------------
# The specific update rule applied in 10C.5 and 10C.6.
# A small 3-gene circuit:
#   g0(t+1) = g1 AND (NOT g2)
#   g1(t+1) = g0 OR  g2
#   g2(t+1) = NOT g1
# ----------------------------------------------------------------------
def update_rule(x):
    g0, g1, g2 = x
    n0 = g1 & (1 - g2)
    n1 = g0 | g2
    n2 = 1 - g1
    return (n0, n1, n2)


n = 3
attractors = boolean_attractors(update_rule, n)


def fmt(state):
    return "".join(str(b) for b in state)


print("Number of genes n:", n)
print("Number of start states enumerated:", 2 ** n)
print("Number of distinct attractors found:", len(attractors))

fixed_points = [a for a in attractors if len(a) == 1]
cyclic = [a for a in attractors if len(a) > 1]
print("Number of fixed-point attractors:", len(fixed_points))
print("Number of cyclic attractors:", len(cyclic))

for i, a in enumerate(attractors):
    kind = "fixed-point" if len(a) == 1 else "cyclic"
    print("Attractor {} ({}, period {}): {}".format(
        i + 1, kind, len(a), " -> ".join(fmt(s) for s in a)))

# ----------------------------------------------------------------------
# Separate check: confirm the routine returns EVERY fixed point and every
# cyclic attractor.  Independently, (1) every state must flow into exactly
# one recorded attractor, and (2) a brute-force fixed-point scan
# (states with update(x) == x) must match the recorded fixed points.
# ----------------------------------------------------------------------

# (1) basin coverage: which recorded attractor each start state ends in
attractor_of_state = {}
for idx, a in enumerate(attractors):
    for s in a:
        attractor_of_state[s] = idx

covered = 0
for k in range(2 ** n):
    state = int_to_state(k, n)
    seen = {}
    while state not in seen:
        seen[state] = True
        state = update_rule(state)
    # 'state' is a state on the attractor this start falls into
    if state in attractor_of_state:
        covered += 1

print("States that provably reach a recorded attractor:",
      covered, "of", 2 ** n)

# (2) brute-force fixed points computed independently of the driver
brute_fixed = [int_to_state(k, n) for k in range(2 ** n)
               if update_rule(int_to_state(k, n)) == int_to_state(k, n)]
recorded_fixed = sorted(a[0] for a in fixed_points)
print("Brute-force fixed points:", [fmt(s) for s in sorted(brute_fixed)])
print("Recorded  fixed points:", [fmt(s) for s in recorded_fixed])
print("Fixed-point sets match:", sorted(brute_fixed) == recorded_fixed)
print("All states covered:", covered == 2 ** n)

# This check confirms the result because a routine that both maps every one
# of the 2^n start states into a recorded attractor AND reproduces the
# independent brute-force fixed-point list cannot have missed any attractor.

# ----------------------------------------------------------------------
# Visualization: state-transition graph laid out so attractors stand out.
# ----------------------------------------------------------------------
import math

fig, ax = plt.subplots(figsize=(7, 7))
states = [int_to_state(k, n) for k in range(2 ** n)]
pos = {}
for j, s in enumerate(states):
    ang = 2 * math.pi * j / len(states)
    pos[s] = (math.cos(ang), math.sin(ang))

for s in states:
    ns = update_rule(s)
    x0, y0 = pos[s]
    x1, y1 = pos[ns]
    on_attr = s in attractor_of_state
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->",
                                color="crimson" if on_attr else "gray",
                                lw=2 if on_attr else 1, alpha=0.9))

for s in states:
    x0, y0 = pos[s]
    ax.plot(x0, y0, "o", ms=22,
            color="gold" if s in attractor_of_state else "lightsteelblue",
            zorder=3)
    ax.text(x0, y0, fmt(s), ha="center", va="center",
            fontsize=9, zorder=4)

ax.set_title("Synchronous Boolean network: state-transition graph\n"
             "(red edges / gold nodes = attractors)")
ax.set_aspect("equal")
ax.axis("off")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.4.1_s3.png")
print("Figure saved.")
