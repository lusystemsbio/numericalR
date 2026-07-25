import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Boolean network attractor enumeration (synchronous updating).
#
# A state is an n-tuple of 0/1 gene values.  The update rule
# f maps a state x(t) -> x(t+1) deterministically.  Because the
# state space is finite (2^n states), iterating from any start
# must eventually revisit a state; the loop of states between the
# first repeat and its earlier occurrence IS the attractor.
# ---------------------------------------------------------------


def boolean_attractors(update, n):
    """Enumerate all attractors of an n-gene Boolean network.

    update: function taking a state tuple -> next state tuple
    Returns a list of attractors, each a tuple of states forming
    the cycle (length 1 => fixed point, >1 => cyclic attractor).
    """
    attractors = []          # list of canonical attractor cycles
    seen_attractor = set()   # canonical keys of attractors already found

    # Run every one of the 2^n initial states to its first repeat.
    for start in range(2 ** n):
        # Decode integer 'start' into a state tuple of 0/1 bits.
        state = tuple((start >> i) & 1 for i in range(n))

        # Walk the trajectory, remembering the order states appear
        # and the step index at which each state was first seen.
        first_seen = {}   # state -> position in trajectory
        trajectory = []   # ordered list of visited states
        t = 0
        while state not in first_seen:
            first_seen[state] = t
            trajectory.append(state)
            state = tuple(update(state))  # synchronous update of all genes
            t += 1

        # 'state' is the first repeated state.  The attractor is the
        # slice of the trajectory from its first appearance onward.
        cycle_start = first_seen[state]
        cycle = tuple(trajectory[cycle_start:])

        # Canonicalize the cycle (rotate so its smallest state leads)
        # so the same attractor found from different starts matches.
        rot = min(range(len(cycle)), key=lambda k: cycle[k:] + cycle[:k])
        canon = cycle[rot:] + cycle[:rot]

        if canon not in seen_attractor:
            seen_attractor.add(canon)
            attractors.append(canon)

    return attractors


# ---------------------------------------------------------------
# Update rule applied in 10C.5 and 10C.6 (a small 3-gene circuit).
#   g0(t+1) = g2
#   g1(t+1) = g0 AND g2
#   g2(t+1) = NOT g1
# ---------------------------------------------------------------
def update(x):
    g0, g1, g2 = x
    return (g2, g0 and g2, 1 - g1)


n = 3
attractors = boolean_attractors(update, n)

print("Number of genes n =", n)
print("Total initial states explored = 2^n =", 2 ** n)
print("Number of attractors found =", len(attractors))

fixed_points = 0
cyclic = 0
for idx, att in enumerate(attractors):
    kind = "fixed-point" if len(att) == 1 else "cyclic"
    if len(att) == 1:
        fixed_points += 1
    else:
        cyclic += 1
    print("Attractor {} ({}, period {}): {}".format(idx + 1, kind, len(att), att))

print("Number of fixed-point attractors =", fixed_points)
print("Number of cyclic attractors =", cyclic)

# ---------------------------------------------------------------
# Separate check: every one of the 2^n states must flow into
# exactly one of the enumerated attractors.  We rebuild each
# state's basin by iterating to a repeat and confirm its cycle
# matches (in canonical form) one of the attractors above.
# This confirms the routine returns EVERY attractor, because if
# any attractor were missing, some state's cycle would fail to
# match the enumerated set.
# ---------------------------------------------------------------
canon_set = set(attractors)
covered = 0
all_ok = True
for start in range(2 ** n):
    state = tuple((start >> i) & 1 for i in range(n))
    first_seen = {}
    trajectory = []
    while state not in first_seen:
        first_seen[state] = len(trajectory)
        trajectory.append(state)
        state = tuple(update(state))
    cycle = tuple(trajectory[first_seen[state]:])
    rot = min(range(len(cycle)), key=lambda k: cycle[k:] + cycle[:k])
    canon = cycle[rot:] + cycle[:rot]
    if canon in canon_set:
        covered += 1
    else:
        all_ok = False

print("States that reach an enumerated attractor =", covered, "of", 2 ** n)
print("Check passed (all states covered) =", all_ok)

# ---------------------------------------------------------------
# Visualize: bar chart of attractor periods.
# ---------------------------------------------------------------
periods = [len(a) for a in attractors]
labels = ["A{}".format(i + 1) for i in range(len(attractors))]
plt.figure()
plt.bar(labels, periods, color="steelblue")
plt.xlabel("Attractor")
plt.ylabel("Period (number of states in cycle)")
plt.title("Boolean network attractors (synchronous update, n={})".format(n))
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.4.1_s5.png")
