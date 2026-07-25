import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Boolean network dynamics with SYNCHRONOUS updating.
# A state is a tuple of 0/1 values, one per gene: x = (x_0,...,x_{n-1}).
# Every gene updates at once via x_i(t+1) = f_i(x(t)).
# From any start, iterating the deterministic map eventually repeats,
# and the repeated part is the attractor (fixed point or cycle).
# ---------------------------------------------------------------


def boolean_attractors(update, n):
    """Generic driver: run every one of the 2^n start states to its first
    repeat, and collect the distinct attractors reached.

    'update' maps a state tuple -> next state tuple (the vector f).
    Returns a list of attractors; each attractor is a tuple of the states
    on the cycle (length 1 == fixed point)."""

    attractors = []          # list of attractor cycles (as tuples of states)
    attractor_of = {}        # state -> index into 'attractors' (its basin)

    # Enumerate ALL 2^n initial states exhaustively.
    for start_int in range(2 ** n):
        # Decode the integer into a Boolean state tuple (bit i -> gene i).
        state = tuple((start_int >> i) & 1 for i in range(n))

        # Walk the trajectory forward, remembering the order we saw states,
        # until we hit a state we have already visited on THIS walk (a repeat)
        # or a state whose attractor we already know from a previous walk.
        seen_order = []          # states visited on this trajectory, in order
        seen_index = {}          # state -> position in seen_order
        while state not in seen_index and state not in attractor_of:
            seen_index[state] = len(seen_order)
            seen_order.append(state)
            state = update(state)  # synchronous step: apply f to the whole vector

        if state in attractor_of:
            # Merged into a previously discovered attractor; label the tail too.
            idx = attractor_of[state]
        else:
            # First repeat is within this walk: the cycle starts at that repeat.
            cycle_start = seen_index[state]
            cycle = tuple(seen_order[cycle_start:])   # the attractor states
            # Canonicalize the cycle (rotate so smallest state is first) so the
            # same attractor found from different starts is recognized as one.
            k = min(range(len(cycle)), key=lambda j: cycle[j])
            cycle = cycle[k:] + cycle[:k]
            idx = len(attractors)
            attractors.append(cycle)
            for s in cycle:
                attractor_of[s] = idx

        # Every transient state on this walk falls into the same attractor.
        for s in seen_order:
            attractor_of[s] = idx

    return attractors, attractor_of


# ---------------------------------------------------------------
# The specific update rule used in 10C.5 and 10C.6 (a 3-gene circuit).
#   gene0(t+1) = NOT gene2
#   gene1(t+1) = gene0 AND gene2
#   gene2(t+1) = gene0 OR  gene1
# ---------------------------------------------------------------
def update(x):
    x0, x1, x2 = x
    n0 = 1 - x2
    n1 = x0 & x2
    n2 = x0 | x1
    return (n0, n1, n2)


n = 3
attractors, attractor_of = boolean_attractors(update, n)

# ---- Report the attractors enumerated for this update rule ----
fixed_points = [a for a in attractors if len(a) == 1]
cycles = [a for a in attractors if len(a) > 1]

print("Number of genes n =", n)
print("Total states (2^n) =", 2 ** n)
print("Number of attractors found =", len(attractors))
print("Number of fixed-point attractors =", len(fixed_points))
print("Number of cyclic attractors =", len(cycles))
for i, a in enumerate(attractors):
    kind = "fixed point" if len(a) == 1 else ("cycle length %d" % len(a))
    print("Attractor %d (%s): %s" % (i, kind, list(a)))

# ---------------------------------------------------------------
# SEPARATE CHECK: confirm the routine returns EVERY fixed point and
# EVERY cyclic attractor. Two independent verifications:
#  (1) coverage: every one of the 2^n states is assigned to some attractor;
#  (2) closure : applying the update to each reported attractor reproduces
#                that same attractor set (so it is genuinely invariant).
# ---------------------------------------------------------------
all_states = set(tuple((s >> i) & 1 for i in range(n)) for s in range(2 ** n))
coverage_ok = set(attractor_of.keys()) == all_states

closure_ok = True
for a in attractors:
    a_set = set(a)
    img = set(update(s) for s in a)      # one synchronous step of the whole cycle
    if img != a_set:
        closure_ok = False

# Independently brute-force every fixed point (f(x)==x) and confirm each is found.
brute_fixed = set(s for s in all_states if update(s) == s)
found_fixed = set(a[0] for a in fixed_points)
fixed_match = brute_fixed == found_fixed

print("Check coverage (all 2^n states land in an attractor):", coverage_ok)
print("Check closure  (each attractor maps onto itself):", closure_ok)
print("Check fixed points match brute-force f(x)==x scan:", fixed_match)
print("All checks passed:", coverage_ok and closure_ok and fixed_match)

# One-sentence explanation of why the check confirms the result:
print("Explanation: because we exhausted all 2^n starts, coverage plus "
      "self-closure guarantees the reported list contains every attractor "
      "and nothing that is not actually an attractor.")

# ---------------------------------------------------------------
# Visualization: basin sizes per attractor.
# ---------------------------------------------------------------
basin_sizes = [0] * len(attractors)
for s, idx in attractor_of.items():
    basin_sizes[idx] += 1

labels = []
for i, a in enumerate(attractors):
    labels.append(("FP " if len(a) == 1 else "Cyc%d " % len(a)) + "#%d" % i)

plt.figure(figsize=(7, 4))
plt.bar(range(len(attractors)), basin_sizes, color="steelblue")
plt.xticks(range(len(attractors)), labels, rotation=30, ha="right")
plt.ylabel("Basin size (number of start states)")
plt.title("Boolean network attractors (synchronous, n=%d)" % n)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/10C.4.1_s4.png")

for i, sz in enumerate(basin_sizes):
    print("Basin size of attractor %d = %d" % (i, sz))
print("Sum of basin sizes =", sum(basin_sizes))
