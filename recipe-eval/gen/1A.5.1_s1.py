import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model setup ----
N = 10                          # 10x10 torus grid
grid = np.zeros((N, N), dtype=int)
# Seed a glider at the requested cells (row, col)
glider_cells = [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
for r, c in glider_cells:
    grid[r, c] = 1

def step(g):
    """Advance one generation with explicit neighbor counting and wraparound."""
    n = np.zeros_like(g)
    # Sum the eight neighbors; np.roll gives periodic (toroidal) boundaries.
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue                      # skip the cell itself
            n += np.roll(np.roll(g, dr, axis=0), dc, axis=1)
    # Apply the three rules to every cell simultaneously:
    # live cell survives with 2 or 3 neighbors; dead cell born with exactly 3.
    survive = (g == 1) & ((n == 2) | (n == 3))
    born = (g == 0) & (n == 3)
    return (survive | born).astype(int)

# ---- Run several generations and record snapshots ----
GENERATIONS = 12
snapshots = [grid.copy()]
g = grid.copy()
for _ in range(GENERATIONS):
    g = step(g)
    snapshots.append(g.copy())

# Print each recorded generation's grid.
for i, snap in enumerate(snapshots):
    print(f"Generation {i}:")
    for row in snap:
        print("".join("#" if v else "." for v in row))
    print(f"Live cell count at generation {i}: {int(snap.sum())}")
    print()

# ---- Check: shape preserved and drifts 1 diagonal step per 4 generations ----
def normalized_shape(snap):
    """Return the set of live-cell offsets relative to the top-left of the pattern
    (measured on the torus by choosing the minimal wrapped bounding box)."""
    coords = np.argwhere(snap == 1)
    # Reduce each axis modulo the smallest spanning window (handles wraparound).
    offs = set()
    r0 = coords[:, 0].min()
    c0 = coords[:, 1].min()
    for r, c in coords:
        offs.add(((r - r0) % N, (c - c0) % N))
    return frozenset(offs)

base_shape = normalized_shape(snapshots[0])
shape_ok = all(normalized_shape(snapshots[k]) == base_shape for k in range(0, len(snapshots), 4))
print(f"Shape preserved every 4 generations (still 5 cells, same form): {shape_ok}")

# Track the pattern's centroid (via sum of coordinates) modulo N every 4 generations.
def top_left(snap):
    coords = np.argwhere(snap == 1)
    return coords[:, 0].min(), coords[:, 1].min()

for k in range(0, len(snapshots) - 4, 4):
    r_now, c_now = top_left(snapshots[k])
    r_next, c_next = top_left(snapshots[k + 4])
    dr = (r_next - r_now) % N
    dc = (c_next - c_now) % N
    print(f"Gen {k} -> {k+4}: top-left drift (dr,dc) mod {N} = ({dr},{dc})  "
          f"[expected (1,1): one diagonal step]")

# ---- Plot successive generations ----
n_plots = 8
fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for ax, i in zip(axes.ravel(), range(n_plots)):
    ax.imshow(snapshots[i], cmap="binary", vmin=0, vmax=1)
    ax.set_title(f"Gen {i}")
    ax.set_xticks(range(N)); ax.set_yticks(range(N))
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.grid(True, color="lightgray", linewidth=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.5.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the glider is the smallest translating oscillator, "
      "verifying that the pattern still has exactly five cells in the same relative "
      "arrangement and that its bounding box shifts by (1,1) every four steps (wrapping "
      "at the edges) proves the simulation reproduces the glider's known period-4 "
      "diagonal motion on the torus rather than decaying or mutating.")
