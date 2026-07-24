import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Grid setup: 10x10, seeded with a glider ----
N = 10
grid = np.zeros((N, N), dtype=int)
# Glider cells (row, col)
glider_cells = [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
for (r, c) in glider_cells:
    grid[r, c] = 1


def count_neighbors(g):
    # Count live neighbors of every cell with wraparound (toroidal boundary).
    # np.roll shifts the whole grid; summing the 8 shifted copies gives,
    # for each cell, how many of its 8 neighbors are alive.
    neighbors = np.zeros_like(g)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue  # skip the cell itself
            neighbors += np.roll(np.roll(g, dr, axis=0), dc, axis=1)
    return neighbors


def step(g):
    # Apply the three rules to all cells at once.
    n = count_neighbors(g)
    new = np.zeros_like(g)
    # Live cell with 2 or 3 live neighbors survives.
    new[(g == 1) & ((n == 2) | (n == 3))] = 1
    # Dead cell with exactly 3 live neighbors becomes alive.
    new[(g == 0) & (n == 3)] = 1
    # Every other cell is dead next generation (already 0 by default).
    return new


def normalized_shape(g):
    # Return the set of live-cell coordinates translated so the minimum
    # row and column are 0, giving a position-independent shape signature.
    coords = np.argwhere(g == 1)
    coords = coords - coords.min(axis=0)
    return {tuple(c) for c in coords}


def centroid(g):
    coords = np.argwhere(g == 1)
    return coords.mean(axis=0)


# ---- Run the simulation and collect generations ----
generations = []
grids = []
g = grid.copy()
n_gens = 16
for gen in range(n_gens + 1):
    generations.append(gen)
    grids.append(g.copy())
    print(f"Generation {gen}: live_cells={int(g.sum())}")
    g = step(g)

# ---- Print several successive generations as text ----
for gen in range(0, 9):
    print(f"\n--- Generation {gen} grid ---")
    for row in grids[gen]:
        print("".join("#" if v else "." for v in row))

# ---- Check 1: glider keeps its five-cell shape ----
base_shape = normalized_shape(grids[0])
all_five = all(int(gg.sum()) == 5 for gg in grids)
print(f"\nGlider always has 5 live cells: {all_five}")

# Shape should repeat (up to translation) every 4 generations, since the
# glider only rotates through 4 phases before returning to its start shape.
shape_matches_every_4 = all(
    normalized_shape(grids[k]) == base_shape for k in range(0, n_gens + 1, 4)
)
print(f"Shape identical (up to translation) every 4 generations: {shape_matches_every_4}")

# ---- Check 2: drifts one step diagonally every four generations ----
# Compare centroid at gen 0 and gen 4, accounting for torus wraparound.
c0 = centroid(grids[0])
c4 = centroid(grids[4])
raw_shift = c4 - c0
# Wrap each component into (-N/2, N/2] to handle edge wraparound.
wrapped_shift = ((raw_shift + N / 2) % N) - N / 2
print(f"\nCentroid gen 0: ({c0[0]:.2f}, {c0[1]:.2f})")
print(f"Centroid gen 4: ({c4[0]:.2f}, {c4[1]:.2f})")
print(f"Diagonal drift over 4 generations (row, col): "
      f"({wrapped_shift[0]:.2f}, {wrapped_shift[1]:.2f})")

# ---- Check that it reappears on the opposite side after wrapping ----
# Track the max row/col over the run to show the glider crosses an edge.
reached_edge = any((gg[N - 1, :].sum() > 0) or (gg[:, N - 1].sum() > 0) for gg in grids)
print(f"Glider reaches a grid edge during the run: {reached_edge}")
wrapped_back = any((gg[0, :].sum() > 0) for gg in grids[5:])
print(f"Glider reappears on opposite (top) side after wrapping: {wrapped_back}")

# ---- Plot several successive generations ----
show_gens = [0, 1, 2, 3, 4, 8, 12, 16]
fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for ax, gen in zip(axes.flat, show_gens):
    ax.imshow(grids[gen], cmap="binary", vmin=0, vmax=1)
    ax.set_title(f"Generation {gen}")
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.grid(True, color="lightgray", linewidth=0.5)
    ax.tick_params(length=0, labelsize=6)
fig.suptitle("Conway's Game of Life: glider on a 10x10 torus")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.5.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: A glider is confirmed if its live-cell count stays 5, "
      "its shape (up to translation) recurs every 4 generations, and its "
      "centroid moves exactly one cell diagonally per 4-generation period "
      "with wraparound, since that is precisely the defining periodic "
      "translating behavior of a glider on a torus.")
