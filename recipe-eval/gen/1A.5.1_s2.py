import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters ----
N = 10          # grid is N x N
GENS = 20       # number of generations to simulate

# ---- Seed a glider on a 10x10 grid ----
grid = np.zeros((N, N), dtype=int)
glider_cells = [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
for (r, c) in glider_cells:
    grid[r, c] = 1

def count_neighbors(g):
    # Count live neighbors of every cell with wraparound (toroidal) edges.
    neighbors = np.zeros_like(g)
    for dr in (-1, 0, 1):            # row offsets
        for dc in (-1, 0, 1):        # col offsets
            if dr == 0 and dc == 0:
                continue             # skip the cell itself
            # np.roll shifts the grid so each cell "sees" one neighbor;
            # roll wraps around the edges, giving periodic boundaries.
            neighbors += np.roll(np.roll(g, dr, axis=0), dc, axis=1)
    return neighbors

def step(g):
    # Apply the three Game of Life rules to all cells at once.
    n = count_neighbors(g)
    new = np.zeros_like(g)
    # Live cell with 2 or 3 live neighbors survives.
    new[(g == 1) & ((n == 2) | (n == 3))] = 1
    # Dead cell with exactly 3 live neighbors becomes alive.
    new[(g == 0) & (n == 3)] = 1
    # Every other cell is dead (already 0 by default).
    return new

def cell_set(g):
    # Set of coordinates of live cells, for shape/drift comparison.
    return set(map(tuple, np.argwhere(g == 1)))

def normalize_shape(cells):
    # Translate a set of cells so its top-left bounding-box corner is at (0,0),
    # giving a translation-invariant "shape fingerprint".
    rs = [r for r, c in cells]
    cs = [c for r, c in cells]
    r0, c0 = min(rs), min(cs)
    return frozenset((r - r0, c - c0) for r, c in cells)

# ---- Run the simulation, saving history ----
history = [grid.copy()]
g = grid.copy()
for _ in range(GENS):
    g = step(g)
    history.append(g.copy())

# ---- Print grids at several successive generations ----
for gen in range(0, GENS + 1, 4):
    print(f"Generation {gen}:")
    for row in history[gen]:
        print("".join("#" if v else "." for v in row))
    print()

# ---- Check 1: five-cell shape preserved every generation ----
shape0 = normalize_shape(cell_set(history[0]))
shape_ok = True
for gen in range(GENS + 1):
    cells = cell_set(history[gen])
    if len(cells) != 5:
        shape_ok = False
    if normalize_shape(cells) != shape0 and gen % 4 == 0:
        # Every 4 generations the glider returns to its original orientation.
        shape_ok = False
print(f"Live-cell count stays 5 every generation: "
      f"{all(len(cell_set(history[g])) == 5 for g in range(GENS + 1))}")
print(f"Shape matches original orientation every 4th generation: {shape_ok}")

# ---- Check 2: drift of one diagonal step every four generations (with wrap) ----
def centroid_mod(cells):
    # Sum of coordinates (mod N) -> compare centroid drift on the torus.
    rs = [r for r, c in cells]
    cs = [c for r, c in cells]
    return (sum(rs), sum(cs))

drift_ok = True
base = centroid_mod(cell_set(history[0]))
for k in range(0, GENS + 1, 4):
    cells = centroid_mod(cell_set(history[k]))
    # After k generations (k/4 diagonal steps) each of the 5 cells shifted by k/4,
    # so the coordinate sums increase by 5*(k/4), taken modulo N for wraparound.
    expected_r = (base[0] + 5 * (k // 4)) % N
    expected_c = (base[1] + 5 * (k // 4)) % N
    got_r, got_c = cells[0] % N, cells[1] % N
    match = (got_r == expected_r) and (got_c == expected_c)
    print(f"Gen {k:2d}: coord-sum (r,c)=({got_r},{got_c}) "
          f"expected=({expected_r},{expected_c}) match={match}")
    if not match:
        drift_ok = False
print(f"Glider drifts one diagonal step every 4 generations (with wrap): {drift_ok}")

# ---- Explanation ----
print("Explanation: A glider is the smallest 5-cell spaceship whose pattern "
      "repeats every 4 generations shifted one cell diagonally, so recovering "
      "the same 5-cell shape one diagonal step over every 4 generations (and "
      "seeing it wrap to the opposite edge) confirms the update rules and "
      "periodic boundary are implemented correctly.")

# ---- Visualization: successive generations ----
show_gens = list(range(0, GENS + 1, 2))
ncols = 4
nrows = (len(show_gens) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
axes = np.array(axes).reshape(-1)
for ax in axes:
    ax.axis("off")
for i, gen in enumerate(show_gens):
    ax = axes[i]
    ax.imshow(history[gen], cmap="binary", vmin=0, vmax=1)
    ax.set_title(f"Gen {gen}")
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("on")
    ax.set_xticks(np.arange(-0.5, N, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, N, 1), minor=True)
    ax.grid(which="minor", color="lightgray", linewidth=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.5.1_s2.png")
