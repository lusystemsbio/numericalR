import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Conway's Game of Life on a 10x10 torus (periodic/wraparound edges)
# ---------------------------------------------------------------

N = 10  # grid size

def step(grid):
    """Advance the grid one generation, implemented explicitly."""
    n = grid.shape[0]
    # neighbor counts for every cell
    counts = np.zeros_like(grid)
    # sum over the eight neighbor offsets, using np.roll for wraparound
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue  # skip the cell itself
            # roll shifts the grid so neighbor lands on the cell; edges wrap
            counts += np.roll(np.roll(grid, di, axis=0), dj, axis=1)
    # apply the three rules to all cells at once
    new = np.zeros_like(grid)
    # live cell with 2 or 3 live neighbors survives
    new[(grid == 1) & ((counts == 2) | (counts == 3))] = 1
    # dead cell with exactly 3 live neighbors becomes alive
    new[(grid == 0) & (counts == 3)] = 1
    # every other cell is dead (already 0 by default)
    return new

# ---- seed a glider ----
grid = np.zeros((N, N), dtype=int)
glider_cells = [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
for (r, c) in glider_cells:
    grid[r, c] = 1

# ---- run and record several generations ----
generations = 20
history = [grid.copy()]
g = grid.copy()
for _ in range(generations):
    g = step(g)
    history.append(g.copy())

# print a few successive generations as text so the glider is visible
for gen in range(0, 9):
    print(f"Generation {gen}:")
    for row in history[gen]:
        print("".join("#" if v else "." for v in row))
    print()

# ---------------------------------------------------------------
# Check: glider keeps its 5-cell shape and drifts (1,1) every 4 gens
# ---------------------------------------------------------------

def cells_of(grid):
    """Return the set of live-cell coordinates."""
    return set(map(tuple, np.argwhere(grid == 1)))

def normalize(cells, n):
    """Shift a shape so its top-left bounding-box corner is at origin.
    Coordinates are taken mod n so a wrapped glider still matches."""
    rs = [r for r, c in cells]
    cs = [c for r, c in cells]
    r0, c0 = min(rs), min(cs)
    return frozenset(((r - r0) % n, (c - c0) % n) for r, c in cells)

base_shape = normalize(cells_of(history[0]), N)
print("Live-cell count each generation:")
for gen in range(generations + 1):
    print(f"  gen {gen:2d}: {int(history[gen].sum())} live cells")
print()

# glider period is 4; verify shape preserved and diagonal drift of (1,1)
period = 4
def centroid_mod(cells, n):
    """Sum of coordinates mod n; used to measure translation."""
    rs = sum(r for r, c in cells)
    cs = sum(c for r, c in cells)
    return rs, cs

print("Shape/drift check (comparing gen k to gen k+4):")
all_ok = True
for gen in range(generations + 1 - period):
    a = cells_of(history[gen])
    b = cells_of(history[gen + period])
    same_shape = (normalize(a, N) == normalize(b, N)) and (len(a) == 5) and (len(b) == 5)
    # measure the translation of every cell (mod N to handle wrap)
    ra, ca = centroid_mod(a, N)
    rb, cb = centroid_mod(b, N)
    drow = ((rb - ra) // 5) % N  # 5 cells, so average shift = sum/5
    dcol = ((cb - ca) // 5) % N
    drift_ok = (drow == 1 and dcol == 1)
    ok = same_shape and drift_ok
    all_ok = all_ok and ok
    print(f"  gen {gen:2d} -> {gen+period:2d}: shape preserved={same_shape}, "
          f"diagonal drift=({drow},{dcol}), pass={ok}")
print()
print(f"Glider maintains 5-cell shape and drifts (1,1) every 4 generations: {all_ok}")

# confirm wraparound: after 40 generations (10 diagonal steps) glider returns home
g = grid.copy()
for _ in range(40):
    g = step(g)
returned = normalize(cells_of(g), N) == base_shape and cells_of(g) == cells_of(history[0])
print(f"After 40 generations glider wraps fully around and returns to start: {returned}")
print()
print("Why this check confirms the result: a glider is defined precisely by its "
      "5-cell shape reappearing translated by (1,1) after its period of 4, so "
      "matching that invariant (with wraparound) is exactly what identifies the "
      "pattern as a correctly simulated, edge-wrapping glider.")

# ---------------------------------------------------------------
# Plot successive generations to show the glider moving
# ---------------------------------------------------------------
fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for i, ax in enumerate(axes.flat):
    ax.imshow(history[i], cmap="binary", vmin=0, vmax=1)
    ax.set_title(f"Gen {i}")
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Conway's Game of Life: glider drifting on a 10x10 torus")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.5.1_s3.png")
