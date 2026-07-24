import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Game of Life on a torus (periodic boundary) ---------------------------

def count_live_neighbors(grid):
    # Sum the eight shifted copies of the grid.
    # np.roll wraps around the edges, giving the periodic (torus) boundary.
    neighbors = np.zeros_like(grid)
    for dr in (-1, 0, 1):          # row shift
        for dc in (-1, 0, 1):      # column shift
            if dr == 0 and dc == 0:
                continue           # skip the cell itself
            neighbors += np.roll(np.roll(grid, dr, axis=0), dc, axis=1)
    return neighbors

def step(grid):
    n = count_live_neighbors(grid)                # live neighbor count for every cell
    survive = (grid == 1) & ((n == 2) | (n == 3)) # live cell with 2 or 3 neighbors survives
    born = (grid == 0) & (n == 3)                 # dead cell with exactly 3 neighbors is born
    return (survive | born).astype(int)           # everything else is dead

# --- Seed a 10x10 grid with a glider ---------------------------------------

N = 10
grid = np.zeros((N, N), dtype=int)
glider = [(1, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
for r, c in glider:
    grid[r, c] = 1

# --- Run and record several generations ------------------------------------

GENERATIONS = 40
history = [grid.copy()]
g = grid.copy()
for _ in range(GENERATIONS):
    g = step(g)
    history.append(g.copy())

# --- Print grids at successive generations ---------------------------------

def print_grid(gen, g):
    print(f"Generation {gen}: live cells = {int(g.sum())}, "
          f"center of mass = {tuple(round(x, 2) for x in np.argwhere(g == 1).mean(axis=0))}")
    for row in g:
        print("".join("#" if c else "." for c in row))
    print()

for gen in range(9):
    print_grid(gen, history[gen])

# --- Check: shape preserved, drifts (1,1) every 4 generations, wraps -------
# Track center of mass modulo N; expected drift is +1 row and +1 col each 4 gens.
print("=== Glider tracking check ===")
com0 = np.argwhere(history[0] == 1).mean(axis=0)
all_five = True
drift_ok = True
for gen in range(0, GENERATIONS + 1, 4):
    g = history[gen]
    live = int(g.sum())
    com = np.argwhere(g == 1).mean(axis=0)
    # expected wrapped shift of the center of mass on the torus
    expected = ((com0 + gen / 4) % N)
    diff = (com - expected + N / 2) % N - N / 2   # signed difference accounting for wrap
    ok = np.allclose(diff, 0, atol=1e-9)
    if live != 5:
        all_five = False
    if not ok:
        drift_ok = False
    print(f"gen {gen:2d}: live = {live}, "
          f"COM = ({com[0]:.2f},{com[1]:.2f}), "
          f"expected (wrapped) = ({expected[0]:.2f},{expected[1]:.2f}), "
          f"matches diagonal drift = {ok}")

print(f"Shape stayed five cells at every 4th generation: {all_five}")
print(f"Diagonal one-step-per-four-generations drift held (with wraparound): {drift_ok}")
print("Why this confirms it: a glider is the smallest spaceship, so a constant "
      "five-cell count plus a steady (+1,+1)-per-4-generations shift of its center "
      "of mass (wrapping at the edges) is exactly the signature of an intact glider "
      "translating across the torus.")

# --- Plot several successive generations -----------------------------------

fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for ax, gen in zip(axes.ravel(), range(8)):
    ax.imshow(history[gen], cmap="binary", vmin=0, vmax=1)
    ax.set_title(f"Gen {gen}")
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Conway's Game of Life: glider on a 10x10 torus")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/1A.5.1_s5.png")
