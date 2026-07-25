import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Kessler-Levine model: excitable cAMP field + discrete Dictyostelium cells
#   dc/dt = a^2 (d2c/dX2 + d2c/dY2) - k c + s
# Cells cycle: inactive -> excited (fire+secrete when local c>c_T) -> refractory -> inactive
# ----------------------------------------------------------------------

# ---- parameters ----
N     = 101      # grid size (N x N)
frac  = 0.15     # fraction of grid sites occupied by a cell
c_T   = 1.0      # firing threshold
dc    = 300.0    # total cAMP secreted by an excited cell over its excited phase
t_e   = 2.0      # duration of the excited (secreting) phase
t_r   = 20.0     # duration of the refractory (recovery) phase
k     = 0.5      # cAMP degradation rate
a2    = 2.0      # diffusion coefficient a^2 (dx=1); dt<=dx^2/(4 a2)=0.125 => stable
dt    = 0.01     # time step
T     = 150.0    # total simulated time
dx    = 1.0

nsteps = int(round(T / dt))
sec_rate = dc / t_e          # constant secretion rate s while a cell is excited
rng = np.random.default_rng(0)

# ---- cAMP field ----
c = np.zeros((N, N))         # cAMP concentration on the grid

# ---- place cells at unique random grid sites ----
ncells = int(frac * N * N)
flat = rng.choice(N * N, size=ncells, replace=False)
cell_r = flat // N           # row (Y) index of each cell
cell_c = flat % N            # col (X) index of each cell

# cell state machine: 0=inactive, 1=excited, 2=refractory ; timer counts time in state
state = np.zeros(ncells, dtype=int)
timer = np.zeros(ncells)

# ---- initial condition: a BROKEN wave -> curls into a spiral ----
# Elevated cAMP block spanning the upper half with a free tip at row 50 (the wave front).
c[50:N, 45:56] = 5.0 * c_T
# Force the region behind/left of the front (upper-left) to be refractory so the
# front cannot propagate leftward; its free tip then rotates -> spiral.
init_refr = (cell_r >= 50) & (cell_c < 45)
state[init_refr] = 2
timer[init_refr] = 0.0

# ---- time integration (explicit) ----
for step in range(nsteps):
    # local cAMP seen by each cell
    local_c = c[cell_r, cell_c]

    # (1) firing: inactive cell with local c above threshold becomes excited
    fire = (state == 0) & (local_c > c_T)
    state[fire] = 1
    timer[fire] = 0.0

    # (2) secretion: every excited cell injects cAMP at its site (source term s*dt)
    exc = (state == 1)
    if np.any(exc):
        np.add.at(c, (cell_r[exc], cell_c[exc]), sec_rate * dt)

    # (3) advance state timers for excited & refractory cells
    timer[state != 0] += dt

    # (4) excited -> refractory after t_e ; refractory -> inactive after t_r
    to_refr = (state == 1) & (timer >= t_e)
    state[to_refr] = 2
    timer[to_refr] = 0.0
    to_inact = (state == 2) & (timer >= t_r)
    state[to_inact] = 0
    timer[to_inact] = 0.0

    # (5) diffusion + degradation of the cAMP field (explicit FD, no-flux BC)
    #     no-flux (Neumann) boundaries via edge-padding => zero gradient at walls
    cp = np.pad(c, 1, mode="edge")
    lap = (cp[2:, 1:-1] + cp[:-2, 1:-1] +
           cp[1:-1, 2:] + cp[1:-1, :-2] - 4.0 * c) / (dx * dx)
    c = c + dt * (a2 * lap - k * c)

# ----------------------------------------------------------------------
# Snapshot: cAMP field with rotating spiral + excited cells, and cell states
# ----------------------------------------------------------------------
exc = (state == 1)
refr = (state == 2)
inact = (state == 0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))

im = ax1.imshow(c, origin="lower", cmap="inferno", vmin=0, vmax=np.percentile(c, 99))
ax1.scatter(cell_c[exc], cell_r[exc], s=6, c="cyan", label="excited cells")
ax1.set_title("cAMP field (spiral wave) + excited cells  t=%.0f" % T)
ax1.set_xlabel("X"); ax1.set_ylabel("Y")
ax1.legend(loc="upper right", fontsize=8)
fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="cAMP c")

ax2.scatter(cell_c[inact], cell_r[inact], s=6, c="lightgray", label="inactive")
ax2.scatter(cell_c[refr],  cell_r[refr],  s=6, c="orange",    label="refractory")
ax2.scatter(cell_c[exc],   cell_r[exc],   s=6, c="red",       label="excited")
ax2.set_xlim(0, N - 1); ax2.set_ylim(0, N - 1); ax2.set_aspect("equal")
ax2.set_title("Cell states"); ax2.set_xlabel("X"); ax2.set_ylabel("Y")
ax2.legend(loc="upper right", fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7D.2.1_s4.png", dpi=130)

# ----------------------------------------------------------------------
# Quantitative check: excited cells ride the high-cAMP crest & activity is sustained
# ----------------------------------------------------------------------
mean_c_global   = float(c.mean())
mean_c_excited  = float(c[cell_r[exc], cell_c[exc]].mean()) if np.any(exc) else float("nan")
mean_c_cells    = float(c[cell_r, cell_c].mean())
n_excited       = int(exc.sum())
n_refractory    = int(refr.sum())
n_inactive      = int(inact.sum())
crest_ratio     = mean_c_excited / mean_c_global if mean_c_global > 0 else float("nan")

print("grid size N x N                 : %d x %d" % (N, N))
print("total cells                     : %d" % ncells)
print("simulated time T                : %.1f  (%d steps, dt=%.3f)" % (T, nsteps, dt))
print("cAMP field min                  : %.6f" % float(c.min()))
print("cAMP field max                  : %.6f" % float(c.max()))
print("cAMP field mean (global)        : %.6f" % mean_c_global)
print("mean cAMP at all cell sites     : %.6f" % mean_c_cells)
print("mean cAMP at EXCITED cell sites : %.6f" % mean_c_excited)
print("crest ratio (excited/global)    : %.4f" % crest_ratio)
print("firing threshold c_T            : %.4f" % c_T)
print("number of excited cells   (t=T) : %d" % n_excited)
print("number of refractory cells(t=T) : %d" % n_refractory)
print("number of inactive cells  (t=T) : %d" % n_inactive)

# Why this check confirms the result (one sentence):
print("EXPLANATION: A rotating cAMP spiral is confirmed because a cohort of cells is "
      "still excited at t=T (self-sustained activity, not a decayed transient) AND their "
      "mean local cAMP far exceeds both the threshold and the global mean (crest ratio >> 1), "
      "i.e. the excited cells sit on the traveling high-cAMP wave crest that keeps re-firing "
      "them as it rotates.")
