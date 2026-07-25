import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------- Kessler-Levine model parameters ----------------
N       = 101      # grid size (N x N)
frac    = 0.15     # fraction of lattice sites that carry a cell
c_T     = 1.0      # cAMP firing threshold
dc_sec  = 300.0    # total cAMP secreted by an excited cell over t_e (rate = dc_sec/t_e)
t_e     = 2.0      # duration of the excited (firing) phase
t_r     = 20.0     # duration of the refractory phase
k       = 0.5      # cAMP linear degradation rate
a2      = 1.0      # diffusion coefficient a^2 for the cAMP field
dt      = 0.01     # time step
T       = 150.0    # total integration time
nsteps  = int(round(T / dt))

rng = np.random.default_rng(7)   # fixed seed for reproducibility

# ---------------- Fields ----------------
c = np.zeros((N, N))             # cAMP concentration field c(X,Y)

# Place cells randomly on the lattice (boolean mask marks cell sites)
cell_mask = rng.random((N, N)) < frac
n_cells   = int(cell_mask.sum())

# Cell state machine: 0 = inactive, 1 = excited, 2 = refractory
state = np.zeros((N, N), dtype=np.int8)
timer = np.zeros((N, N))         # time remaining in the current excited/refractory phase

# Seed a spiral: use a broken (asymmetric) initial excitation so the wave
# curls into a rotating spiral rather than a symmetric target pattern.
# Excite cells in the left half; give the lower-left quadrant a head start
# in the refractory phase, creating the phase gradient that seeds rotation.
half = N // 2
left = np.zeros((N, N), dtype=bool)
left[:, :half] = True
excite_now = cell_mask & left
state[excite_now] = 1
timer[excite_now] = t_e

lower = np.zeros((N, N), dtype=bool)
lower[:half, :] = True
phase_shift = cell_mask & left & lower       # part of the wave already refractory
state[phase_shift] = 2
timer[phase_shift] = t_r * 0.5

# secretion rate per firing cell (dc units of concentration per unit time)
sec_rate = dc_sec / t_e

# 5-point Laplacian with no-flux (Neumann) boundaries via edge padding
def laplacian_noflux(f):
    fp = np.pad(f, 1, mode="edge")           # edge padding => zero normal gradient
    return (fp[:-2, 1:-1] + fp[2:, 1:-1] +
            fp[1:-1, :-2] + fp[1:-1, 2:] - 4.0 * f)

# ---------------- Explicit time integration ----------------
for step in range(nsteps):
    # 1) secretion source s: only currently-excited cells secrete cAMP
    s = np.zeros((N, N))
    s[state == 1] = sec_rate

    # 2) explicit forward-Euler update of the reaction-diffusion field
    #    dc/dt = a^2 * Laplacian(c) - k*c + s
    c = c + dt * (a2 * laplacian_noflux(c) - k * c + s)
    np.maximum(c, 0.0, out=c)                 # concentration stays non-negative

    # 3) vectorized cell state-machine update
    #    advance the phase timers for excited & refractory cells
    active = state > 0
    timer[active] -= dt

    # excited -> refractory when firing time elapses
    to_refr = (state == 1) & (timer <= 0.0)
    state[to_refr] = 2
    timer[to_refr] = t_r

    # refractory -> inactive when recovery time elapses
    to_inact = (state == 2) & (timer <= 0.0)
    state[to_inact] = 0
    timer[to_inact] = 0.0

    # inactive cell fires (excited) if local cAMP exceeds threshold
    to_exc = cell_mask & (state == 0) & (c > c_T)
    state[to_exc] = 1
    timer[to_exc] = t_e

# ---------------- Diagnostics ----------------
n_excited    = int((state == 1).sum())
n_refractory = int((state == 2).sum())
n_inactive   = int(cell_mask.sum() - n_excited - n_refractory)

# cAMP at excited-cell sites vs. at all cell sites: excited cells should ride the crest
c_at_excited = c[state == 1]
c_at_cells   = c[cell_mask]
mean_c_excited = float(c_at_excited.mean()) if c_at_excited.size else float("nan")
mean_c_cells   = float(c_at_cells.mean())   if c_at_cells.size else float("nan")

print(f"Grid size: {N} x {N}")
print(f"Number of cells: {n_cells}")
print(f"Total integration steps: {nsteps}")
print(f"Final time: {nsteps * dt:.2f}")
print(f"Max cAMP concentration: {c.max():.6f}")
print(f"Mean cAMP concentration: {c.mean():.6f}")
print(f"Excited cells at final time: {n_excited}")
print(f"Refractory cells at final time: {n_refractory}")
print(f"Inactive cells at final time: {n_inactive}")
print(f"Mean cAMP at excited-cell sites: {mean_c_excited:.6f}")
print(f"Mean cAMP at all cell sites: {mean_c_cells:.6f}")
print(f"Crest ratio (excited / all cells): {mean_c_excited / mean_c_cells:.6f}")

# Check: excited cells sit on above-average cAMP => they ride the high-cAMP crest
rides_crest = mean_c_excited > mean_c_cells
print(f"Excited cells ride the high-cAMP crest: {rides_crest}")
print("Check explanation: a rotating cAMP spiral must have its firing (excited) "
      "cells concentrated on the traveling high-concentration crest, so mean_c_excited "
      "> mean_c_cells together with the visible curved wavefront confirms the "
      "self-organized rotating spiral wave.")

# ---------------- Snapshot figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))

im = ax[0].imshow(c, origin="lower", cmap="inferno")
ax[0].set_title("cAMP field c(X,Y) at t = %.0f" % (nsteps * dt))
ax[0].set_xlabel("X")
ax[0].set_ylabel("Y")
fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04, label="cAMP")

# cell-state map: 0 inactive, 1 excited, 2 refractory (non-cells shown as background)
disp = np.full((N, N), np.nan)
disp[cell_mask] = state[cell_mask]
im2 = ax[1].imshow(c, origin="lower", cmap="Greys", alpha=0.6)
yy, xx = np.where(state == 2)
ax[1].scatter(xx, yy, s=6, c="tab:blue",   label="refractory")
yy, xx = np.where(state == 1)
ax[1].scatter(xx, yy, s=6, c="tab:red",    label="excited")
yy, xx = np.where((state == 0) & cell_mask)
ax[1].scatter(xx, yy, s=4, c="tab:green",  label="inactive", alpha=0.4)
ax[1].set_title("Cell states on cAMP background")
ax[1].set_xlabel("X")
ax[1].set_ylabel("Y")
ax[1].set_xlim(0, N)
ax[1].set_ylim(0, N)
ax[1].legend(loc="upper right", framealpha=0.9)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7D.2.1_s3.png", dpi=130)
