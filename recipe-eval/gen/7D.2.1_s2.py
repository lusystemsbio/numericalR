import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------
# Kessler-Levine model of cAMP spiral waves in a Dictyostelium colony.
# Continuous cAMP field  dc/dt = a^2 * Laplacian(c) - k*c + s
# coupled to discrete excitable cells cycling:
#   inactive -> excited (fires/secretes when local c > c_T) -> refractory -> inactive
# -------------------------------------------------------------------

# ---- Parameters ----
N       = 101      # grid size (N x N)
frac    = 0.15     # fraction of grid sites that are cells
c_T     = 1.0      # excitation threshold on local cAMP
dc      = 300.0    # total cAMP secreted by a firing cell
t_e     = 2.0      # excited (firing) duration
t_r     = 20.0     # refractory (recovery) duration
k       = 0.5      # cAMP degradation rate
a2      = 1.0      # diffusion coefficient a^2 (a=1)
dx      = 1.0      # grid spacing
dt      = 0.01     # time step
t_end   = 150.0    # final time
nsteps  = int(round(t_end / dt))

s_rate  = dc / t_e  # secretion RATE so that total over t_e equals dc

rng = np.random.default_rng(0)

# ---- Fields and cell arrays (full-grid, masked by is_cell) ----
c      = np.zeros((N, N))              # cAMP concentration field
is_cell = rng.random((N, N)) < frac    # which sites hold a cell
state  = np.zeros((N, N), dtype=int)   # 0=inactive, 1=excited, 2=refractory
timer  = np.zeros((N, N))              # countdown timer for current state

# ---- Seed a BROKEN planar wave to nucleate a spiral -------------------
# A wave front lives only in the upper half (Y>=50); its free end at Y=50
# curls around into a rotating spiral (classic excitable-media nucleation).
X, Y = np.meshgrid(np.arange(N), np.arange(N), indexing='xy')
front = is_cell & (X >= 40) & (X <= 50) & (Y >= 50)   # excited wavefront
tail  = is_cell & (X < 40)            & (Y >= 50)      # refractory tail behind it

state[front] = 1
timer[front] = rng.uniform(0.1, t_e, size=front.sum())   # partway through firing
c[front]     = 5.0                                        # kickstart the field

state[tail]  = 2
timer[tail]  = rng.uniform(0.1, t_r, size=tail.sum())     # recovering behind front


def laplacian(f):
    # 5-point Laplacian with no-flux (reflecting) boundaries via edge padding
    fp = np.pad(f, 1, mode='edge')
    lap = (fp[2:, 1:-1] + fp[:-2, 1:-1] +
           fp[1:-1, 2:] + fp[1:-1, :-2] - 4.0 * fp[1:-1, 1:-1]) / (dx * dx)
    return lap


# ---- Rotation tracking: angle of the excited-cell centroid vs grid center ----
cx = cy = (N - 1) / 2.0
track_t, track_ang = [], []
record_every = int(round(1.0 / dt))   # sample once per unit time

# ---- Explicit time integration loop ----
for step in range(nsteps):
    # 1) build source array: only currently-excited cells secrete cAMP
    source = np.where(state == 1, s_rate, 0.0)

    # 2) explicit Euler update of the cAMP field
    c += dt * (a2 * laplacian(c) - k * c + source)

    # 3) vectorized cell state-machine update
    # 3a) excited cells count down; when done -> refractory
    exc = (state == 1)
    timer[exc] -= dt
    done_e = exc & (timer <= 0)
    state[done_e] = 2
    timer[done_e] = t_r

    # 3b) refractory cells count down; when done -> inactive
    ref = (state == 2)
    timer[ref] -= dt
    done_r = ref & (timer <= 0)
    state[done_r] = 0
    timer[done_r] = 0.0

    # 3c) inactive cells above threshold -> excited (fire)
    fire = is_cell & (state == 0) & (c > c_T)
    state[fire] = 1
    timer[fire] = t_e

    # 4) record excited-cell centroid angle for rotation check
    if (step + 1) % record_every == 0:
        ex_mask = (state == 1)
        if ex_mask.sum() > 0:
            ex_y, ex_x = np.nonzero(ex_mask)
            ang = np.arctan2(ex_y.mean() - cy, ex_x.mean() - cx)
            track_t.append((step + 1) * dt)
            track_ang.append(ang)

# -------------------------------------------------------------------
# Analysis / verification
# -------------------------------------------------------------------
excited_mask   = (state == 1)
refract_mask   = (state == 2)
n_excited      = int(excited_mask.sum())
n_refractory   = int(refract_mask.sum())
mean_c_all     = c.mean()
mean_c_excited = c[excited_mask].mean() if n_excited > 0 else float('nan')
crest_ratio    = mean_c_excited / mean_c_all if mean_c_all != 0 else float('nan')

# Rotation: unwrap centroid angle over the final part of the run and sum winding
track_t   = np.array(track_t)
track_ang = np.array(track_ang)
tail_sel  = track_t >= (t_end - 60.0)
unwrapped = np.unwrap(track_ang[tail_sel])
total_rotation_deg = np.degrees(unwrapped[-1] - unwrapped[0]) if unwrapped.size > 1 else 0.0

print(f"Grid size:                              {N} x {N}")
print(f"Number of cells:                        {int(is_cell.sum())}")
print(f"Final excited cells:                    {n_excited}")
print(f"Final refractory cells:                 {n_refractory}")
print(f"Mean cAMP over whole field:             {mean_c_all:.4f}")
print(f"Mean cAMP at excited cell locations:    {mean_c_excited:.4f}")
print(f"Crest ratio (excited mean / field mean):{crest_ratio:.4f}")
print(f"Max cAMP in field:                      {c.max():.4f}")
print(f"Total centroid rotation (last 60 t):    {total_rotation_deg:.2f} deg")
print(f"Wave self-sustained at t=150:           {n_excited > 0}")
print(f"Excited cells ride the cAMP crest:      {crest_ratio > 1.0}")
print(f"Rotating spiral confirmed:              {(n_excited > 0) and (crest_ratio > 1.0) and (abs(total_rotation_deg) > 90.0)}")

# -------------------------------------------------------------------
# Snapshot figure
# -------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 6))

im = axes[0].imshow(c, origin='lower', cmap='inferno')
axes[0].set_title(f"cAMP field at t = {t_end:g} (rotating spiral wave)")
axes[0].set_xlabel("X"); axes[0].set_ylabel("Y")
fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="cAMP c")

# Cell states over a faint cAMP backdrop: excited on the high-cAMP crest
axes[1].imshow(c, origin='lower', cmap='Greys', alpha=0.6)
ry, rx = np.nonzero(refract_mask)
ey, ex = np.nonzero(excited_mask)
axes[1].scatter(rx, ry, s=6, c='tab:blue', label='refractory')
axes[1].scatter(ex, ey, s=8, c='red',      label='excited (firing)')
axes[1].set_title("Cell states on the cAMP field")
axes[1].set_xlabel("X"); axes[1].set_ylabel("Y")
axes[1].set_xlim(0, N - 1); axes[1].set_ylim(0, N - 1)
axes[1].legend(loc='upper right', framealpha=0.9)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7D.2.1_s2.png", dpi=130)

# One-sentence explanation of the check:
print("Check rationale: a persistent set of excited cells whose local cAMP exceeds "
      "the field mean (crest_ratio > 1) together with a steadily winding excited-cell "
      "centroid (>90 deg rotation) confirms a self-sustained rotating cAMP spiral with "
      "cells riding the high-cAMP crest.")
