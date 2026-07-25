import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# Kessler-Levine model of spiral cAMP waves in a Dictyostelium colony.
#   Field:  dc/dt = a^2 * (d2c/dX2 + d2c/dY2) - k*c + s
#   Cells:  discrete state machine  inactive -> excited -> refractory -> inactive
#           an excited cell secretes cAMP (source s); an inactive cell fires
#           when the local cAMP exceeds the threshold c_T.
# Everything (Laplacian, time stepping, state machine) is done by hand.
# -------------------------------------------------------------------------

# ---- parameters ----
N       = 101      # grid points per side (101 x 101)
dx      = 1.0      # grid spacing
a       = 1.0      # diffusion length scale  (diffusion coeff = a^2)
frac    = 0.15     # fraction of grid sites that carry a cell
c_T     = 1.0      # firing threshold
dc      = 300.0    # total cAMP secreted during one excitation
t_e     = 2.0      # excited (secreting) duration
t_r     = 20.0     # refractory duration
k       = 0.5      # cAMP degradation rate
dt      = 0.01     # time step
t_end   = 150.0    # total simulated time
sec_rate = dc / t_e            # secretion rate so that rate*t_e = dc
nsteps   = int(round(t_end/dt))

# FTCS stability check for 2D diffusion (must be <= 1/4)
print(f"diffusion_stability_number = {dt*a**2/dx**2*4:.4f}  (needs <= 0.25 for stability)")

rng = np.random.default_rng(0)

# ---- geometry: integer row (Y) and column (X) indices ----
rows, cols = np.indices((N, N))          # rows = Y index, cols = X index
cy = cx = N // 2                         # colony centre

# ---- place cells on a random subset of grid sites ----
mask = rng.random((N, N)) < frac         # True where a cell sits
n_cells = int(mask.sum())
print(f"number_of_cells = {n_cells}")

# ---- cell state arrays ( -1 = no cell, 0 = inactive, 1 = excited, 2 = refractory ) ----
state = np.full((N, N), -1, dtype=int)
state[mask] = 0
timer = np.zeros((N, N))                 # time spent in current excited/refractory phase

# ---- cAMP field ----
c = np.zeros((N, N))

# ---- spiral seed (phase-singularity / broken-front initialisation) ----
# Put supra-threshold cAMP in the right half -> a wavefront along X = cx.
c[:, cx:] = 5.0
# Make the upper half of cells refractory so the front has a FREE END at the
# centre; a graded refractory timer staggers recovery so the free end curls
# into a rotating spiral instead of relaxing to a plane wave.
top = mask & (rows > cy)
state[top] = 2
timer[top] = t_r * (1.0 - cols[top] / (N - 1))   # graded 0..t_r across X

# ---- diagnostics collected during the run ----
angle_times, angles = [], []             # centroid angle of excited cells (rotation check)

def excited_centroid_angle():
    exc = (state == 1)
    if exc.sum() == 0:
        return None
    my = rows[exc].mean(); mx = cols[exc].mean()
    return np.arctan2(my - cy, mx - cx)

# -------------------------------------------------------------------------
# time integration (explicit FTCS diffusion + explicit state machine)
# -------------------------------------------------------------------------
for step in range(nsteps):
    # --- Laplacian with no-flux (Neumann) boundaries via edge padding ---
    cp = np.pad(c, 1, mode="edge")       # replicate edge => zero normal gradient
    lap = (cp[2:, 1:-1] + cp[:-2, 1:-1] +
           cp[1:-1, 2:] + cp[1:-1, :-2] - 4.0*c) / dx**2

    # --- source: excited cells secrete cAMP ---
    s = np.zeros((N, N))
    s[state == 1] = sec_rate

    # --- explicit Euler update of the cAMP field ---
    c = c + dt*(a**2*lap - k*c + s)
    np.maximum(c, 0.0, out=c)            # cAMP concentration stays non-negative

    # --- vectorized cell state-machine update ---
    excited     = (state == 1)
    refractory  = (state == 2)
    inactive    = (state == 0)
    timer[excited | refractory] += dt    # advance phase clocks

    # excited -> refractory once it has secreted for t_e
    to_ref = excited & (timer >= t_e)
    state[to_ref] = 2; timer[to_ref] = 0.0

    # refractory -> inactive once it has recovered for t_r
    to_inact = refractory & (timer >= t_r)
    state[to_inact] = 0; timer[to_inact] = 0.0

    # inactive -> excited when local cAMP exceeds the threshold
    fire = inactive & (c > c_T)
    state[fire] = 1; timer[fire] = 0.0

    # --- record spiral-rotation diagnostic every 1000 steps ---
    if step % 1000 == 0:
        ang = excited_centroid_angle()
        if ang is not None:
            angle_times.append(step*dt); angles.append(ang)

# -------------------------------------------------------------------------
# CHECK: (1) excited cells sit on the high-cAMP crest,
#        (2) the excited region rotates about the core.
# -------------------------------------------------------------------------
exc = (state == 1)
mean_c_all       = c[mask].mean()
mean_c_excited   = c[exc].mean() if exc.sum() else float("nan")
print(f"mean_cAMP_all_cells      = {mean_c_all:.4f}")
print(f"mean_cAMP_excited_cells  = {mean_c_excited:.4f}")
print(f"crest_ratio_excited_over_all = {mean_c_excited/mean_c_all:.4f}")
print(f"num_excited_cells_final  = {int(exc.sum())}")
print(f"num_refractory_cells_final = {int((state==2).sum())}")
print(f"num_inactive_cells_final = {int((state==0).sum())}")

# unwrap centroid angle and report total rotation (evidence of a rotating spiral)
if len(angles) > 1:
    unwrapped = np.unwrap(np.array(angles))
    total_rotation = unwrapped[-1] - unwrapped[0]
    print(f"total_excited_centroid_rotation_rad = {total_rotation:.4f}")
    print(f"total_excited_centroid_rotation_turns = {total_rotation/(2*np.pi):.4f}")
print(f"cAMP_field_min = {c.min():.4f}")
print(f"cAMP_field_max = {c.max():.4f}")

# -------------------------------------------------------------------------
# snapshot: cAMP field + cell states
# -------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))

im = ax[0].imshow(c, origin="lower", cmap="magma")
ax[0].set_title(f"cAMP field at t = {t_end:g} (rotating spiral)")
ax[0].set_xlabel("X"); ax[0].set_ylabel("Y")
fig.colorbar(im, ax=ax[0], label="cAMP concentration")

# right panel: cAMP crest as background + coloured cell states on top
ax[1].imshow(c, origin="lower", cmap="Greys", alpha=0.6)
for st, color, lab in [(0, "tab:blue", "inactive"),
                       (2, "tab:orange", "refractory"),
                       (1, "red", "excited")]:
    sel = (state == st)
    ax[1].scatter(cols[sel], rows[sel], s=6, c=color, label=lab)
ax[1].set_title("cell states (excited ride the crest)")
ax[1].set_xlabel("X"); ax[1].set_ylabel("Y")
ax[1].legend(loc="upper right", framealpha=0.9, markerscale=2)
ax[1].set_xlim(-0.5, N-0.5); ax[1].set_ylim(-0.5, N-0.5)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7D.2.1_s1.png", dpi=130)

# -------------------------------------------------------------------------
# One-sentence explanation of the check:
# The check confirms the result because a persistently rotating excited-cell
# centroid (nonzero net rotation) combined with excited cells whose mean cAMP
# far exceeds the colony average (crest_ratio >> 1) is the defining signature
# of a self-organized rotating spiral wave, distinguishing it from a static
# or purely outward-radial pattern.
# -------------------------------------------------------------------------
print("check_explanation = A net centroid rotation plus excited cells sitting "
      "at cAMP far above the mean shows a persistent rotating high-cAMP spiral "
      "with excited cells riding its crest, not a static or radial pattern.")
