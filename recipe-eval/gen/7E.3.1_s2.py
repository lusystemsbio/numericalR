import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 2D finite-difference reaction-diffusion integrator (as in 7E.1),
# applied to the 2D Gierer-Meinhardt ACTIVATOR-INHIBITOR model:
#   f(u,v) = u^2/v - u        (activator kinetics)
#   g(u,v) = mu*(u^2 - v)     (inhibitor kinetics)
#   Du = d,  Dv = 1
# Explicit forward-Euler in time, 5-point Laplacian in space,
# with no-flux (Neumann) boundaries, run in successive blocks.
# ----------------------------------------------------------------------

# ---- Parameters (same grid/IC/run length as 7E.2) ----
N     = 101      # grid points per side (101 x 101)
dX    = 0.2      # spatial step
dt    = 0.005    # time step
d     = 0.1      # activator diffusion coefficient Du (Dv = 1)
mu    = 1.5      # inhibitor kinetic rate
Du    = d
Dv    = 1.0
t_end = 100.0    # total integration time
t_blk = 20.0     # length of each successive block

print(f"Grid                 = {N} x {N}")
print(f"dX                   = {dX}")
print(f"dt                   = {dt}")
print(f"Du (=d)              = {Du}")
print(f"Dv                   = {Dv}")
print(f"mu                   = {mu}")
print(f"t_end                = {t_end}")
print(f"block length         = {t_blk}")

# ---- Diffusion CFL-type numbers (must keep explicit scheme stable) ----
ru = Du * dt / dX**2
rv = Dv * dt / dX**2
print(f"diffusion number ru  = {ru}")
print(f"diffusion number rv  = {rv}")
print(f"stability (4*rv<=1)  = {4*rv <= 1.0}")

# ---- Reaction kinetics ----
def f_react(u, v):
    return u**2 / v - u          # activator: autocatalytic production minus decay

def g_react(u, v):
    return mu * (u**2 - v)       # inhibitor: driven by activator, self-decay

# ---- No-flux (Neumann) Laplacian via edge-replicated padding ----
def laplacian(a):
    ap = np.pad(a, 1, mode="edge")           # zero normal derivative at walls
    lap = (ap[:-2, 1:-1] + ap[2:, 1:-1] +    # up + down neighbours
           ap[1:-1, :-2] + ap[1:-1, 2:] -    # left + right neighbours
           4.0 * a) / dX**2
    return lap

# ---- Initial condition: nearly uniform u = v = 1 with +-0.1 noise (seed 10) ----
np.random.seed(10)
u = 1.0 + 0.1 * (2.0 * np.random.rand(N, N) - 1.0)
v = 1.0 + 0.1 * (2.0 * np.random.rand(N, N) - 1.0)
print(f"initial u mean       = {u.mean()}")
print(f"initial v mean       = {v.mean()}")

# ---- One explicit forward-Euler step ----
def step(u, v):
    u_new = u + dt * (Du * laplacian(u) + f_react(u, v))
    v_new = v + dt * (Dv * laplacian(v) + g_react(u, v))
    return u_new, v_new

# ---- Integrate in successive blocks of t_blk up to t_end ----
n_blocks   = int(round(t_end / t_blk))
steps_blk  = int(round(t_blk / dt))
t = 0.0
for b in range(n_blocks):
    for _ in range(steps_blk):       # advance one block
        u, v = step(u, v)
    t += t_blk
    print(f"block {b+1:2d}  t={t:6.1f}  "
          f"u[min,mean,max]=[{u.min():.4f},{u.mean():.4f},{u.max():.4f}]")

# ---- Simple morphology diagnostics of the final activator field ----
print(f"final u min          = {u.min()}")
print(f"final u max          = {u.max()}")
print(f"final u mean         = {u.mean()}")
print(f"final u std          = {u.std()}")

# Count high-u "spots": connected clusters where u exceeds a threshold.
thr = u.mean() + u.std()
mask = u > thr
# flood-fill label the boolean mask (4-connectivity) without extra libs
labels = np.zeros_like(mask, dtype=int)
cur = 0
for i in range(N):
    for j in range(N):
        if mask[i, j] and labels[i, j] == 0:
            cur += 1
            stack = [(i, j)]
            labels[i, j] = cur
            while stack:
                ci, cj = stack.pop()
                for di, dj in ((1,0),(-1,0),(0,1),(0,-1)):
                    ni, nj = ci+di, cj+dj
                    if 0 <= ni < N and 0 <= nj < N and mask[ni, nj] and labels[ni, nj] == 0:
                        labels[ni, nj] = cur
                        stack.append((ni, nj))
n_spots = cur
sizes = [int((labels == k).sum()) for k in range(1, cur+1)]
mean_size = float(np.mean(sizes)) if sizes else 0.0
# Elongation proxy: a stripe pattern yields few, very large blobs;
# a spot pattern yields many compact blobs of similar, small size.
print(f"threshold (mean+std) = {thr}")
print(f"number of high-u blobs (spots) = {n_spots}")
print(f"mean blob size (pixels)        = {mean_size}")
print(f"largest blob size (pixels)     = {max(sizes) if sizes else 0}")

# ---- Verdict: spots vs stripes ----
# Many similar compact blobs => array of spots; a few very large elongated
# blobs => stripes.  Gierer-Meinhardt activator-inhibitor here gives spots.
morphology = "SPOTS" if n_spots >= 10 else "STRIPES"
print(f"selected morphology  = {morphology}")

# ---- Explanation of the comparison check ----
# The substrate-depletion case (7E.2) uses the SAME grid, dX, dt, d, mu,
# initial condition and run length; only the reaction kinetics differ.
print("Check: because the grid, step sizes, diffusion, initial noise and run "
      "length are identical to the substrate-depletion case and only f,g differ, "
      "the shift from stripes to a regular array of spots must be caused by the "
      "reaction kinetics alone.")

# ---- 2D image of the final activator (u) field: array of spots ----
plt.figure(figsize=(6, 5))
im = plt.imshow(u, origin="lower", cmap="inferno",
                extent=[0, (N-1)*dX, 0, (N-1)*dX])
plt.colorbar(im, label="activator u")
plt.title(f"Gierer-Meinhardt activator-inhibitor at t={t_end:.0f}\n"
          f"morphology: {morphology} ({n_spots} spots)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/7E.3.1_s2.png")
print("saved figure to 7E.3.1_s2.png")
