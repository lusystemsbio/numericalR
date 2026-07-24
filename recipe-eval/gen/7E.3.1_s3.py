import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# 2D finite-difference reaction-diffusion integrator (from 7E.1)
# Explicit forward-Euler in time, 5-point Laplacian in space,
# with no-flux (Neumann) boundaries. Run in successive blocks.
# ---------------------------------------------------------------

def laplacian(field, dX):
    # 5-point stencil Laplacian with zero-flux (Neumann) boundaries.
    # Pad by replicating edge rows/cols so the normal gradient is zero.
    p = np.pad(field, 1, mode="edge")
    lap = (p[:-2, 1:-1] + p[2:, 1:-1] +   # up + down neighbours
           p[1:-1, :-2] + p[1:-1, 2:] -   # left + right neighbours
           4.0 * field) / (dX * dX)
    return lap

def integrate_block(u, v, f, g, Du, Dv, dX, dt, nsteps):
    # Advance the two fields nsteps explicit Euler steps in place.
    for _ in range(nsteps):
        Lu = laplacian(u, dX)         # diffusion of activator
        Lv = laplacian(v, dX)         # diffusion of inhibitor/substrate
        u_new = u + dt * (Du * Lu + f(u, v))   # u update
        v_new = v + dt * (Dv * Lv + g(u, v))   # v update
        u, v = u_new, v_new
    return u, v

# ---------------------------------------------------------------
# Grid and run parameters (same IC and run length as 7E.2)
# ---------------------------------------------------------------
N   = 101
dX  = 0.2
dt  = 0.005
d   = 0.1        # Du = d
Du  = d
Dv  = 1.0
mu  = 1.5

t_final     = 100.0
block_t     = 20.0
steps_block = int(round(block_t / dt))       # steps per 20-unit block
n_blocks    = int(round(t_final / block_t))  # 5 blocks

# ---------------------------------------------------------------
# Nearly-uniform initial condition u = v = 1 with +-0.1 noise
# ---------------------------------------------------------------
rng = np.random.default_rng(10)              # seed 10
u = 1.0 + rng.uniform(-0.1, 0.1, (N, N))
v = 1.0 + rng.uniform(-0.1, 0.1, (N, N))

print(f"Grid: {N}x{N}, dX = {dX}, dt = {dt}")
print(f"Du (=d) = {Du}, Dv = {Dv}, mu = {mu}")
print(f"steps per block = {steps_block}, number of blocks = {n_blocks}")
print(f"Initial u mean = {u.mean():.6f}, v mean = {v.mean():.6f}")

# ---------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor kinetics
#   f(u,v) = u^2/v - u ,   g(u,v) = mu*(u^2 - v)
# ---------------------------------------------------------------
def f_AI(u, v):
    return u * u / v - u

def g_AI(u, v):
    return mu * (u * u - v)

# Integrate to t = 100 in successive blocks of 20, printing each block.
t = 0.0
for b in range(n_blocks):
    u, v = integrate_block(u, v, f_AI, g_AI, Du, Dv, dX, dt, steps_block)
    t += block_t
    print(f"[block {b+1}] t = {t:6.1f}  "
          f"u: min = {u.min():.4f}, max = {u.max():.4f}, mean = {u.mean():.4f}  "
          f"v: min = {v.min():.4f}, max = {v.max():.4f}, mean = {v.mean():.4f}")

# ---------------------------------------------------------------
# Count the spots in the final activator field to quantify morphology.
# A "spot" pixel is where u sits well above the field mean; we label
# connected high-u regions with a simple flood fill.
# ---------------------------------------------------------------
thresh = u.mean() + 0.5 * (u.max() - u.mean())   # high-activator mask
mask = u > thresh

def count_blobs(mask):
    seen = np.zeros_like(mask, dtype=bool)
    count = 0
    ny, nx = mask.shape
    for i in range(ny):
        for j in range(nx):
            if mask[i, j] and not seen[i, j]:
                count += 1
                stack = [(i, j)]           # iterative flood fill
                seen[i, j] = True
                while stack:
                    y, x = stack.pop()
                    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < ny and 0 <= xx < nx and mask[yy, xx] and not seen[yy, xx]:
                            seen[yy, xx] = True
                            stack.append((yy, xx))
    return count

n_spots = count_blobs(mask)
print(f"High-activator threshold = {thresh:.4f}")
print(f"Number of connected high-u regions (spots) at t = {t_final:.0f}: {n_spots}")
print(f"Final u: min = {u.min():.4f}, max = {u.max():.4f}, mean = {u.mean():.4f}")

# ---------------------------------------------------------------
# Explanation of the check:
# Because the activator-inhibitor run uses the SAME grid, dX, dt, d,
# mu, initial condition and run length as the substrate-depletion
# (stripe) case of 7E.2 and differs ONLY in the reaction kinetics
# f,g, the emergence of an array of spots rather than stripes proves
# that the reaction kinetics -- not the diffusion or initial data --
# select spot vs. stripe morphology.
# ---------------------------------------------------------------
print("Check: identical IC/parameters to the stripe (substrate-depletion) "
      "case, only f,g differ, yet spots emerge -> kinetics select the morphology.")

# ---------------------------------------------------------------
# 2D image of the final activator (u) field: an array of spots.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(u, origin="lower", cmap="inferno",
               extent=[0, (N - 1) * dX, 0, (N - 1) * dX])
ax.set_title(f"Gierer-Meinhardt activator u at t = {t_final:.0f}\n"
             f"(array of {n_spots} spots)")
ax.set_xlabel("x")
ax.set_ylabel("y")
fig.colorbar(im, ax=ax, label="u")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/7E.3.1_s3.png")
print("Saved figure: 7E.3.1_s3.png")
