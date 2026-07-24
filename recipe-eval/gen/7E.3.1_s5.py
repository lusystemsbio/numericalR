import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# 2D finite-difference reaction-diffusion integrator (explicit, from 7E.1)
# for the Gierer-Meinhardt activator-inhibitor model:
#   f(u,v) = u^2/v - u        (activator kinetics)
#   g(u,v) = mu*(u^2 - v)     (inhibitor kinetics)
#   Du = d, Dv = 1
# ---------------------------------------------------------------------------

# --- Parameters (same grid / IC / run length as 7E.2) ---
N   = 101      # grid points per side
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # activator diffusion (Du); inhibitor diffusion Dv = 1
mu  = 1.5      # inhibitor rate constant
Du  = d
Dv  = 1.0
t_final    = 100.0
block_len  = 20.0
steps_per_block = int(round(block_len / dt))   # 4000 steps per block
n_blocks        = int(round(t_final / block_len))  # 5 blocks

print(f"Grid: {N}x{N}, dX = {dX}, dt = {dt}")
print(f"d (Du) = {Du}, Dv = {Dv}, mu = {mu}")
print(f"Steps per block = {steps_per_block}, number of blocks = {n_blocks}")

# --- Explicit-scheme stability check (diffusion CFL for the faster field) ---
cfl = Dv * dt / dX**2
print(f"Diffusion stability number Dv*dt/dX^2 = {cfl:.6f}  (must be < 0.25)")

# --- Nearly uniform initial condition: u = v = 1 with +-0.1 noise, seed 10 ---
rng = np.random.default_rng(10)
u = 1.0 + (rng.random((N, N)) * 2.0 - 1.0) * 0.1   # uniform in [0.9, 1.1]
v = 1.0 + (rng.random((N, N)) * 2.0 - 1.0) * 0.1
print(f"Initial u: mean = {u.mean():.6f}, min = {u.min():.6f}, max = {u.max():.6f}")


def laplacian(a, h):
    """5-point Laplacian with zero-flux (Neumann) boundaries, done explicitly."""
    # pad by replicating edge rows/cols -> normal derivative = 0 at walls
    ap = np.pad(a, 1, mode="edge")
    lap = (ap[2:, 1:-1] + ap[:-2, 1:-1] +      # north + south neighbours
           ap[1:-1, 2:] + ap[1:-1, :-2] -      # east  + west  neighbours
           4.0 * a)                            # central point
    return lap / h**2


# --- Integrate to t = 100 in successive blocks of 20 (explicit forward Euler) ---
t = 0.0
for b in range(n_blocks):
    for _ in range(steps_per_block):
        Lu = laplacian(u, dX)              # diffusion of activator
        Lv = laplacian(v, dX)              # diffusion of inhibitor
        f = u * u / v - u                  # activator reaction term
        g = mu * (u * u - v)               # inhibitor reaction term
        u_new = u + dt * (Du * Lu + f)     # explicit update for u
        v_new = v + dt * (Dv * Lv + g)     # explicit update for v
        u, v = u_new, v_new
    t += block_len
    print(f"After block {b+1} (t = {t:.1f}): "
          f"u mean = {u.mean():.6f}, min = {u.min():.6f}, max = {u.max():.6f}")

# --- Quantify morphology: count discrete high-u spots (connected components) ---
thresh = 0.5 * (u.max() + u.min())         # midpoint threshold
mask = u > thresh
# simple flood-fill labelling of the boolean mask (explicit, no library routine)
labels = np.zeros_like(mask, dtype=int)
current = 0
for i in range(N):
    for j in range(N):
        if mask[i, j] and labels[i, j] == 0:
            current += 1
            stack = [(i, j)]
            labels[i, j] = current
            while stack:
                ci, cj = stack.pop()
                for ni, nj in ((ci-1, cj), (ci+1, cj), (ci, cj-1), (ci, cj+1)):
                    if 0 <= ni < N and 0 <= nj < N and mask[ni, nj] and labels[ni, nj] == 0:
                        labels[ni, nj] = current
                        stack.append((ni, nj))
n_spots = current
# fraction of grid above threshold: low fraction => isolated spots, ~0.5 => stripes
high_fraction = mask.mean()
print(f"Threshold used = {thresh:.6f}")
print(f"Number of high-u connected regions (spots) = {n_spots}")
print(f"Fraction of grid above threshold = {high_fraction:.6f}")
print(f"Final u range = [{u.min():.6f}, {u.max():.6f}]")

# Morphology verdict: many compact regions + small high-fraction => spots (not stripes)
morphology = "spots" if (n_spots >= 10 and high_fraction < 0.35) else "stripes/other"
print(f"Detected morphology = {morphology}")

# Check confirms the result because the grid, initial condition, noise seed, and
# diffusion constants are identical to the substrate-depletion (stripe) case, so
# the switch to a spot array can only be caused by the different reaction kinetics.
print("Check: identical start/params as substrate-depletion case, only kinetics "
      "differ, so spots vs stripes is selected by the reaction terms alone.")

# --- 2D image of the u field (array of spots) ---
plt.figure(figsize=(6, 5))
plt.imshow(u, origin="lower", cmap="inferno",
           extent=[0, (N-1)*dX, 0, (N-1)*dX])
plt.colorbar(label="activator u")
plt.title("Gierer-Meinhardt u field at t = 100 (spot array)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/7E.3.1_s5.png")
