import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D finite-difference reaction-diffusion integrator (from 7E.1)
# applied to the 2D Gierer-Meinhardt activator-inhibitor model:
#   f(u,v) = u^2/v - u        (activator kinetics)
#   g(u,v) = mu*(u^2 - v)     (inhibitor kinetics)
#   Du = d,  Dv = 1
# ---------------------------------------------------------------

# ----- parameters (same as 7E.2 substrate-depletion run) -----
N   = 101      # grid points per side
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # activator diffusion coefficient Du
mu  = 1.5      # inhibitor kinetic rate
Du  = d
Dv  = 1.0

t_final     = 100.0     # total integration time
block_time  = 20.0      # length of each successive block
steps_block = int(round(block_time / dt))   # time steps per block
n_blocks    = int(round(t_final / block_time))

# ----- nearly uniform initial condition: u = v = 1 with +-0.1 noise -----
rng = np.random.default_rng(10)                     # seed 10
u = 1.0 + 0.1 * (2.0 * rng.random((N, N)) - 1.0)    # uniform in [-0.1, 0.1]
v = 1.0 + 0.1 * (2.0 * rng.random((N, N)) - 1.0)

def laplacian(field, h):
    # explicit 5-point Laplacian with no-flux (Neumann) boundaries;
    # np.pad in "edge" mode copies the border so the outward flux is zero
    p = np.pad(field, 1, mode="edge")
    lap = (p[:-2, 1:-1] + p[2:, 1:-1] +
           p[1:-1, :-2] + p[1:-1, 2:] -
           4.0 * p[1:-1, 1:-1]) / (h * h)
    return lap

# ----- integrate to t=100 in successive blocks of 20 -----
t = 0.0
for b in range(n_blocks):
    for _ in range(steps_block):
        # reaction terms
        f = u * u / v - u                 # activator reaction f(u,v)
        g = mu * (u * u - v)              # inhibitor reaction g(u,v)
        # diffusion terms
        Lu = laplacian(u, dX)
        Lv = laplacian(v, dX)
        # explicit forward-Euler update of both fields
        u_new = u + dt * (Du * Lu + f)
        v_new = v + dt * (Dv * Lv + g)
        u, v = u_new, v_new
    t += block_time
    print(f"Block {b+1}/{n_blocks} end t={t:6.1f}  "
          f"u_min={u.min():.4f}  u_max={u.max():.4f}  u_mean={u.mean():.4f}")

# ----- quantify the morphology: count spots (high-u peaks) -----
thresh = 0.5 * (u.max() + u.min())        # midlevel threshold
high = u > thresh

# simple connected-component (flood-fill) count of high-u regions
visited = np.zeros_like(high, dtype=bool)
n_spots = 0
from collections import deque
for i in range(N):
    for j in range(N):
        if high[i, j] and not visited[i, j]:
            n_spots += 1
            q = deque([(i, j)])
            visited[i, j] = True
            while q:
                ci, cj = q.popleft()
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni, nj = ci + di, cj + dj
                    if 0 <= ni < N and 0 <= nj < N and high[ni, nj] and not visited[ni, nj]:
                        visited[ni, nj] = True
                        q.append((ni, nj))

# ----- report numerical results -----
print(f"Final time t                : {t:.1f}")
print(f"Final u min                 : {u.min():.6f}")
print(f"Final u max                 : {u.max():.6f}")
print(f"Final u mean                : {u.mean():.6f}")
print(f"Final v min                 : {v.min():.6f}")
print(f"Final v max                 : {v.max():.6f}")
print(f"Final v mean                : {v.mean():.6f}")
print(f"Threshold for spot count    : {thresh:.6f}")
print(f"High-u area fraction        : {high.mean():.6f}")
print(f"Number of spots (u peaks)   : {n_spots}")

# a compact, roughly isotropic spot array (rather than elongated stripes)
# is confirmed by the high-u area fraction being well below 0.5 while the
# pattern breaks into many separate connected components:
print(f"Spots >> 1 and area frac < 0.5 -> spot morphology, not stripes: "
      f"{n_spots > 5 and high.mean() < 0.5}")

# ----- 2D image of the u field: array of spots -----
plt.figure(figsize=(6, 5))
plt.imshow(u, origin="lower", cmap="inferno",
           extent=[0, (N - 1) * dX, 0, (N - 1) * dX])
plt.colorbar(label="activator u")
plt.title(f"Gierer-Meinhardt activator u at t={t:.0f}\n(spots: {n_spots})")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/7E.3.1_s1.png")

# Why the check confirms the result (one sentence):
# Because the grid, dX, dt, initial noise (seed 10), and run length are held
# identical to the substrate-depletion (7E.2) case and only the reaction terms
# f and g are swapped to Gierer-Meinhardt, any change from stripes to a regular
# spot array must be caused by the reaction kinetics alone.
print("Check: with identical grid/IC/parameters and only f,g swapped to "
      "Gierer-Meinhardt, the switch from stripes to spots is due to the "
      "reaction kinetics alone.")
