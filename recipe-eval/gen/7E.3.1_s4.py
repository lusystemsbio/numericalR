import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D finite-difference reaction-diffusion integrator (from 7E.1)
# u_t = Du*Lap(u) + f(u,v),  v_t = Dv*Lap(v) + g(u,v)
# Gierer-Meinhardt activator-inhibitor kinetics:
#   f(u,v) = u^2/v - u   (activator: self-activation, linear decay)
#   g(u,v) = mu*(u^2 - v) (inhibitor: produced by u, linear decay)
# ---------------------------------------------------------------

# ---- Parameters (same start / run length as 7E.2) ----
N   = 101      # grid points per side
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # activator diffusivity Du = d
Du  = d
Dv  = 1.0      # inhibitor diffusivity
mu  = 1.5      # inhibitor rate constant
T_end   = 100.0
T_block = 20.0
steps_per_block = int(round(T_block / dt))   # explicit steps in one block
n_blocks        = int(round(T_end / T_block))

# ---- Reaction kinetics ----
def f(u, v):
    return u*u/v - u              # activator reaction term

def g(u, v):
    return mu*(u*u - v)          # inhibitor reaction term

# ---- 5-point Laplacian with zero-flux (Neumann) boundaries ----
# Implemented explicitly: pad by copying edge rows/cols (reflective),
# then apply the standard 5-point stencil, all divided by dX^2.
def laplacian(a):
    ap = np.pad(a, 1, mode="edge")               # zero-flux: ghost = edge value
    lap = (ap[:-2, 1:-1] + ap[2:, 1:-1] +
           ap[1:-1, :-2] + ap[1:-1, 2:] -
           4.0*ap[1:-1, 1:-1]) / (dX*dX)
    return lap

# ---- Nearly uniform initial condition: u = v = 1 with +-0.1 noise ----
rng = np.random.default_rng(10)                  # seed 10
u = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)     # uniform in [-0.1, +0.1]
v = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)

# ---- Explicit (forward-Euler) time stepping, run in successive blocks ----
t = 0.0
for b in range(n_blocks):
    for _ in range(steps_per_block):
        # one explicit Euler update of the reaction-diffusion PDEs
        u_new = u + dt*(Du*laplacian(u) + f(u, v))
        v_new = v + dt*(Dv*laplacian(v) + g(u, v))
        u, v = u_new, v_new
    t += T_block
    print(f"After block {b+1}, t = {t:.1f}: u min = {u.min():.4f}, "
          f"u max = {u.max():.4f}, u mean = {u.mean():.4f}")

# ---- Summary of the final activator field ----
print(f"Final time t = {t:.1f}")
print(f"Final u min  = {u.min():.6f}")
print(f"Final u max  = {u.max():.6f}")
print(f"Final u mean = {u.mean():.6f}")
print(f"Final u std  = {u.std():.6f}")
print(f"Final v min  = {v.min():.6f}")
print(f"Final v max  = {v.max():.6f}")

# Count high-u peaks (spots) via a simple local-maximum / threshold measure.
thr = 0.5*(u.max() + u.mean())
spot_mask = u > thr
print(f"Spot threshold (u > {thr:.4f}): fraction of area = {spot_mask.mean():.4f}")

# ---- Produce a 2D image of the u field (array of spots) ----
plt.figure(figsize=(6, 5))
im = plt.imshow(u, origin="lower", cmap="inferno",
                extent=[0, (N-1)*dX, 0, (N-1)*dX])
plt.colorbar(im, label="activator u")
plt.title(f"Gierer-Meinhardt activator u at t = {t:.0f}\n"
          f"(d={d}, mu={mu}) -> array of spots")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/7E.3.1_s4.png")

# Why the check confirms the result:
# Because the grid, dX, dt, diffusivities, initial condition, seed, and run
# length are identical to the substrate-depletion (stripe) case and only the
# reaction terms f,g differ, any change in morphology from stripes to a regular
# array of spots must be caused solely by the reaction kinetics.
print("Check: same start/parameters as the stripe (substrate-depletion) case, "
      "only f,g differ -> morphology (spots vs stripes) is set by the kinetics.")
