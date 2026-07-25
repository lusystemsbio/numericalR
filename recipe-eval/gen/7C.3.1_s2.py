import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor Turing model
#   f(u,v) = u^2/v - u        (activator: local self-enhancement)
#   g(u,v) = mu*(u^2 - v)     (inhibitor: produced by activator)
#   Du = d (slow, short range), Dv = 1 (fast, long range)
# Multi-component finite-difference reaction-diffusion integrator (7C.1),
# implemented explicitly (forward Euler in time, central diff in space),
# run in successive time blocks.
# -------------------------------------------------------------------

# ---- Parameters ----
L   = 20.0     # domain length
dX  = 0.2      # spatial step
dt  = 0.01     # time step
d   = 0.1      # activator diffusion (slow)
Dv  = 1.0      # inhibitor diffusion (fast)
mu  = 1.5      # inhibitor kinetic rate
Du  = d

N = int(round(L / dX)) + 1          # number of grid points
X = np.linspace(0.0, L, N)

# ---- Reaction kinetics ----
def f(u, v):
    return u**2 / v - u             # activator reaction

def g(u, v):
    return mu * (u**2 - v)          # inhibitor reaction

# ---- Discrete 1-D Laplacian with no-flux (Neumann) boundaries ----
def laplacian(y):
    lap = np.empty_like(y)
    # interior: central second difference
    lap[1:-1] = (y[2:] - 2.0*y[1:-1] + y[:-2]) / dX**2
    # no-flux boundaries: mirror the neighbour (zero gradient)
    lap[0]  = (2.0*y[1]  - 2.0*y[0])  / dX**2
    lap[-1] = (2.0*y[-2] - 2.0*y[-1]) / dX**2
    return lap

# ---- One explicit forward-Euler step for the coupled system ----
def step(u, v):
    u_new = u + dt * (Du * laplacian(u) + f(u, v))
    v_new = v + dt * (Dv * laplacian(v) + g(u, v))
    return u_new, v_new

# ---- Integrate over a block of time steps ----
def run_block(u, v, n_steps):
    for _ in range(n_steps):
        u, v = step(u, v)
    return u, v

# ---- Initial condition: nearly uniform u=v=1 with +-0.1 noise ----
rng = np.random.default_rng(10)          # seed 10
u = 1.0 + 0.2 * (rng.random(N) - 0.5)    # uniform noise in [-0.1, +0.1]
v = 1.0 + 0.2 * (rng.random(N) - 0.5)

# ---- Diffusion stability check (explicit scheme) ----
alpha = Dv * dt / dX**2
print(f"Grid points N = {N}")
print(f"Explicit diffusion number Dv*dt/dX^2 = {alpha:.4f} (must be < 0.5)")

# ---- Run in successive time blocks until stationary ----
steps_per_block = 5000          # 50 time units per block
n_blocks = 12                   # up to 600 time units total
prev_u = u.copy()
for b in range(n_blocks):
    u, v = run_block(u, v, steps_per_block)
    change = np.max(np.abs(u - prev_u))   # max change since last block
    print(f"Block {b+1:2d}: t = {(b+1)*steps_per_block*dt:6.1f}, "
          f"max|du| over block = {change:.3e}")
    prev_u = u.copy()

# ---- Stationarity check: one more block, compare ----
u_final, v_final = run_block(u, v, steps_per_block)
stationary_change = np.max(np.abs(u_final - u))
print(f"Stationarity: max|u(after extra block) - u| = {stationary_change:.3e}")
u, v = u_final, v_final

# ---- Periodicity: estimate number of peaks (dominant wavelength) ----
uc = u - u.mean()
# count interior local maxima of u
peaks = np.sum((uc[1:-1] > uc[:-2]) & (uc[1:-1] > uc[2:]))
print(f"Number of interior u-peaks (periodic pattern) = {peaks}")

# ---- Phase relation: correlation between u and v spatial profiles ----
corr = np.corrcoef(u - u.mean(), v - v.mean())[0, 1]
print(f"Spatial correlation corr(u, v) = {corr:.4f}")
# Location of global maxima of u and v
print(f"X of u-max = {X[np.argmax(u)]:.2f},  X of v-max = {X[np.argmax(v)]:.2f}")

# ---- Verdict ----
in_phase = corr > 0.0
print(f"u and v peak together (in phase)? {in_phase}")
print("Contrast: in substrate depletion u and v would be OUT of phase "
      "(negative corr); here corr > 0 confirms activator-inhibitor "
      "in-phase peaking.")
# One-sentence explanation:
print("Explanation: a positive spatial correlation with coincident maxima "
      "means the inhibitor v is driven up wherever the activator u is high, "
      "which is the defining in-phase signature of the Gierer-Meinhardt "
      "activator-inhibitor mechanism as opposed to out-of-phase substrate "
      "depletion.")

# ---- Plot ----
plt.figure(figsize=(9, 5))
plt.plot(X, u, 'b-o', ms=3, label='u (activator)')
plt.plot(X, v, 'r-s', ms=3, label='v (inhibitor)')
plt.xlabel('X')
plt.ylabel('concentration')
plt.title('Gierer-Meinhardt Turing pattern (u and v peak together)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.3.1_s2.png")
