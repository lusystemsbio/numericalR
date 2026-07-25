import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor Turing pattern (1D)
# Multi-component explicit finite-difference reaction-diffusion
# integrator, run in successive time blocks (as in 7C.1).
# ---------------------------------------------------------------

# --- Parameters ---
L      = 20.0      # domain length
dX     = 0.2       # spatial step
dt     = 0.01      # time step
d      = 0.1       # activator diffusion Du (slow)
Dv     = 1.0       # inhibitor diffusion (fast)
Du     = d
mu     = 1.5       # kinetic parameter
N      = int(round(L / dX)) + 1   # number of grid points
X      = np.linspace(0.0, L, N)

# --- Reaction kinetics (Gierer-Meinhardt) ---
def f(u, v):   # activator: local self-enhancement u^2/v minus decay
    return u**2 / v - u
def g(u, v):   # inhibitor: produced by activator (u^2), decays (-v)
    return mu * (u**2 - v)

# --- Initial condition: nearly uniform u=v=1 with +-0.1 noise ---
np.random.seed(10)
u = 1.0 + 0.2 * (np.random.rand(N) - 0.5)   # uniform noise in [-0.1, 0.1]
v = 1.0 + 0.2 * (np.random.rand(N) - 0.5)

# --- Laplacian with no-flux (Neumann) boundary conditions ---
def laplacian(y):
    lap = np.empty_like(y)
    lap[1:-1] = (y[2:] - 2.0 * y[1:-1] + y[:-2]) / dX**2   # interior
    lap[0]    = (2.0 * y[1] - 2.0 * y[0]) / dX**2          # left  no-flux
    lap[-1]   = (2.0 * y[-2] - 2.0 * y[-1]) / dX**2        # right no-flux
    return lap

# --- One explicit forward-Euler reaction-diffusion step ---
def step(u, v):
    u_new = u + dt * (Du * laplacian(u) + f(u, v))   # activator update
    v_new = v + dt * (Dv * laplacian(v) + g(u, v))   # inhibitor update
    return u_new, v_new

# stability check for the fast (inhibitor) component
r = Dv * dt / dX**2
print(f"Diffusion number r = Dv*dt/dX^2 = {r:.4f} (must be < 0.5 for stability)")

# --- Run in successive time blocks until the pattern stops changing ---
steps_per_block = 20000        # 200 time units per block
n_blocks        = 15
u_prev_block    = u.copy()
for b in range(n_blocks):
    for _ in range(steps_per_block):
        u, v = step(u, v)      # integrate one dt forward
    # measure how much u changed across this block (stationarity monitor)
    change = np.max(np.abs(u - u_prev_block))
    u_prev_block = u.copy()
    print(f"Block {b+1:2d}: total time = {(b+1)*steps_per_block*dt:7.1f}, "
          f"max |du| over block = {change:.3e}")

total_time = n_blocks * steps_per_block * dt
print(f"Total integration time = {total_time:.1f}")

# --- Report the stationary pattern ---
print(f"u range: min = {u.min():.4f}, max = {u.max():.4f}")
print(f"v range: min = {v.min():.4f}, max = {v.max():.4f}")

# ---------------------------------------------------------------
# CHECK: stationary periodic Turing pattern with u and v IN PHASE
# ---------------------------------------------------------------
# 1) Stationarity: change over the final block should be tiny.
final_change = np.max(np.abs(u - u_prev_block))  # u_prev_block == u now -> 0
# recompute change of the last block explicitly from stored value above:
print(f"Stationarity: max |du| over final block = {change:.3e} (near zero => stationary)")

# 2) Periodicity: count interior peaks of u (multiple peaks => periodic).
def count_peaks(y):
    return int(np.sum((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])))
n_peaks_u = count_peaks(u)
n_peaks_v = count_peaks(v)
print(f"Number of interior u peaks = {n_peaks_u} (>1 => spatially periodic)")
print(f"Number of interior v peaks = {n_peaks_v}")

# 3) Phase relation: correlation of the spatial profiles about their means.
uc = u - u.mean()
vc = v - v.mean()
corr = np.sum(uc * vc) / np.sqrt(np.sum(uc**2) * np.sum(vc**2))
print(f"Spatial correlation corr(u, v) = {corr:.4f}")
print("Positive correlation => u and v peak TOGETHER (in phase, activator-inhibitor),")
print("whereas the substrate-depletion case gives negative correlation (out of phase).")

# 4) Confirm peak locations coincide: location of global u max vs v max.
print(f"Location of global u max: X = {X[np.argmax(u)]:.2f}")
print(f"Location of global v max: X = {X[np.argmax(v)]:.2f}")

# Explanation of why this confirms the result:
print("Check rationale: a near-zero final-block change plus multiple peaks proves a "
      "STATIONARY PERIODIC pattern, and a POSITIVE u-v correlation proves the activator "
      "and inhibitor peak together, distinguishing it from out-of-phase substrate depletion.")

# --- Plot the stationary periodic pattern ---
plt.figure(figsize=(9, 5))
plt.plot(X, u, 'b-o', ms=3, label='u (activator)')
plt.plot(X, v, 'r-s', ms=3, label='v (inhibitor)')
plt.xlabel('X')
plt.ylabel('concentration')
plt.title('Gierer-Meinhardt Turing pattern: u and v peak together (in phase)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.3.1_s3.png")
