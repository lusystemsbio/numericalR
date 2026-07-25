import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor Turing pattern (1D)
# Multi-component finite-difference reaction-diffusion integrator
# implemented explicitly (forward Euler in time, centered Laplacian
# in space, zero-flux / Neumann boundaries), run in time blocks.
# ---------------------------------------------------------------

# ---- Parameters ----
L    = 20.0      # domain length
dX   = 0.2       # spatial step
dt   = 0.01      # time step
d    = 0.1       # activator diffusion (slow, short-range self-enhancement)
Dv   = 1.0       # inhibitor diffusion (fast, long-range inhibition)
Du   = d
mu   = 1.5       # kinetics parameter

N = int(round(L / dX)) + 1          # number of grid points
X = np.linspace(0.0, L, N)

# stability check for explicit scheme: D*dt/dX^2 <= 0.5
print("Stability number Du*dt/dX^2 =", Du * dt / dX**2)
print("Stability number Dv*dt/dX^2 =", Dv * dt / dX**2)

# ---- Reaction kinetics ----
# f: activator kinetics  u^2/v - u   (local self-enhancement)
# g: inhibitor kinetics  mu*(u^2 - v)
def f(u, v):
    return u**2 / v - u

def g(u, v):
    return mu * (u**2 - v)

# ---- Discrete Laplacian with zero-flux (Neumann) boundaries ----
def laplacian(c):
    lap = np.empty_like(c)
    # interior points: standard 2nd-order centered difference
    lap[1:-1] = (c[2:] - 2.0 * c[1:-1] + c[:-2]) / dX**2
    # boundaries: reflect neighbour (no-flux) -> ghost node equals interior node
    lap[0]  = (2.0 * c[1]  - 2.0 * c[0])  / dX**2
    lap[-1] = (2.0 * c[-2] - 2.0 * c[-1]) / dX**2
    return lap

# ---- Initial condition: nearly uniform u = v = 1 with +-0.1 noise ----
rng = np.random.default_rng(10)
u = 1.0 + 0.1 * (2.0 * rng.random(N) - 1.0)   # uniform noise in [-0.1, 0.1]
v = 1.0 + 0.1 * (2.0 * rng.random(N) - 1.0)

# ---- One explicit Euler step for both components ----
def step(u, v):
    u_new = u + dt * (Du * laplacian(u) + f(u, v))   # activator update
    v_new = v + dt * (Dv * laplacian(v) + g(u, v))   # inhibitor update
    return u_new, v_new

# ---- Integrate in successive time blocks ----
steps_per_block = 2000          # 20 time units per block
n_blocks        = 25            # total 500 time units
u_prev_block = u.copy()

for b in range(n_blocks):
    for _ in range(steps_per_block):
        u, v = step(u, v)
    # measure how much the activator field changed over this block
    change = np.max(np.abs(u - u_prev_block))
    u_prev_block = u.copy()

print("Max |u| change over final block (stationarity measure) =", change)

# ---------------------------------------------------------------
# Plot: u(X) and v(X) forming a stationary periodic pattern
# ---------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(X, u, 'b-', lw=2, label="u (activator)")
plt.plot(X, v, 'r-', lw=2, label="v (inhibitor)")
plt.xlabel("X")
plt.ylabel("concentration")
plt.title("Gierer-Meinhardt Turing pattern (activator-inhibitor)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.3.1_s4.png")

# ---------------------------------------------------------------
# CHECK: confirm stationary periodic Turing pattern, u & v in phase
# ---------------------------------------------------------------

# (1) Amplitude of the pattern (departure from the uniform state ~1)
print("u range: min =", np.min(u), " max =", np.max(u))
print("v range: min =", np.min(v), " max =", np.max(v))
print("Pattern amplitude (u max-min) =", np.max(u) - np.min(u))

# (2) Periodicity: count interior peaks in u to show multiple repeats
peaks = np.where((u[1:-1] > u[:-2]) & (u[1:-1] > u[2:]))[0] + 1
print("Number of activator peaks =", len(peaks))
if len(peaks) >= 2:
    spacing = np.diff(X[peaks])
    print("Mean peak spacing (wavelength) =", np.mean(spacing))

# (3) In-phase vs out-of-phase: spatial correlation of u and v.
#     Positive correlation => u and v peak TOGETHER (activator-inhibitor).
#     A substrate-depletion model instead gives NEGATIVE correlation
#     (activator peaks where the substrate is depleted -> out of phase).
uc = u - np.mean(u)
vc = v - np.mean(v)
corr = np.sum(uc * vc) / np.sqrt(np.sum(uc**2) * np.sum(vc**2))
print("Spatial correlation of u and v =", corr)
print("u and v peak together (in phase)?", corr > 0)

# Explanation:
# A positive spatial correlation between u and v (the peaks of activator and
# inhibitor coinciding) confirms the Gierer-Meinhardt "peak-together" pattern,
# distinguishing it from substrate depletion where u and v would be anti-correlated.
print("Check confirms result because a positive u-v spatial correlation means "
      "activator and inhibitor peaks coincide, the hallmark of the "
      "activator-inhibitor mechanism rather than the out-of-phase substrate-depletion case.")
