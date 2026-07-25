import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor Turing model (multi-component RD)
#   u_t = Du*u_xx + f(u,v),   f = u^2/v - u          (activator, slow diffusion)
#   v_t = Dv*v_xx + g(u,v),   g = mu*(u^2 - v)       (inhibitor, fast diffusion)
# Local self-enhancement (u^2/v) + long-range inhibition (Dv >> Du) => Turing.
# ----------------------------------------------------------------------

# ---- parameters ----
L   = 20.0     # domain length
dX  = 0.2      # spatial step
dt  = 0.01     # time step
d   = 0.1      # activator diffusivity Du (slow)
Du  = d
Dv  = 1.0      # inhibitor diffusivity (fast)
mu  = 1.5      # kinetic parameter

nX  = int(round(L / dX)) + 1        # number of grid points
X   = np.linspace(0.0, L, nX)

# ---- reaction terms ----
def f(u, v):   # activator kinetics
    return u**2 / v - u

def g(u, v):   # inhibitor kinetics
    return mu * (u**2 - v)

# ---- Laplacian with no-flux (Neumann) boundaries, computed explicitly ----
def laplacian(y):
    lap = np.empty_like(y)
    # interior: standard second difference
    lap[1:-1] = (y[2:] - 2.0*y[1:-1] + y[:-2]) / dX**2
    # zero-flux boundaries: ghost point equals neighbor
    lap[0]  = (y[1]  - y[0])  * 2.0 / dX**2
    lap[-1] = (y[-2] - y[-1]) * 2.0 / dX**2
    return lap

# ---- initial condition: nearly uniform u=v=1 with +-0.1 noise ----
np.random.seed(10)
u = 1.0 + 0.1 * (2.0*np.random.rand(nX) - 1.0)
v = 1.0 + 0.1 * (2.0*np.random.rand(nX) - 1.0)

# ---- explicit Euler integration, run in successive time blocks ----
block_steps = 5000          # steps per block
n_blocks    = 40            # total time = block_steps*n_blocks*dt = 2000
prev_u = u.copy()
for b in range(n_blocks):
    for _ in range(block_steps):
        # explicit forward-Euler update of both components at once
        u_new = u + dt * (Du * laplacian(u) + f(u, v))
        v_new = v + dt * (Dv * laplacian(v) + g(u, v))
        u, v = u_new, v_new
    # measure how much the pattern still changes between blocks (approach to steady state)
    change = np.max(np.abs(u - prev_u))
    prev_u = u.copy()
    print(f"Block {b+1:2d}  t={((b+1)*block_steps*dt):7.1f}  max|du| per block = {change:.3e}")

# ----------------------------------------------------------------------
# Check that the pattern is a stationary periodic Turing pattern with
# u and v PEAKING TOGETHER (in phase), unlike substrate depletion.
# ----------------------------------------------------------------------

# stationarity: residual of the RHS (should be ~0 if stationary)
res_u = np.max(np.abs(Du * laplacian(u) + f(u, v)))
res_v = np.max(np.abs(Dv * laplacian(v) + g(u, v)))
print(f"Stationarity residual max|du/dt| = {res_u:.3e}")
print(f"Stationarity residual max|dv/dt| = {res_v:.3e}")

# periodicity: dominant spatial wavelength from the FFT of u (drop k=0 mean)
u_fluct = u - u.mean()
spec = np.abs(np.fft.rfft(u_fluct))
freqs = np.fft.rfftfreq(nX, d=dX)
kdom = freqs[1:][np.argmax(spec[1:])]         # dominant spatial frequency
wavelength = 1.0 / kdom if kdom > 0 else np.inf
print(f"Number of interior peaks in u = {int(np.sum((u[1:-1] > u[:-2]) & (u[1:-1] > u[2:])))}")
print(f"Dominant spatial wavelength    = {wavelength:.3f}")

# in-phase check: spatial correlation between u and v fluctuations
v_fluct = v - v.mean()
corr = np.corrcoef(u_fluct, v_fluct)[0, 1]
print(f"Spatial correlation corr(u,v)  = {corr:+.3f}")
print(f"u and v in phase (peak together)? {'YES' if corr > 0 else 'NO'}")

# report peak locations to show they coincide
u_peaks = X[1:-1][(u[1:-1] > u[:-2]) & (u[1:-1] > u[2:])]
v_peaks = X[1:-1][(v[1:-1] > v[:-2]) & (v[1:-1] > v[2:])]
print(f"u peak locations X = {np.round(u_peaks, 2).tolist()}")
print(f"v peak locations X = {np.round(v_peaks, 2).tolist()}")

# One-sentence explanation of the check:
# A positive spatial correlation between u and v (both patterned, both stationary)
# confirms activator-inhibitor behavior: the activator and inhibitor peak at the
# same locations (in phase), whereas substrate depletion would give corr < 0
# (activator peaks where the depleted substrate is in a trough, i.e. out of phase).
print("Explanation: corr(u,v) > 0 means the activator u and inhibitor v peak at the "
      "same places (in phase), the activator-inhibitor signature, unlike the "
      "out-of-phase (corr < 0) substrate-depletion pattern.")

# ---- plot the stationary periodic Turing pattern ----
plt.figure(figsize=(9, 5))
plt.plot(X, u, 'b-o', ms=3, label='activator u(X)')
plt.plot(X, v, 'r-s', ms=3, label='inhibitor v(X)')
plt.xlabel('X')
plt.ylabel('concentration')
plt.title('Gierer-Meinhardt Turing pattern (u and v peak together)')
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.3.1_s1.png")
