import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
L = 100          # domain length
dX = 1.0         # spatial step
dt = 0.01        # time step
D = 1.0          # diffusion coefficient

# --- Spatial grid centered on X=0 (from -L/2 to +L/2) ---
X = np.arange(-L / 2, L / 2 + dX, dX)   # grid points
N = X.size

# --- Stability check for the explicit scheme (FTCS) ---
alpha = D * dt / dX**2                   # diffusion number
print(f"Diffusion number alpha = D*dt/dX^2 = {alpha:.6f}  (stable if <= 0.5)")

# --- Initial condition: all probability concentrated at X = 0 ---
P = np.zeros(N)
i0 = np.argmin(np.abs(X))                # index of X=0
P[i0] = 1.0 / dX                         # unit probability mass (density = 1/dX)

# --- Time stepping setup ---
t_end = 100.0
n_steps = int(round(t_end / dt))
record_times = [1.0, 5.0, 20.0, 50.0, 100.0]   # snapshot times for P(X)
record_steps = {int(round(rt / dt)): rt for rt in record_times}

# storage for diagnostics vs time
t_hist = []
var_hist = []
mean_hist = []
mass_hist = []
snapshots = {}

def moments(P, X, dX):
    """Return total probability, mean, and variance of the distribution P(X)."""
    mass = np.sum(P) * dX                          # total probability
    mean = np.sum(X * P) * dX / mass               # normalized mean
    var = np.sum((X - mean)**2 * P) * dX / mass    # normalized variance
    return mass, mean, var

# record initial state
m0, mu0, v0 = moments(P, X, dX)
t_hist.append(0.0); mass_hist.append(m0); mean_hist.append(mu0); var_hist.append(v0)

# --- Explicit finite-difference integration (FTCS), done step-by-step ---
for step in range(1, n_steps + 1):
    # second spatial derivative via central differences on interior points
    d2P = np.zeros(N)
    d2P[1:-1] = (P[2:] - 2.0 * P[1:-1] + P[:-2]) / dX**2
    # forward Euler update in time: P_new = P + D*dt*d2P
    P = P + D * dt * d2P
    # Dirichlet (absorbing) boundaries: ends held at zero, so mass leaks out
    P[0] = 0.0
    P[-1] = 0.0

    # record diagnostics every so often (and always at record times)
    if step % 100 == 0 or step in record_steps:
        t = step * dt
        mass, mean, var = moments(P, X, dX)
        t_hist.append(t); mass_hist.append(mass)
        mean_hist.append(mean); var_hist.append(var)
    if step in record_steps:
        snapshots[record_steps[step]] = P.copy()

t_hist = np.array(t_hist)
var_hist = np.array(var_hist)
mean_hist = np.array(mean_hist)
mass_hist = np.array(mass_hist)

# --- Print numerical results at the snapshot times ---
for rt in record_times:
    mass, mean, var = moments(snapshots[rt], X, dX)
    print(f"t = {rt:6.1f} | total_prob = {mass:.6f} | mean = {mean:+.4f} | "
          f"variance = {var:.4f} | theory 2Dt = {2*D*rt:.4f}")

# early-time comparison (where boundaries not yet reached)
print("\n--- Early-time variance vs theory (2*D*t) ---")
for t, v in zip(t_hist, var_hist):
    if t in record_times[:3]:
        print(f"t = {t:6.1f} | measured var = {v:.4f} | 2*D*t = {2*D*t:.4f} | "
              f"ratio = {v/(2*D*t):.4f}")

print(f"\nFinal total probability at t={t_end:.0f}: {mass_hist[-1]:.6f} "
      f"(started at 1.0; deficit leaked through absorbing ends)")

# --- Plot: P(X) spreading over time, and variance vs t ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for rt in record_times:
    ax1.plot(X, snapshots[rt], label=f"t = {rt:g}")
ax1.set_xlabel("X")
ax1.set_ylabel("P(X)")
ax1.set_title("Diffusion of a point distribution")
ax1.set_xlim(-40, 40)
ax1.legend()

ax2.plot(t_hist, var_hist, "b-", lw=2, label="measured variance")
ax2.plot(t_hist, 2 * D * t_hist, "r--", lw=2, label="theory 2*D*t")
ax2.set_xlabel("t")
ax2.set_ylabel("variance")
ax2.set_title("Variance growth vs theory")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.2.1_s2.png")

# --- One-sentence explanation of why the check confirms the result ---
print("\nExplanation: A pure diffusion process must turn a delta spike into a "
      "Gaussian whose variance rises linearly as 2*D*t, so seeing the simulated "
      "spread match 2*D*t at early times (before mass reaches the absorbing ends "
      "and slowly leaks) confirms the integrator reproduces the correct physics.")
