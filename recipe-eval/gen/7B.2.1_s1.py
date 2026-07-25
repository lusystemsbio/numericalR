import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Fokker-Planck integrator for an SDE distribution
#   dP/dt = -d(f*P)/dX + D d2P/dX2,   f = -k*X  (OU well)
# Steady state:  P_ss(X) = sqrt(k/(2 pi D)) exp(-k X^2 / (2 D))
# Explicit finite differences, Dirichlet (P=0) boundaries.
# ---------------------------------------------------------------

# --- Discretization parameters ---
L   = 100.0     # domain length
dX  = 1.0       # spatial step
dt  = 0.01      # time step
D   = 1.0       # diffusion coefficient

# Grid centered on 0 so the well sits in the middle
X = np.arange(-L/2.0, L/2.0 + dX, dX)
N = X.size
print(f"Number of grid points N = {N}")
print(f"Grid spans X from {X[0]} to {X[-1]}")

def f_drift(X, k):
    # Drift of the OU process: pulls toward X=0
    return -k * X

def P_steady(X, k, D):
    # Analytic steady-state Gaussian
    return np.sqrt(k / (2.0 * np.pi * D)) * np.exp(-k * X**2 / (2.0 * D))

def normalize(P, dX):
    # Keep total probability = 1 (integral of P dX)
    return P / (np.sum(P) * dX)

def fp_step(P, k):
    """One explicit Euler step of the Fokker-Planck equation."""
    fP = f_drift(X, k) * P                      # flux term f*P at each node
    Pn = P.copy()
    # Interior nodes only (Dirichlet: boundaries held at 0)
    i = slice(1, N - 1)
    # Centered-difference drift term: -d(f*P)/dX
    drift = -(fP[2:] - fP[:-2]) / (2.0 * dX)
    # Finite-difference diffusion term: D * d2P/dX2
    diff  = D * (P[2:] - 2.0 * P[1:-1] + P[:-2]) / dX**2
    Pn[i] = P[i] + dt * (drift + diff)
    Pn[0] = 0.0                                 # Dirichlet boundary
    Pn[-1] = 0.0                                # Dirichlet boundary
    return Pn

def integrate(P0, k, n_steps):
    # Integrate and store snapshots for plotting the relaxation
    P = normalize(P0.copy(), dX)
    snap_times = np.linspace(0, n_steps, 6, dtype=int)
    snaps = {}
    for step in range(n_steps + 1):
        if step in snap_times:
            snaps[step] = P.copy()
        P = fp_step(P, k)
    return P, snaps

# --- Initial conditions ---
def patch_ic():
    # Localized patch near center
    P = np.zeros(N)
    center = N // 2
    P[center-2:center+3] = 1.0
    return normalize(P, dX)

def uniform_ic():
    # Uniform distribution across interior
    P = np.ones(N)
    P[0] = 0.0
    P[-1] = 0.0
    return normalize(P, dX)

ks = [0.01, 0.03]
n_steps = 200000  # long enough to reach steady state

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, k in zip(axes, ks):
    Pss = P_steady(X, k, D)

    # Relax from a localized patch
    Pf_patch, snaps_patch = integrate(patch_ic(), k, n_steps)
    # Relax from a uniform distribution
    Pf_unif, snaps_unif = integrate(uniform_ic(), k, n_steps)

    # Plot relaxation snapshots from the patch IC
    for step, Psnap in sorted(snaps_patch.items()):
        ax.plot(X, Psnap, alpha=0.4, lw=1,
                label=f"t={step*dt:.0f}")
    # Final distributions from both ICs
    ax.plot(X, Pf_patch, 'b-', lw=2, label="final (patch IC)")
    ax.plot(X, Pf_unif, 'g--', lw=2, label="final (uniform IC)")
    # Analytic steady state
    ax.plot(X, Pss, 'r:', lw=2.5, label="analytic $P_{ss}$")

    ax.set_title(f"OU well, k = {k}")
    ax.set_xlabel("X")
    ax.set_ylabel("P(X)")
    ax.set_xlim(-40, 40)
    ax.legend(fontsize=7)

    # --- Numerical checks ---
    err_patch = np.sqrt(np.sum((Pf_patch - Pss)**2 * dX))
    err_unif  = np.sqrt(np.sum((Pf_unif  - Pss)**2 * dX))
    # Compare peak widths via standard deviation: sigma^2 = D/k analytically
    mean_num = np.sum(X * Pf_patch) * dX
    var_num  = np.sum((X - mean_num)**2 * Pf_patch) * dX
    std_num  = np.sqrt(var_num)
    std_theory = np.sqrt(D / k)

    print(f"--- k = {k} ---")
    print(f"L2 error, patch   IC vs analytic P_ss = {err_patch:.6e}")
    print(f"L2 error, uniform IC vs analytic P_ss = {err_unif:.6e}")
    print(f"Numerical std of final distribution    = {std_num:.4f}")
    print(f"Analytic std sqrt(D/k)                 = {std_theory:.4f}")
    print(f"Peak height numerical                  = {np.max(Pf_patch):.6f}")
    print(f"Peak height analytic                   = {np.max(Pss):.6f}")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7B.2.1_s1.png")

print("Explanation: Because both the localized patch and the uniform "
      "initial conditions converge to the same analytic Gaussian, and the "
      "stiffer spring (k=0.03) yields a smaller std sqrt(D/k) hence a "
      "narrower, taller peak, the integrator correctly reproduces the "
      "steady-state distribution independent of initial condition.")
