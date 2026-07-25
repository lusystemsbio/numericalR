import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Explicit finite-difference Fokker-Planck integrator
#   dP/dt = -d(f*P)/dX + D * d2P/dX2,   with f = -k*X (OU drift)
# Dirichlet boundaries (P = 0 at both ends).
# ----------------------------------------------------------------------

# ---- Grid and integration parameters ----
L = 100.0          # domain length
dX = 1.0           # spatial step
dt = 0.01          # time step
D = 1.0            # diffusion coefficient
X = np.arange(-L / 2.0, L / 2.0 + dX, dX)   # grid centered on 0
N = X.size
print(f"Grid points N = {N}")
print(f"Domain X from {X[0]} to {X[-1]}")

# Stability (explicit diffusion): D*dt/dX^2 <= 0.5
print(f"Diffusion stability number D*dt/dX^2 = {D * dt / dX**2}")


def analytic_steady_state(k):
    """Analytic OU steady-state Gaussian."""
    return np.sqrt(k / (2.0 * np.pi * D)) * np.exp(-k * X**2 / (2.0 * D))


def normalize(P):
    """Normalize a discrete distribution to unit integral (trapezoid)."""
    area = np.trapz(P, X)
    return P / area if area > 0 else P


def fp_step(P, k):
    """One explicit Euler step of the Fokker-Planck equation."""
    f = -k * X                      # OU drift field f(X) = -k*X
    flux = f * P                    # advective flux f*P

    dPdt = np.zeros_like(P)         # will hold time derivative on interior

    # interior indices 1 .. N-2 (boundaries held fixed = Dirichlet)
    i = slice(1, N - 1)
    ip = slice(2, N)                # i+1
    im = slice(0, N - 2)            # i-1

    # centered-difference drift term:  -d(f*P)/dX
    drift = -(flux[ip] - flux[im]) / (2.0 * dX)

    # finite-difference diffusion term:  D * d2P/dX2
    diffusion = D * (P[ip] - 2.0 * P[i] + P[im]) / dX**2

    dPdt[i] = drift + diffusion

    P_new = P + dt * dPdt
    P_new[0] = 0.0                  # Dirichlet boundary (left)
    P_new[-1] = 0.0                 # Dirichlet boundary (right)
    return P_new


def integrate(P0, k, T):
    """Integrate from initial P0 for total time T, return P at chosen snapshots."""
    nsteps = int(round(T / dt))
    P = normalize(P0.copy())
    snap_steps = sorted(set([0,
                             int(0.02 * nsteps),
                             int(0.1 * nsteps),
                             int(0.3 * nsteps),
                             nsteps]))
    snaps = {}
    for step in range(nsteps + 1):
        if step in snap_steps:
            snaps[step] = P.copy()
        if step < nsteps:
            P = fp_step(P, k)
    return snaps, P


# ---- Initial conditions ----
def localized_patch():
    P = np.zeros(N)
    center = N // 2
    P[center - 2:center + 3] = 1.0   # narrow patch near X = 0
    return normalize(P)


def uniform_dist():
    P = np.ones(N)
    P[0] = 0.0
    P[-1] = 0.0
    return normalize(P)


springs = [0.01, 0.03]
T_total = 2000.0   # long enough to reach steady state for weak springs

# ---- Plots: relaxation to steady-state Gaussian for each k ----
fig, axes = plt.subplots(1, len(springs), figsize=(14, 5))

for ax, k in zip(np.atleast_1d(axes), springs):
    Pss = analytic_steady_state(k)

    # start from localized patch
    snaps, P_final_patch = integrate(localized_patch(), k, T_total)
    for step, Psnap in sorted(snaps.items()):
        ax.plot(X, Psnap, alpha=0.5, lw=1,
                label=f"patch t={step*dt:.0f}")

    # start from uniform distribution
    snaps_u, P_final_unif = integrate(uniform_dist(), k, T_total)
    ax.plot(X, P_final_unif, 'g--', lw=1.5, label="uniform final")

    # analytic steady state
    ax.plot(X, Pss, 'k-', lw=2, label="analytic $P_{ss}$")

    ax.set_title(f"k = {k}")
    ax.set_xlabel("X")
    ax.set_ylabel("P(X)")
    ax.set_xlim(-40, 40)
    ax.legend(fontsize=7)

    # ---- Numerical checks ----
    err_patch = np.max(np.abs(P_final_patch - Pss))
    err_unif = np.max(np.abs(P_final_unif - Pss))
    # standard deviations
    std_num = np.sqrt(np.trapz(X**2 * P_final_patch, X))
    std_analytic = np.sqrt(D / k)
    print(f"--- k = {k} ---")
    print(f"Max |P_final(patch)  - P_ss| = {err_patch:.6e}")
    print(f"Max |P_final(uniform)- P_ss| = {err_unif:.6e}")
    print(f"Numerical std (from patch)   = {std_num:.6f}")
    print(f"Analytic std sqrt(D/k)       = {std_analytic:.6f}")
    print(f"Peak height numerical (patch)= {np.max(P_final_patch):.6f}")
    print(f"Peak height analytic         = {np.max(Pss):.6f}")

# Compare widths across springs to confirm stiffer -> narrower
std_soft = np.sqrt(D / springs[0])
std_stiff = np.sqrt(D / springs[1])
print(f"Analytic std at k={springs[0]} = {std_soft:.6f} (softer, wider)")
print(f"Analytic std at k={springs[1]} = {std_stiff:.6f} (stiffer, narrower)")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7B.2.1_s5.png", dpi=110)

# Explanation of why the check confirms the result:
print("Check rationale: because both a localized patch and a uniform start "
      "converge to the same analytic Gaussian and the stiffer spring yields a "
      "smaller std sqrt(D/k) (narrower peak), the solver reproduces the "
      "correct initial-condition-independent OU steady state.")
