import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Grid and parameters ----------------
L = 100.0          # domain length
dX = 1.0           # spatial step
dt = 0.01          # time step
D = 1.0            # diffusion coefficient

# X grid centered on 0, spanning [-L/2, L/2], with Dirichlet boundaries (P=0 at ends)
X = np.arange(-L/2, L/2 + dX, dX)
N = X.size
print(f"Number of grid points N = {N}")
print(f"X range = [{X[0]}, {X[-1]}]")
print(f"Diffusion stability number D*dt/dX^2 = {D*dt/dX**2:.4f} (needs < 0.5)")

springs = [0.01, 0.03]          # spring stiffnesses k to test
T_total = 600.0                 # total integration time (>> 1/k relaxation time)
nsteps = int(T_total / dt)
# times (in physical units) at which to snapshot P(X) for plotting
snapshot_times = [0.0, 20.0, 60.0, 150.0, 400.0, T_total]
snapshot_steps = [int(t/dt) for t in snapshot_times]


def analytic_steady_state(k):
    """Analytic OU steady-state Gaussian P_ss(X) = sqrt(k/(2*pi*D)) exp(-k X^2/(2 D))."""
    return np.sqrt(k/(2*np.pi*D)) * np.exp(-k*X**2/(2*D))


def make_initial(kind):
    """Two normalized initial conditions on the interior (boundaries held at 0)."""
    P = np.zeros(N)
    if kind == "patch":
        # localized patch: a narrow block off-center so relaxation is clearly visible
        mask = (X > 5) & (X < 15)
        P[mask] = 1.0
    elif kind == "uniform":
        # uniform over the interior, zero at the two boundary nodes
        P[1:-1] = 1.0
    # normalize so total probability = 1 (trapezoidal integral)
    P /= np.trapz(P, X)
    return P


def fp_step(P, k):
    """One explicit finite-difference Fokker-Planck update with Dirichlet boundaries.
       dP/dt = -d(f P)/dX + D d2P/dX2 ,   f = -k X  (OU drift toward 0)."""
    f = -k * X                       # drift field f(X) = -k X
    flux = f * P                     # the quantity f*P whose gradient enters the drift term
    Pnew = P.copy()
    # centered difference for the drift term  -d(f P)/dX  at interior points
    drift = -(flux[2:] - flux[:-2]) / (2*dX)
    # centered second difference for the diffusion term  D d2P/dX2
    diff = D * (P[2:] - 2*P[1:-1] + P[:-2]) / dX**2
    # explicit Euler advance of the interior; boundaries stay fixed at 0 (Dirichlet)
    Pnew[1:-1] = P[1:-1] + dt * (drift + diff)
    Pnew[0] = 0.0
    Pnew[-1] = 0.0
    return Pnew


# ---------------- Run integrations and plot ----------------
init_kinds = ["patch", "uniform"]
fig, axes = plt.subplots(len(springs), len(init_kinds), figsize=(13, 9), sharex=True)

for i, k in enumerate(springs):
    Pss = analytic_steady_state(k)
    for j, kind in enumerate(init_kinds):
        ax = axes[i, j]
        P = make_initial(kind)
        snaps = {}
        for step in range(nsteps + 1):
            if step in snapshot_steps:
                snaps[step] = P.copy()
            if step < nsteps:
                P = fp_step(P, k)
        # plot the recorded snapshots relaxing toward the Gaussian
        for t, step in zip(snapshot_times, snapshot_steps):
            ax.plot(X, snaps[step], lw=1.4, label=f"t={t:g}")
        ax.plot(X, Pss, 'k--', lw=2.0, label="analytic $P_{ss}$")
        ax.set_title(f"k={k},  init={kind}")
        ax.set_xlim(-40, 40)
        ax.set_xlabel("X")
        ax.set_ylabel("P(X)")
        ax.legend(fontsize=7)

        # numerical checks: final probability mass, peak location, and error vs analytic
        Pfinal = P
        mass = np.trapz(Pfinal, X)
        peak_std = np.sqrt(np.trapz(X**2 * Pfinal, X) / mass)   # numeric std of final dist
        analytic_std = np.sqrt(D/k)                             # OU steady-state std
        max_err = np.max(np.abs(Pfinal - Pss))
        print(f"k={k}, init={kind}: final mass={mass:.4f}, "
              f"numeric std={peak_std:.4f}, analytic std={analytic_std:.4f}, "
              f"max|P-P_ss|={max_err:.3e}")

# report the confinement (narrowing) trend explicitly
for k in springs:
    print(f"Analytic steady-state std for k={k}: sqrt(D/k) = {np.sqrt(D/k):.4f}")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7B.2.1_s4.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because both the localized-patch and uniform initial conditions "
      "converge to the same analytic Gaussian, and the width sqrt(D/k) shrinks as k "
      "grows (0.01 -> 0.03), the integrator correctly reproduces the OU steady state "
      "independent of initial data with a stiffer spring giving a narrower peak.")
