import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Gillespie SSA for a self-repressing gene with promoter binding.
# Promoter is either unbound (D) or bound by a dimer (C).
# Reactions:
#   R1: D -> D + M          (production from unbound promoter, rate g0)
#   R2: C -> C + M          (production from bound  promoter, rate g1)
#   R3: M -> 0              (degradation, rate k per molecule)
#   R4: 2M + D -> C         (dimer binding,   propensity kon*M*(M-1)/2, only if D)
#   R5: C -> D + 2M         (dimer unbinding, rate koff,               only if C)
# ---------------------------------------------------------------

def gillespie(g0, g1, k, kon, koff, T, dt, rng):
    t = 0.0
    M = 0            # molecule count
    bound = False    # promoter state: False = D (unbound), True = C (bound)

    # time grid on which we record the trajectory (sample-and-hold)
    grid = np.arange(0.0, T + dt, dt)
    rec = np.empty_like(grid)
    gi = 0

    while t < T:
        # --- compute propensities for the current state ---
        a_prod = g0 if not bound else g1          # production depends on promoter state
        a_deg  = k * M                            # degradation scales with copy number
        a_bind = kon * M * (M - 1) / 2.0 if not bound else 0.0  # need a dimer + free D
        a_unb  = koff if bound else 0.0           # unbinding only if currently bound
        a_tot  = a_prod + a_deg + a_bind + a_unb

        if a_tot <= 0.0:
            break

        # --- draw waiting time to next reaction (exponential) ---
        tau = rng.exponential(1.0 / a_tot)
        t_next = t + tau

        # --- record state on every grid point crossed before the jump ---
        while gi < len(grid) and grid[gi] <= t_next:
            rec[gi] = M
            gi += 1

        t = t_next
        if t >= T:
            break

        # --- pick which reaction fires (proportional to propensity) ---
        r = rng.random() * a_tot
        if r < a_prod:
            M += 1                       # R1/R2: produce one molecule
        elif r < a_prod + a_deg:
            M -= 1                       # R3: degrade one molecule
        elif r < a_prod + a_deg + a_bind:
            bound = True; M -= 2         # R4: two monomers bind promoter as a dimer
        else:
            bound = False; M += 2        # R5: dimer releases, returning two monomers

    # fill any remaining grid points at the final state
    while gi < len(grid):
        rec[gi] = M
        gi += 1
    return grid, rec

# ---------------------------------------------------------------
# Parameters: k = 0.1 in all regimes; seed 101.
# ---------------------------------------------------------------
k = 0.1
T = 2000.0     # long run so mean/std are well sampled
dt = 0.5
burn = int(200.0 / dt)   # discard transient before computing statistics
rng = np.random.default_rng(101)

regimes = {
    "self-inhibition": dict(g0=55, g1=5,  kon=0.002, koff=90),
    "no feedback":     dict(g0=30, g1=30, kon=0.0,   koff=90),
    "low copy number": dict(g0=5.5, g1=0.5, kon=0.2, koff=90),
}

results = {}
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

for ax, (name, p) in zip(axes, regimes.items()):
    grid, M = gillespie(g0=p["g0"], g1=p["g1"], k=k,
                        kon=p["kon"], koff=p["koff"], T=T, dt=dt, rng=rng)
    m_ss = M[burn:]                     # steady-state portion
    mean = m_ss.mean()
    std  = m_ss.std()
    cv   = std / mean if mean > 0 else float("nan")
    results[name] = dict(mean=mean, std=std, cv=cv)

    print(f"[{name}] mean M = {mean:.3f}")
    print(f"[{name}] std  M = {std:.3f}")
    print(f"[{name}] CV (std/mean) = {cv:.3f}")

    ax.plot(grid, M, lw=0.6, color="steelblue")
    ax.axhline(mean, color="crimson", lw=1.2, label=f"mean={mean:.1f}")
    ax.fill_between(grid, mean - std, mean + std, color="crimson", alpha=0.15,
                    label=f"±std={std:.1f}")
    ax.set_title(f"{name}: g0={p['g0']}, g1={p['g1']}, kon={p['kon']}, koff={p['koff']}")
    ax.set_ylabel("molecule count M")
    ax.legend(loc="upper right", fontsize=8)

axes[-1].set_xlabel("time")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.4.1_s2.png", dpi=130)

# ---------------------------------------------------------------
# Checks
# ---------------------------------------------------------------
std_self = results["self-inhibition"]["std"]
std_none = results["no feedback"]["std"]
suppresses = std_self < std_none
print(f"CHECK noise suppression: std(self-inhibition)={std_self:.3f} < std(no feedback)={std_none:.3f} -> {suppresses}")

cv_low  = results["low copy number"]["cv"]
cv_none = results["no feedback"]["cv"]
bursty = cv_low > cv_none
print(f"CHECK low-copy burstiness: CV(low copy)={cv_low:.3f} > CV(no feedback)={cv_none:.3f} -> {bursty}")

# Explanation: a self-repressing gene with the same average output but a smaller
# standard deviation than the feedback-free gene is, by definition, expressing
# with tighter fluctuations around its mean, which is exactly what noise
# suppression means; the elevated CV at low copy number shows the large relative
# swings (bursts) expected when a handful of molecules control promoter binding.
print("EXPLANATION: A lower standard deviation at the same mean means tighter "
      "fluctuations about the average, so negative autoregulation directly reduces "
      "expression noise, while the larger CV at low copy number reflects the strong "
      "relative bursts caused by few-molecule promoter binding.")
