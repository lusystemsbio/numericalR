import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Gillespie SSA for a self-repressing gene circuit.
# State: promoter s in {0 = D unbound, 1 = C bound}; M = molecule count.
# Reactions:
#   R1 production from D : D -> D + M   rate g0        (only if s==0)
#   R2 production from C : C -> C + M   rate g1        (only if s==1)
#   R3 degradation       : M -> 0       rate k*M
#   R4 binding 2M+D->C   : s:0->1, M-=2  rate kon*M*(M-1)  (only if s==0, M>=2)
#   R5 unbinding C->D+2M : s:1->0, M+=2  rate koff             (only if s==1)
# ---------------------------------------------------------------

def gillespie(g0, g1, kon, koff, k, T, rng, M0=0, s0=0):
    t = 0.0
    M = M0
    s = s0
    ts = [t]
    Ms = [M]
    while t < T:
        # --- compute propensities for current state ---
        a_prod = g0 if s == 0 else g1          # production (D or C)
        a_deg  = k * M                          # degradation
        a_bind = kon * M * (M - 1) if s == 0 else 0.0   # dimer binds free promoter
        a_unb  = koff if s == 1 else 0.0        # promoter releases dimer
        a = np.array([a_prod, a_deg, a_bind, a_unb])
        a0 = a.sum()
        if a0 <= 0:
            break
        # --- draw waiting time (exponential) ---
        tau = rng.exponential(1.0 / a0)
        t += tau
        # --- choose which reaction fires (proportional to propensity) ---
        r = rng.random() * a0
        cum = np.cumsum(a)
        j = int(np.searchsorted(cum, r))
        if j == 0:        # production
            M += 1
        elif j == 1:      # degradation
            M -= 1
        elif j == 2:      # binding: consume 2 M, promoter -> bound
            s = 1
            M -= 2
        else:             # unbinding: release 2 M, promoter -> free
            s = 0
            M += 2
        ts.append(t)
        Ms.append(M)
    return np.array(ts), np.array(Ms)

def time_weighted_stats(ts, Ms, burnin):
    # piecewise-constant M between events -> time-weighted mean/std
    mask = ts[:-1] >= burnin
    dt = np.diff(ts)[mask]
    vals = Ms[:-1][mask]
    W = dt.sum()
    mean = np.sum(vals * dt) / W
    var = np.sum((vals - mean) ** 2 * dt) / W
    return mean, np.sqrt(var)

# --- parameters ---
k = 0.1
T = 2000.0
burnin = 200.0
rng = np.random.default_rng(101)

regimes = {
    "self-inhibition": dict(g0=55, g1=5,  kon=0.002, koff=90),
    "no feedback":     dict(g0=30, g1=30, kon=0.0,   koff=90),
    "low copy number": dict(g0=5.5, g1=0.5, kon=0.2, koff=90),
}

results = {}
fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
for ax, (name, p) in zip(axes, regimes.items()):
    ts, Ms = gillespie(p["g0"], p["g1"], p["kon"], p["koff"], k, T, rng)
    mean, std = time_weighted_stats(ts, Ms, burnin)
    results[name] = (mean, std)
    ax.step(ts, Ms, where="post", lw=0.6)
    ax.axhline(mean, color="C1", lw=1.2, label=f"mean={mean:.2f}")
    ax.axhline(mean + std, color="C3", ls="--", lw=0.9, label=f"std={std:.2f}")
    ax.axhline(mean - std, color="C3", ls="--", lw=0.9)
    ax.set_ylabel("M count")
    ax.set_title(name)
    ax.legend(loc="upper right", fontsize=8)
axes[-1].set_xlabel("time")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.4.1_s5.png")

# --- report ---
for name, (mean, std) in results.items():
    print(f"{name}: mean = {mean:.4f}, std = {std:.4f}, CV = {std/mean:.4f}")

std_self = results["self-inhibition"][1]
std_none = results["no feedback"][1]
cv_low   = results["low copy number"][1] / results["low copy number"][0]
cv_none  = results["no feedback"][1] / results["no feedback"][0]

print(f"std(self-inhibition) < std(no feedback): {std_self < std_none} "
      f"({std_self:.4f} < {std_none:.4f})")
print(f"CV(low copy number) = {cv_low:.4f} vs CV(no feedback) = {cv_none:.4f}, "
      f"burstier: {cv_low > cv_none}")

# One-sentence explanation:
print("Explanation: a lower standard deviation in the self-inhibiting gene at a "
      "comparable mean shows the repression feedback damps fluctuations, while the "
      "much larger coefficient of variation at low copy number confirms that "
      "few-molecule binding produces strongly bursty dynamics.")
