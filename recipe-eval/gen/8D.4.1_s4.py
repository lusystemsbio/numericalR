import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Gillespie SSA for a self-repressing gene.
# Promoter is a single copy in one of two states:
#   D = 1 (unbound, produces mRNA/protein M at fast rate g0)
#   C = 1 (dimer-bound, produces at slower rate g1)   with D + C = 1.
# Reactions:
#   R1: production from D    D -> D + M        rate  g0 * D
#   R2: production from C    C -> C + M        rate  g1 * C
#   R3: degradation          M -> 0            rate  k  * M
#   R4: binding   2M + D -> C                  rate  kon * M*(M-1) * D
#   R5: unbinding C -> 2M + D                  rate  koff * C
# ----------------------------------------------------------------------


def gillespie(g0, g1, k, kon, koff, T, rng, M0=0):
    """Explicit Gillespie SSA; returns arrays of times and molecule counts M(t)."""
    t = 0.0
    M = M0          # molecule count
    D = 1           # promoter starts unbound
    C = 0           # promoter starts not dimer-bound
    ts = [t]
    Ms = [M]
    while t < T:
        # --- compute the five reaction propensities from the current state ---
        a1 = g0 * D                 # production while unbound
        a2 = g1 * C                 # production while bound
        a3 = k * M                  # degradation of one molecule
        a4 = kon * M * (M - 1) * D  # two molecules bind the free promoter (dimer)
        a5 = koff * C               # bound complex releases the dimer
        a0 = a1 + a2 + a3 + a4 + a5
        if a0 <= 0.0:               # no reaction possible -> jump to end
            break
        # --- draw the waiting time from an exponential with rate a0 ---
        tau = rng.exponential(1.0 / a0)
        t += tau
        if t > T:                   # do not record past the horizon
            break
        # --- pick which reaction fires, proportional to its propensity ---
        r = rng.random() * a0
        if r < a1:                  # R1
            M += 1
        elif r < a1 + a2:           # R2
            M += 1
        elif r < a1 + a2 + a3:      # R3
            M -= 1
        elif r < a1 + a2 + a3 + a4: # R4: binding consumes 2 M, flips promoter
            M -= 2
            D, C = 0, 1
        else:                       # R5: unbinding releases 2 M, flips promoter
            M += 2
            D, C = 1, 0
        ts.append(t)
        Ms.append(M)
    return np.array(ts), np.array(Ms)


def sampled_stats(ts, Ms, T, burn_frac=0.2, n=2000):
    """Sample the piecewise-constant trajectory on a uniform grid (after burn-in)
    and return the grid, sampled values, mean and standard deviation."""
    grid = np.linspace(burn_frac * T, T, n)
    idx = np.searchsorted(ts, grid, side="right") - 1   # last event at or before each grid time
    idx = np.clip(idx, 0, len(Ms) - 1)
    vals = Ms[idx]
    return grid, vals, vals.mean(), vals.std()


# ----------------------------------------------------------------------
# Three regimes; k = 0.1 throughout, seed 101.
# ----------------------------------------------------------------------
k = 0.1
T = 4000.0
rng = np.random.default_rng(101)

regimes = {
    "self-inhibition": dict(g0=55.0, g1=5.0,  kon=0.002, koff=90.0),
    "no feedback":     dict(g0=30.0, g1=30.0, kon=0.002, koff=90.0),
    "low copy number": dict(g0=5.5,  g1=0.5,  kon=0.2,   koff=90.0),
}

results = {}
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

for ax, (name, p) in zip(axes, regimes.items()):
    ts, Ms = gillespie(p["g0"], p["g1"], k, p["kon"], p["koff"], T, rng)
    grid, vals, mean, std = sampled_stats(ts, Ms, T)
    results[name] = dict(mean=mean, std=std)

    ax.step(ts, Ms, where="post", lw=0.6, color="steelblue")
    ax.axhline(mean, color="black", ls="-", lw=1.2, label=f"mean = {mean:.2f}")
    ax.axhline(mean + std, color="red", ls="--", lw=1.0, label=f"mean +/- sd (sd = {std:.2f})")
    ax.axhline(mean - std, color="red", ls="--", lw=1.0)
    ax.set_ylabel("molecule count M")
    ax.set_title(f"{name}: g0={p['g0']}, g1={p['g1']}, kon={p['kon']}, koff={p['koff']}")
    ax.legend(loc="upper right", fontsize=8)

axes[-1].set_xlabel("time")
fig.suptitle("Gillespie SSA: negative autoregulation suppresses gene-expression noise")
fig.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.4.1_s4.png", dpi=130)

# ----------------------------------------------------------------------
# Report numerical results.
# ----------------------------------------------------------------------
for name in regimes:
    r = results[name]
    cv = r["std"] / r["mean"] if r["mean"] != 0 else float("nan")
    print(f"{name}: mean = {r['mean']:.4f}")
    print(f"{name}: std  = {r['std']:.4f}")
    print(f"{name}: coefficient of variation (std/mean) = {cv:.4f}")

# ----------------------------------------------------------------------
# Separate checks.
# ----------------------------------------------------------------------
sd_self = results["self-inhibition"]["std"]
sd_none = results["no feedback"]["std"]
check_noise = sd_self < sd_none
print(f"\nCheck 1 - self-inhibition std ({sd_self:.4f}) < no-feedback std ({sd_none:.4f}): {check_noise}")

cv_low = results["low copy number"]["std"] / results["low copy number"]["mean"]
cv_none = results["no feedback"]["std"] / results["no feedback"]["mean"]
check_bursty = cv_low > cv_none
print(f"Check 2 - low-copy CV ({cv_low:.4f}) > no-feedback CV ({cv_none:.4f}), i.e. strongly bursty: {check_bursty}")

print("\nExplanation: A lower standard deviation for the self-inhibiting gene at a comparable mean "
      "means fluctuations around the average are smaller, which is exactly what noise suppression by "
      "negative autoregulation predicts (while the low-copy regime shows the opposite, an inflated "
      "coefficient of variation, because binding/unbinding of just a few molecules produces large bursts).")
