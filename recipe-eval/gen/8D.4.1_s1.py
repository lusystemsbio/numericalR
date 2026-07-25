import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Gillespie SSA for a self-repressing gene with promoter binding.
# Single promoter is either UNBOUND (D) or BOUND by a dimer (C): D + C = 1.
# Molecules M are produced from D at g0 and from C at the slower g1,
# degrade at rate k, and a dimer (2M) binds/unbinds the promoter.
# ---------------------------------------------------------------

def gillespie(g0, g1, k, kon, koff, T, seed):
    rng = np.random.default_rng(seed)   # dedicated RNG so each run is reproducible
    t = 0.0
    M = 0                               # start with zero molecules
    b = 0                               # promoter bound? 0 = unbound (D), 1 = bound (C)
    ts = [t]                            # record times ...
    Ms = [M]                            # ... and molecule counts
    while t < T:
        D = 1 - b                       # promoter is unbound only if not bound
        C = b
        # --- propensities of the five reaction channels ---
        a_prod0 = g0 * D                # production from unbound promoter
        a_prod1 = g1 * C                # production from bound (repressed) promoter
        a_deg   = k * M                 # first-order degradation M -> 0
        a_bind  = kon * M * (M - 1) * D # dimer forms and binds: 2M + D -> C
        a_unbind = koff * C             # dimer unbinds and releases: C -> D + 2M
        a = np.array([a_prod0, a_prod1, a_deg, a_bind, a_unbind])
        a0 = a.sum()
        if a0 <= 0:                     # no reaction possible; stop
            break
        # --- time to next reaction: exponential with rate a0 ---
        tau = rng.exponential(1.0 / a0)
        t += tau
        if t > T:
            break
        # --- choose which reaction fires, weighted by propensity ---
        r = rng.random() * a0
        cum = np.cumsum(a)
        idx = np.searchsorted(cum, r)
        # --- apply the chosen reaction's state change ---
        if idx == 0 or idx == 1:        # production: add one molecule
            M += 1
        elif idx == 2:                  # degradation: remove one molecule
            M -= 1
        elif idx == 3:                  # binding: consume a dimer, bind promoter
            M -= 2
            b = 1
        else:                           # unbinding: release a dimer, free promoter
            M += 2
            b = 0
        ts.append(t)
        Ms.append(M)
    return np.array(ts), np.array(Ms)

def time_avg_stats(ts, Ms, t_burn):
    # Time-weighted mean/std over t >= t_burn (states persist between events).
    keep = ts >= t_burn
    ts_k = ts[keep]
    Ms_k = Ms[keep].astype(float)
    dt = np.diff(ts_k)                  # duration each state is held
    vals = Ms_k[:-1]                    # value held during each interval
    w = dt / dt.sum()
    mean = np.sum(w * vals)
    var = np.sum(w * (vals - mean) ** 2)
    return mean, np.sqrt(var)

# --- simulation settings ---
k = 0.1
T = 2000.0
t_burn = 500.0                          # discard transient before computing stats
seed = 101

regimes = {
    "self-inhibition": dict(g0=55, g1=5,  kon=0.002, koff=90),
    "no feedback":     dict(g0=30, g1=30, kon=0.002, koff=90),
    "low copy number": dict(g0=5.5, g1=0.5, kon=0.2, koff=90),
}

results = {}
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
for ax, (name, p) in zip(axes, regimes.items()):
    ts, Ms = gillespie(k=k, T=T, seed=seed, **p)
    mean, std = time_avg_stats(ts, Ms, t_burn)
    results[name] = (mean, std)
    ax.step(ts, Ms, where="post", lw=0.7)
    ax.axhline(mean, color="red", ls="--", lw=1)
    ax.set_title(f"{name}: mean={mean:.2f}, std={std:.2f}")
    ax.set_ylabel("molecule count M")
    print(f"{name}: mean = {mean:.4f}")
    print(f"{name}: std  = {std:.4f}")
    print(f"{name}: CV   = {std / mean:.4f}")
axes[-1].set_xlabel("time")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.4.1_s1.png")

# --- separate check ---
std_self = results["self-inhibition"][1]
std_none = results["no feedback"][1]
print(f"CHECK: std(self-inhibition) = {std_self:.4f}")
print(f"CHECK: std(no feedback)     = {std_none:.4f}")
print(f"CHECK: self-inhibition std < no-feedback std ? {std_self < std_none}")

cv_low = results["low copy number"][1] / results["low copy number"][0]
print(f"CHECK: low-copy-number CV (std/mean) = {cv_low:.4f}  (large CV => bursty)")

# Explanation: a lower standard deviation for the self-inhibiting gene than for
# the no-feedback gene at comparable mean directly demonstrates noise suppression,
# because negative autoregulation pushes production down when M is high and up
# when M is low, tightening fluctuations around the mean.
print("EXPLANATION: The self-inhibiting gene's smaller std at a comparable mean shows "
      "negative feedback tightens fluctuations, while the large CV at low copy number "
      "shows few-molecule binding/unbinding drives strongly bursty dynamics.")
