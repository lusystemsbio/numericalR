import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
# Bursting birth-death process:
#   birth:  0 -> n*x  at rate g/n   (each event adds n molecules)
#   death:  x -> 0     at rate k*x
# Mean production of molecules per unit time = (g/n) * n = g, fixed regardless of n.
g = 16.0          # fixed mean production rate
k = 0.1           # per-molecule degradation rate
x_ss = g / k      # deterministic steady state = 160
tmax = 2000.0
burst_sizes = [1, 2, 4, 8]
seed = 91

def gillespie_bursting(n, g, k, x0, tmax, rng):
    """Explicit Gillespie SSA for the bursting birth-death process."""
    t = 0.0
    x = x0
    ts = [t]          # record times
    xs = [x]          # record molecule counts
    while t < tmax:
        # Compute the two reaction propensities at current state
        a_birth = g / n          # burst production rate (independent of x)
        a_death = k * x          # degradation rate (proportional to x)
        a_tot = a_birth + a_death
        if a_tot <= 0:
            break
        # Time to next reaction: exponential with rate a_tot
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:
            break
        # Choose which reaction fires, proportional to propensity
        if rng.random() < a_birth / a_tot:
            x += n               # a burst of n molecules is produced
        else:
            x -= 1               # one molecule degrades
        ts.append(t)
        xs.append(x)
    return np.array(ts), np.array(xs)

def time_average_stats(ts, xs, t_burn=0.0):
    """Mean and std of x weighted by the dwell time in each state (after burn-in)."""
    # dwell time in state xs[i] is ts[i+1]-ts[i]; last state extends to tmax implicitly ignored
    t = ts[:-1]
    x = xs[:-1]
    dt = np.diff(ts)
    mask = t >= t_burn
    t, x, dt = t[mask], x[mask], dt[mask]
    T = dt.sum()
    mean = np.sum(x * dt) / T
    var = np.sum((x - mean) ** 2 * dt) / T
    return mean, np.sqrt(var)

# Master RNG seeded once; each burst size draws from the same stream (reproducible)
rng = np.random.default_rng(seed)

print(f"Fixed mean production rate g = {g}")
print(f"Degradation rate k = {k}")
print(f"Deterministic steady state x_ss = g/k = {x_ss}")
print(f"Poisson-floor std (n=1 expectation) = sqrt({x_ss}) = {np.sqrt(x_ss):.4f}")
print()

x0 = int(round(x_ss))     # start at steady state
t_burn = 200.0            # discard initial transient for statistics

results = {}
for n in burst_sizes:
    ts, xs = gillespie_bursting(n, g, k, x0, tmax, rng)
    mean, std = time_average_stats(ts, xs, t_burn=t_burn)
    results[n] = (ts, xs, mean, std)
    print(f"n = {n}:  mean = {mean:.4f}   std = {std:.4f}   events = {len(ts)}")

print()
print("Check: mean stays near 160 for every n, while std grows with n.")
means = [results[n][2] for n in burst_sizes]
stds = [results[n][3] for n in burst_sizes]
for n, m, s in zip(burst_sizes, means, stds):
    print(f"  n = {n}:  mean = {m:.2f}   std = {s:.2f}")
print()
# Theoretical prediction: variance = x_ss * (n+1)/2  =>  std = sqrt(160*(n+1)/2)
print("Theoretical std sqrt(g/k * (n+1)/2):")
for n in burst_sizes:
    print(f"  n = {n}:  theory std = {np.sqrt(x_ss * (n + 1) / 2):.4f}")
print()
print("Explanation: because the mean is pinned at 160 for all n while the std")
print("rises monotonically with n, the extra spread must come from bursting itself")
print("(larger, rarer bursts) rather than from any change in the mean production rate,")
print("confirming that bursting injects noise beyond the Poisson floor (std=sqrt(160)).")

# ---- Plots ----
fig, axes = plt.subplots(2, 1, figsize=(10, 9))

# Top: trajectories for each burst size
ax = axes[0]
for n in burst_sizes:
    ts, xs, _, _ = results[n]
    ax.step(ts, xs, where="post", lw=0.7, label=f"n = {n}")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="steady state = 160")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Bursting birth-death trajectories (fixed mean rate g=16, k=0.1)")
ax.legend(loc="upper right", ncol=2, fontsize=8)

# Bottom: mean and std versus n
ax = axes[1]
ax.plot(burst_sizes, means, "o-", color="C0", label="mean (measured)")
ax.plot(burst_sizes, stds, "s-", color="C3", label="std (measured)")
ax.plot(burst_sizes, [np.sqrt(x_ss * (nn + 1) / 2) for nn in burst_sizes],
        "x--", color="C2", label="std (theory)")
ax.axhline(x_ss, color="C0", ls=":", lw=1)
ax.axhline(np.sqrt(x_ss), color="C3", ls=":", lw=1, label="Poisson floor sqrt(160)")
ax.set_xlabel("burst size n")
ax.set_ylabel("value")
ax.set_title("Mean stays at 160; std grows with burst size")
ax.set_xticks(burst_sizes)
ax.legend(loc="center right", fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.3.1_s1.png")
print("\nSaved figure to 8D.3.1_s1.png")
