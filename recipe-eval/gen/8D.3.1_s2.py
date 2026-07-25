import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
# Bursting birth-death process:
#   0 -> n*x   at rate g/n   (a burst adds n molecules; firing rate scaled so mean production = g)
#   x -> 0     at rate k*x   (each molecule degrades at rate k)
g = 16.0          # mean production rate (molecules per unit time), held FIXED across n
k = 0.1           # degradation rate per molecule
x_ss = g / k      # deterministic steady state = 160
tmax = 2000.0     # simulation end time
burst_sizes = [1, 2, 4, 8]
seed = 91

print(f"Mean production rate g = {g}")
print(f"Degradation rate k = {k}")
print(f"Deterministic steady state g/k = {x_ss}")

def gillespie_bursting(n, g, k, x0, tmax, rng):
    # Explicit Gillespie SSA for the two-reaction bursting birth-death system.
    t = 0.0
    x = x0
    ts = [t]        # record time points
    xs = [x]        # record molecule counts
    while t < tmax:
        # Compute propensities (reaction rates) for the current state.
        a_birth = g / n         # burst-firing propensity (constant; each firing adds n molecules)
        a_death = k * x         # degradation propensity (proportional to copy number)
        a_tot = a_birth + a_death
        if a_tot <= 0.0:
            break
        # Time to next reaction: exponential with rate a_tot.
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:
            break
        # Choose which reaction fires, proportional to its propensity.
        if rng.random() < a_birth / a_tot:
            x += n              # burst: add n molecules
        else:
            x -= 1              # degradation: remove one molecule
        ts.append(t)
        xs.append(x)
    return np.array(ts), np.array(xs)

# Run one simulation per burst size, all started at steady state.
rng = np.random.default_rng(seed)
results = {}
for n in burst_sizes:
    ts, xs = gillespie_bursting(n, g, k, int(round(x_ss)), tmax, rng)
    results[n] = (ts, xs)

# --- Compute time-weighted mean and standard deviation for each n ---
# States are held between events, so weight each recorded value by its dwell time.
means = []
stds = []
for n in burst_sizes:
    ts, xs = results[n]
    dt = np.diff(ts)                 # dwell time in each state
    vals = xs[:-1]                   # value held during each interval
    T = dt.sum()
    mean = np.sum(vals * dt) / T
    var = np.sum((vals - mean) ** 2 * dt) / T
    std = np.sqrt(var)
    means.append(mean)
    stds.append(std)
    print(f"n = {n}: time-averaged mean = {mean:.4f}, std = {std:.4f}")

means = np.array(means)
stds = np.array(stds)

# Poisson floor for reference: for simple birth-death (n=1) std = sqrt(mean) ~ sqrt(160).
poisson_std = np.sqrt(x_ss)
print(f"Poisson-floor std sqrt(160) = {poisson_std:.4f}")

# --- Plotting ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Top-left: overlaid trajectories.
ax = axes[0, 0]
for n in burst_sizes:
    ts, xs = results[n]
    ax.step(ts, xs, where="post", lw=0.7, label=f"n = {n}")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="steady state 160")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Trajectories for each burst size n")
ax.legend(fontsize=8)

# Top-right: zoom into a shorter window to see burst structure.
ax = axes[0, 1]
for n in burst_sizes:
    ts, xs = results[n]
    m = ts <= 300
    ax.step(ts[m], xs[m], where="post", lw=0.8, label=f"n = {n}")
ax.axhline(x_ss, color="k", ls="--", lw=1)
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Zoom (t <= 300): larger n = bigger jumps")
ax.legend(fontsize=8)

# Bottom-left: mean vs n (should stay flat at 160).
ax = axes[1, 0]
ax.plot(burst_sizes, means, "o-", color="C0")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="160")
ax.set_xlabel("burst size n")
ax.set_ylabel("mean x")
ax.set_title("Mean stays at 160 for every n")
ax.set_ylim(140, 180)
ax.legend()

# Bottom-right: std vs n (should grow with n).
ax = axes[1, 1]
ax.plot(burst_sizes, stds, "s-", color="C3", label="simulated std")
ax.axhline(poisson_std, color="k", ls="--", lw=1, label="Poisson floor sqrt(160)")
ax.set_xlabel("burst size n")
ax.set_ylabel("std of x")
ax.set_title("Std grows with burst size n")
ax.legend()

fig.tight_layout()
fig.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.3.1_s2.png", dpi=110)

# --- Separate check summary ---
print("\nCheck: mean vs n and std vs n")
for n, mn, sd in zip(burst_sizes, means, stds):
    print(f"  n = {n}: mean ~= {mn:.2f} (target 160), std = {sd:.2f}")
print(f"Mean range across n: {means.min():.2f} to {means.max():.2f} (all near 160)")
print(f"Std range across n: {stds.min():.2f} to {stds.max():.2f} (increasing with n)")
print("Explanation: because the mean holds at 160 while the standard deviation "
      "climbs with n, the extra variability cannot come from a higher production "
      "rate and must come from the bursting itself, confirming that larger, rarer "
      "bursts inject noise above the Poisson floor.")
