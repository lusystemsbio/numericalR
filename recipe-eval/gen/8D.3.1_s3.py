import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Bursting birth-death process:
#   production:  0 -> n*x   at rate g/n   (each event adds n molecules)
#   decay:       x -> 0     at rate k*x
# Mean production rate = (g/n) * n = g, held fixed as n varies.
# ----------------------------------------------------------------------

def gillespie_bursting(g, k, n, x0, tmax, rng):
    """Explicit Gillespie SSA for the bursting birth-death model."""
    t = 0.0                       # current time
    x = x0                        # current molecule count
    ts = [t]                      # recorded times
    xs = [x]                      # recorded states
    while t < tmax:
        # 1. compute the two reaction propensities
        a_prod = g / n            # burst production propensity (adds n)
        a_dec = k * x             # decay propensity (removes 1)
        a_tot = a_prod + a_dec
        if a_tot <= 0:            # no reaction possible
            break
        # 2. time to next event: exponential with rate a_tot
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:              # do not step past tmax
            break
        # 3. pick which reaction fires, proportional to propensity
        if rng.random() < a_prod / a_tot:
            x += n                # production burst
        else:
            x -= 1                # single-molecule decay
        # 4. record the new state
        ts.append(t)
        xs.append(x)
    return np.array(ts), np.array(xs)

def time_average_moments(ts, xs):
    """Mean and std of x weighted by the time spent in each state."""
    dt = np.diff(ts)              # dwell time in each state
    xseg = xs[:-1]                # state held during each interval
    T = dt.sum()
    mean = np.sum(xseg * dt) / T
    var = np.sum((xseg - mean)**2 * dt) / T
    return mean, np.sqrt(var)

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
g = 16
k = 0.1
x_ss = int(round(g / k))          # steady-state mean = 160
burst_sizes = [1, 2, 4, 8]
tmax = 2000
seed = 91

rng = np.random.default_rng(seed)

results = {}                      # n -> (ts, xs, mean, std)
for n in burst_sizes:
    ts, xs = gillespie_bursting(g, k, n, x_ss, tmax, rng)
    mean, std = time_average_moments(ts, xs)
    results[n] = (ts, xs, mean, std)

# ----------------------------------------------------------------------
# Print numerical results
# ----------------------------------------------------------------------
print(f"Steady-state mean (g/k): {x_ss}")
for n in burst_sizes:
    _, _, mean, std = results[n]
    print(f"n = {n}: time-averaged mean = {mean:.3f}")
for n in burst_sizes:
    _, _, mean, std = results[n]
    print(f"n = {n}: time-averaged std  = {std:.3f}")

# Poisson floor for comparison: std = sqrt(mean) at n = 1
poisson_std = np.sqrt(x_ss)
print(f"Poisson floor std (sqrt(160)): {poisson_std:.3f}")

# Theoretical std for bursting: var = mean*(n+1)/2  ->  std = sqrt(160*(n+1)/2)
for n in burst_sizes:
    theo_std = np.sqrt(x_ss * (n + 1) / 2.0)
    print(f"n = {n}: theoretical std = {theo_std:.3f}")

# ----------------------------------------------------------------------
# Plot: trajectories + mean/std versus n
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (a) trajectories
ax = axes[0, 0]
for n in burst_sizes:
    ts, xs, _, _ = results[n]
    ax.step(ts, xs, where="post", lw=0.8, label=f"n = {n}")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="mean = 160")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Bursting birth-death trajectories")
ax.legend(fontsize=8)

# (b) mean versus n
ax = axes[0, 1]
means = [results[n][2] for n in burst_sizes]
ax.plot(burst_sizes, means, "o-", color="tab:blue")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="expected 160")
ax.set_xlabel("burst size n")
ax.set_ylabel("mean count")
ax.set_title("Mean stays fixed at 160")
ax.set_ylim(0, 220)
ax.legend(fontsize=8)

# (c) std versus n
ax = axes[1, 0]
stds = [results[n][3] for n in burst_sizes]
theo = [np.sqrt(x_ss * (n + 1) / 2.0) for n in burst_sizes]
ax.plot(burst_sizes, stds, "o-", color="tab:red", label="simulated std")
ax.plot(burst_sizes, theo, "s--", color="tab:green", label="theory sqrt(160(n+1)/2)")
ax.axhline(poisson_std, color="k", ls=":", lw=1, label="Poisson floor sqrt(160)")
ax.set_xlabel("burst size n")
ax.set_ylabel("standard deviation")
ax.set_title("Noise grows with burst size")
ax.legend(fontsize=8)

# (d) histograms of the states
ax = axes[1, 1]
for n in burst_sizes:
    ts, xs, _, _ = results[n]
    dt = np.diff(ts)
    ax.hist(xs[:-1], bins=40, weights=dt, density=True, histtype="step",
            lw=1.2, label=f"n = {n}")
ax.axvline(x_ss, color="k", ls="--", lw=1)
ax.set_xlabel("molecule count x")
ax.set_ylabel("time-weighted density")
ax.set_title("Distribution widens with n")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.3.1_s3.png", dpi=130)

# ----------------------------------------------------------------------
# One-sentence explanation
# ----------------------------------------------------------------------
print("Explanation: Because the mean stays pinned at 160 for every n while the "
      "standard deviation rises with n, the extra spread cannot come from a "
      "change in average expression and must instead be the added burst noise, "
      "confirming that larger, rarer bursts inject variance beyond the Poisson "
      "floor sqrt(160).")
