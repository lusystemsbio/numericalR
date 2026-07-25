import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model parameters -----
# Bursting birth-death process:
#   0 -> n*x  at rate g/n   (production: bursts of size n, fixed mean rate g)
#   x -> 0    at rate k*x    (degradation)
g = 16.0        # mean production rate (molecules per unit time)
k = 0.1         # degradation rate per molecule
x_ss = g / k    # deterministic steady state = 160
tmax = 2000.0
burst_sizes = [1, 2, 4, 8]
seed = 91

print(f"Deterministic steady state x_ss = g/k = {x_ss:.1f}")

def gillespie_bursting(n, g, k, x0, tmax, rng):
    """Explicit Gillespie SSA for the bursting birth-death process."""
    t = 0.0
    x = x0
    ts = [t]          # record times
    xs = [x]          # record molecule counts
    while t < tmax:
        # Compute the two reaction propensities
        a_prod = g / n          # production propensity (constant, burst of size n)
        a_deg  = k * x          # degradation propensity (proportional to x)
        a_tot  = a_prod + a_deg
        if a_tot <= 0.0:
            break
        # Draw time to next event from an exponential distribution
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:
            break
        # Choose which reaction fires, weighted by propensity
        if rng.random() < a_prod / a_tot:
            x += n              # production: add n molecules
        else:
            x -= 1              # degradation: remove one molecule
        ts.append(t)
        xs.append(x)
    return np.array(ts), np.array(xs)

def time_average_stats(ts, xs, t_burn=200.0):
    """Time-weighted mean and std of a piecewise-constant trajectory, after burn-in."""
    # Each recorded state xs[i] persists over [ts[i], ts[i+1]); weight by duration.
    dt = np.diff(ts)
    xv = xs[:-1]
    tv = ts[:-1]
    mask = tv >= t_burn         # discard initial burn-in interval
    w = dt[mask]
    x = xv[mask]
    W = w.sum()
    mean = np.sum(w * x) / W
    var  = np.sum(w * (x - mean) ** 2) / W
    return mean, np.sqrt(var)

rng = np.random.default_rng(seed)

# ----- Run one simulation per burst size, starting at steady state -----
trajectories = {}
means = []
stds = []
for n in burst_sizes:
    ts, xs = gillespie_bursting(n, g, k, int(round(x_ss)), tmax, rng)
    trajectories[n] = (ts, xs)
    m, s = time_average_stats(ts, xs)
    means.append(m)
    stds.append(s)
    print(f"n = {n}: mean = {m:.3f}, std = {s:.3f}, Poisson floor sqrt(mean) = {np.sqrt(m):.3f}")

# ----- Theoretical expectation -----
# For this bursting model the stationary variance grows with burst size:
#   var = x_ss * (1 + n) / 2  ->  std = sqrt(x_ss*(1+n)/2)
for n in burst_sizes:
    print(f"n = {n}: theoretical std = sqrt(x_ss*(1+n)/2) = {np.sqrt(x_ss*(1+n)/2):.3f}")

# ----- Plotting -----
fig, axes = plt.subplots(2, 1, figsize=(10, 9))

# Top: trajectories for each burst size
ax = axes[0]
for n in burst_sizes:
    ts, xs = trajectories[n]
    ax.step(ts, xs, where="post", lw=0.7, label=f"n = {n}")
ax.axhline(x_ss, color="k", ls="--", lw=1, label="steady state 160")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Bursting birth-death trajectories (fixed mean rate g=16, k=0.1)")
ax.legend(loc="upper right", ncol=2, fontsize=8)

# Bottom: mean and std vs n
ax = axes[1]
ax.plot(burst_sizes, means, "o-", label="mean")
ax.plot(burst_sizes, stds, "s-", label="std (simulated)")
ax.plot(burst_sizes, [np.sqrt(x_ss*(1+n)/2) for n in burst_sizes], "x--",
        color="gray", label="std (theory)")
ax.axhline(x_ss, color="k", ls=":", lw=1)
ax.axhline(np.sqrt(x_ss), color="r", ls=":", lw=1, label="Poisson floor sqrt(160)")
ax.set_xlabel("burst size n")
ax.set_ylabel("value")
ax.set_title("Mean stays fixed while std grows with burst size")
ax.set_xticks(burst_sizes)
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.3.1_s4.png")

# ----- Summary of the check -----
print("\nSummary (mean +/- std by burst size):")
for n, m, s in zip(burst_sizes, means, stds):
    print(f"n = {n}: mean = {m:.2f}, std = {s:.2f}")
print("Poisson floor (std if Poisson) = sqrt(160) =", f"{np.sqrt(x_ss):.3f}")

# One-sentence explanation:
print("\nExplanation: Because the mean holds at ~160 for every n while the standard "
      "deviation rises well above the Poisson value sqrt(160)~12.6 as n grows, the "
      "extra spread cannot come from the mean production rate but only from the "
      "bursty, rarer-but-larger production events, confirming that bursting injects "
      "noise beyond the Poisson floor.")
