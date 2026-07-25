import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Bursting birth-death process via explicit Gillespie SSA
#   Production:   0 -> n*x   at rate g/n   (adds n molecules per event)
#   Degradation:  x -> 0     at rate k*x   (each molecule dies at rate k)
# Because the production event rate is g/n and each adds n molecules,
# the MEAN production flux is (g/n)*n = g, independent of burst size n.
# ----------------------------------------------------------------------

g = 16.0        # mean production rate (molecules per unit time)
k = 0.1         # per-molecule degradation rate
x_ss = g / k    # deterministic steady state = 160
tmax = 2000.0
burst_sizes = [1, 2, 4, 8]
seed = 91

def gillespie_burst(n, g, k, x0, tmax, rng):
    # Explicit stochastic simulation of the bursting birth-death process.
    t = 0.0
    x = x0
    times = [t]            # event times
    states = [x]           # molecule count after each event
    dwell_x = []           # molecule count held during each interval
    dwell_dt = []          # duration of each interval (for time-weighted stats)
    while t < tmax:
        a_birth = g / n            # propensity of a burst production event
        a_death = k * x            # propensity of a degradation event
        a_total = a_birth + a_death
        if a_total <= 0.0:         # absorbing (x=0 and no births) - shouldn't happen here
            break
        # 1) draw waiting time to next event from exponential(a_total)
        dt = -np.log(rng.random()) / a_total
        # record the state held over this interval BEFORE the event fires
        dwell_x.append(x)
        # clip the last interval so it does not run past tmax
        dwell_dt.append(min(dt, tmax - t))
        t += dt
        if t > tmax:
            break
        # 2) pick which reaction fires, proportional to its propensity
        if rng.random() * a_total < a_birth:
            x += n                 # a burst: add n molecules at once
        else:
            x -= 1                 # one molecule degrades
        times.append(t)
        states.append(x)
    return np.array(times), np.array(states), np.array(dwell_x), np.array(dwell_dt)

# time-weighted (dwell-time) mean and standard deviation of the trajectory
def time_weighted_stats(dwell_x, dwell_dt):
    T = dwell_dt.sum()
    mean = np.sum(dwell_x * dwell_dt) / T
    var = np.sum(dwell_x**2 * dwell_dt) / T - mean**2
    return mean, np.sqrt(var)

rng = np.random.default_rng(seed)

results = {}
for n in burst_sizes:
    times, states, dwell_x, dwell_dt = gillespie_burst(n, g, k, x_ss, tmax, rng)
    mean, std = time_weighted_stats(dwell_x, dwell_dt)
    results[n] = dict(times=times, states=states, mean=mean, std=std)

# ----------------------------------------------------------------------
# Print numerical results
# ----------------------------------------------------------------------
print(f"Parameters: g = {g}, k = {k}, steady state g/k = {x_ss}, tmax = {tmax}, seed = {seed}")
means = []
stds = []
for n in burst_sizes:
    m = results[n]["mean"]
    s = results[n]["std"]
    means.append(m)
    stds.append(s)
    # theoretical Fano factor for fixed-size bursts is (n+1)/2, so std_theory = sqrt(mean*(n+1)/2)
    fano_theory = (n + 1) / 2.0
    std_theory = np.sqrt(x_ss * fano_theory)
    print(f"burst n = {n}:  time-weighted mean = {m:.3f}   std = {s:.3f}   "
          f"(theory std ~ {std_theory:.3f}, Fano = {fano_theory:.2f})")

print("\nCheck: mean should stay ~160 for every n (fixed mean production g)")
for n, m in zip(burst_sizes, means):
    print(f"  n = {n}: mean = {m:.3f}")
print("\nCheck: standard deviation should grow with burst size n")
for n, s in zip(burst_sizes, stds):
    print(f"  n = {n}: std  = {s:.3f}")
print(f"\nPoisson floor (n=1) std = sqrt(160) = {np.sqrt(x_ss):.3f}")
print("Explanation: because the mean is pinned at 160 for all n while the "
      "standard deviation rises with n, the extra spread cannot come from a "
      "changed production level, so it must come from the burstiness itself "
      "(larger, rarer jumps), proving bursting adds noise beyond the Poisson floor.")

# ----------------------------------------------------------------------
# Plots: trajectories (top) and mean/std vs n (bottom)
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8))

ax1 = fig.add_subplot(2, 1, 1)
for n in burst_sizes:
    r = results[n]
    ax1.step(r["times"], r["states"], where="post", lw=0.7, label=f"n = {n}")
ax1.axhline(x_ss, color="k", ls="--", lw=1, label="steady state 160")
ax1.set_xlabel("time")
ax1.set_ylabel("molecules x(t)")
ax1.set_title("Bursting birth-death trajectories (Gillespie SSA), fixed mean rate g=16")
ax1.legend(loc="upper right", ncol=3, fontsize=8)

ax2 = fig.add_subplot(2, 2, 3)
ax2.plot(burst_sizes, means, "o-", color="C0")
ax2.axhline(x_ss, color="k", ls="--", lw=1)
ax2.set_xlabel("burst size n")
ax2.set_ylabel("mean")
ax2.set_title("Mean stays at 160")
ax2.set_ylim(0, 220)

ax3 = fig.add_subplot(2, 2, 4)
ax3.plot(burst_sizes, stds, "s-", color="C3", label="simulated std")
ax3.plot(burst_sizes, [np.sqrt(x_ss * (n + 1) / 2.0) for n in burst_sizes],
         "x--", color="gray", label="theory sqrt(160*(n+1)/2)")
ax3.set_xlabel("burst size n")
ax3.set_ylabel("standard deviation")
ax3.set_title("Std grows with burst size")
ax3.legend(fontsize=8)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.3.1_s5.png", dpi=120)
