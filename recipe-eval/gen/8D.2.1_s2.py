import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model: birth-death 0 -> x at rate g, x -> 0 at rate k*x ---
# Deterministic steady state x_bar = g/k. Poisson prediction: std = sqrt(x_bar).

k = 0.1
g_values = [100.0, 10.0, 1.0]      # production rates
x_bar_values = [g / k for g in g_values]  # 1000, 100, 10
tmax = 100.0
seed = 101


def gillespie_birth_death(g, k, x0, tmax, rng):
    """Explicit Gillespie SSA for the birth-death process.

    Returns arrays of jump times and molecule counts (step trajectory).
    """
    t = 0.0
    x = x0
    times = [t]
    counts = [x]
    while t < tmax:
        # Propensities for the two reactions
        a_birth = g          # 0 -> x, constant
        a_death = k * x      # x -> 0, proportional to x
        a_total = a_birth + a_death
        if a_total <= 0.0:    # no possible reaction (only if x==0 and g==0)
            break
        # Time to next reaction: exponential with rate a_total
        tau = rng.exponential(1.0 / a_total)
        t += tau
        if t > tmax:          # do not step past the observation window
            break
        # Choose which reaction fires, proportional to its propensity
        if rng.random() < a_birth / a_total:
            x += 1            # birth
        else:
            x -= 1            # death
        times.append(t)
        counts.append(x)
    # Append final segment out to tmax so the last state has proper duration
    times.append(tmax)
    counts.append(counts[-1])
    return np.array(times, dtype=float), np.array(counts, dtype=float)


def time_weighted_stats(times, counts):
    """Time-weighted mean and std of a piecewise-constant trajectory.

    Each state 'counts[i]' persists for duration dt[i] = times[i+1]-times[i];
    weight the moments by these durations (not by number of jumps).
    """
    dt = np.diff(times)                 # duration of each held state
    vals = counts[:-1]                  # value held during each interval
    T = dt.sum()
    mean = np.sum(vals * dt) / T        # time-weighted mean
    var = np.sum((vals - mean) ** 2 * dt) / T  # time-weighted variance
    return mean, np.sqrt(var)


# Use a single seeded generator so the whole run is reproducible
rng = np.random.default_rng(seed)

results_ss = []    # started at steady state
results_zero = []  # started from zero
traj_ss = {}
traj_zero = {}

for g, x_bar in zip(g_values, x_bar_values):
    # Start at deterministic steady state
    t_ss, x_ss = gillespie_birth_death(g, k, int(round(x_bar)), tmax, rng)
    m_ss, s_ss = time_weighted_stats(t_ss, x_ss)
    results_ss.append((x_bar, m_ss, s_ss))
    traj_ss[g] = (t_ss, x_ss)

    # Start from zero (transient climb toward x_bar)
    t_z, x_z = gillespie_birth_death(g, k, 0, tmax, rng)
    m_z, s_z = time_weighted_stats(t_z, x_z)
    results_zero.append((x_bar, m_z, s_z))
    traj_zero[g] = (t_z, x_z)

# --- Print numerical results ---
print("=== Started at steady state ===")
for (x_bar, m, s) in results_ss:
    print(f"x_bar={x_bar:7.1f} | time-weighted mean={m:9.3f} | std={s:8.3f} | "
          f"sqrt(x_bar)={np.sqrt(x_bar):8.3f} | rel_noise(std/mean)={s/m:.4f} | "
          f"1/sqrt(x_bar)={1/np.sqrt(x_bar):.4f}")

print("\n=== Started from zero ===")
for (x_bar, m, s) in results_zero:
    print(f"x_bar={x_bar:7.1f} | time-weighted mean={m:9.3f} | std={s:8.3f} | "
          f"sqrt(x_bar)={np.sqrt(x_bar):8.3f}")

# --- Plots ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (a) Trajectories started at steady state
ax = axes[0, 0]
for g, x_bar in zip(g_values, x_bar_values):
    t, x = traj_ss[g]
    ax.step(t, x, where="post", lw=0.8, label=f"g={g:g}, x_bar={x_bar:g}")
    ax.axhline(x_bar, color="k", ls=":", lw=0.6)
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_yscale("log")
ax.set_title("Trajectories from steady state")
ax.legend(fontsize=8)

# (b) Trajectories started from zero, with x_bar*(1-exp(-k t)) overlay
ax = axes[0, 1]
tt = np.linspace(0, tmax, 400)
for g, x_bar in zip(g_values, x_bar_values):
    t, x = traj_zero[g]
    line, = ax.step(t, x, where="post", lw=0.8, label=f"g={g:g}")
    ax.plot(tt, x_bar * (1 - np.exp(-k * tt)), color=line.get_color(),
            ls="--", lw=1.5)
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("From zero: climb along x_bar*(1-exp(-k t)) (dashed)")
ax.legend(fontsize=8)

# (c) std vs mean compared to sqrt(x_bar)
ax = axes[1, 0]
means = [m for (_, m, _) in results_ss]
stds = [s for (_, _, s) in results_ss]
ax.loglog(means, stds, "o", ms=8, label="Gillespie std (from steady state)")
xline = np.linspace(min(x_bar_values) * 0.7, max(x_bar_values) * 1.3, 100)
ax.loglog(xline, np.sqrt(xline), "-", label="sqrt(x_bar) (Poisson)")
ax.set_xlabel("mean count")
ax.set_ylabel("std of count")
ax.set_title("std vs mean tracks sqrt(x_bar)")
ax.legend(fontsize=8)

# (d) relative noise vs x_bar compared to 1/sqrt(x_bar)
ax = axes[1, 1]
rel = [s / m for (_, m, s) in results_ss]
ax.loglog(x_bar_values, rel, "o", ms=8, label="std/mean (Gillespie)")
ax.loglog(xline, 1 / np.sqrt(xline), "-", label="1/sqrt(x_bar)")
ax.set_xlabel("x_bar")
ax.set_ylabel("relative noise (std/mean)")
ax.set_title("Relative noise falls as 1/sqrt(x_bar)")
ax.legend(fontsize=8)

# One-sentence explanation of why the check confirms the result:
# Because a constitutively transcribed gene's birth-death process has a
# Poisson steady state, so std = sqrt(x_bar) exactly; observing the simulated
# std lie on the sqrt(x_bar) line (and relative noise on 1/sqrt(x_bar))
# confirms the intrinsic noise is Poissonian.
print("\nWhy the check confirms it: the birth-death steady state is Poisson, "
      "so std=sqrt(mean); the simulated std landing on the sqrt(x_bar) line "
      "(and relative noise on 1/sqrt(x_bar)) confirms intrinsic noise is Poissonian.")

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/8D.2.1_s2.png", dpi=130)
