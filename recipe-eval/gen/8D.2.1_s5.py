import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Birth-death process:  0 -> x at rate g (constant production)
#                       x -> 0 at rate k*x (degradation)
# Deterministic steady state: x_bar = g/k
# ---------------------------------------------------------------

k = 0.1
gs = [100.0, 10.0, 1.0]          # production rates
x_bars = [g / k for g in gs]     # steady states: 1000, 100, 10
tmax = 100.0
seed = 101


def gillespie_birth_death(g, k, x0, tmax, rng):
    """Explicit Gillespie SSA for the birth-death process.
    Returns arrays of event times and molecule counts (piecewise-constant)."""
    t = 0.0
    x = x0
    times = [t]
    counts = [x]
    while t < tmax:
        a_birth = g            # propensity of production (constant)
        a_death = k * x        # propensity of degradation (linear in x)
        a_total = a_birth + a_death
        if a_total <= 0.0:     # no possible reaction (x==0 and g==0)
            break
        # 1) time to next reaction: exponential with rate a_total
        tau = rng.exponential(1.0 / a_total)
        t += tau
        if t > tmax:
            break
        # 2) choose which reaction fires, proportional to propensity
        if rng.random() * a_total < a_birth:
            x += 1             # birth
        else:
            x -= 1             # death
        times.append(t)
        counts.append(x)
    return np.array(times), np.array(counts)


def time_weighted_stats(times, counts, tmax):
    """Mean and std of a piecewise-constant trajectory, weighted by the
    dwell time in each state (the correct weighting for SSA output)."""
    dt = np.diff(np.append(times, tmax))     # dwell time in each state
    total = dt.sum()
    mean = np.sum(counts * dt) / total
    var = np.sum((counts - mean) ** 2 * dt) / total
    return mean, np.sqrt(var)


# --- Run simulations from steady state and from zero ---
rng = np.random.default_rng(seed)

results_ss = {}   # started at steady state
results_zero = {} # started at zero

for g, xb in zip(gs, x_bars):
    t_ss, c_ss = gillespie_birth_death(g, k, int(round(xb)), tmax, rng)
    t_z, c_z = gillespie_birth_death(g, k, 0, tmax, rng)
    results_ss[g] = (t_ss, c_ss)
    results_zero[g] = (t_z, c_z)

# --- Compute time-weighted stats (use steady-state runs, past a burn-in) ---
means = []
stds = []
print("=== Time-weighted statistics (started at steady state) ===")
for g, xb in zip(gs, x_bars):
    t_ss, c_ss = results_ss[g]
    # discard a short burn-in to measure fluctuations about steady state
    burn = 10.0
    mask = t_ss >= burn
    if mask.sum() < 2:
        tt, cc = t_ss, c_ss
    else:
        tt, cc = t_ss[mask], c_ss[mask]
        tt = tt - tt[0]
    m, s = time_weighted_stats(tt, cc, tt[-1] if len(tt) else tmax)
    means.append(m)
    stds.append(s)
    print(f"g = {g:6.1f}  x_bar = {xb:7.1f}  mean = {m:9.3f}  std = {s:8.3f}  "
          f"sqrt(x_bar) = {np.sqrt(xb):7.3f}  relative_noise(std/mean) = {s/m:7.4f}")

means = np.array(means)
stds = np.array(stds)
x_bars_arr = np.array(x_bars)

print("\n=== Poisson check: std vs sqrt(x_bar) ===")
for g, xb, s in zip(gs, x_bars, stds):
    print(f"g = {g:6.1f}  measured std = {s:8.3f}  predicted sqrt(x_bar) = {np.sqrt(xb):8.3f}  "
          f"ratio = {s/np.sqrt(xb):7.4f}")

print("\n=== Relative noise falls as 1/sqrt(x_bar) ===")
for g, xb, m, s in zip(gs, x_bars, means, stds):
    print(f"g = {g:6.1f}  std/mean = {s/m:7.4f}  1/sqrt(x_bar) = {1.0/np.sqrt(xb):7.4f}")

# --- Plotting ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (top-left) trajectories started at steady state
ax = axes[0, 0]
for g, xb in zip(gs, x_bars):
    t_ss, c_ss = results_ss[g]
    ax.step(t_ss, c_ss, where="post", lw=0.8, label=f"g={g:.0f}, x_bar={xb:.0f}")
    ax.axhline(xb, color="k", ls=":", lw=0.6)
ax.set_yscale("log")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Trajectories started at steady state")
ax.legend(fontsize=8)

# (top-right) trajectories started at zero, climbing along x_bar*(1-exp(-k t))
ax = axes[0, 1]
tfine = np.linspace(0, tmax, 400)
for g, xb in zip(gs, x_bars):
    t_z, c_z = results_zero[g]
    line, = ax.step(t_z, c_z, where="post", lw=0.8, label=f"g={g:.0f}")
    ax.plot(tfine, xb * (1 - np.exp(-k * tfine)), color=line.get_color(),
            ls="--", lw=1.5)
ax.set_yscale("symlog")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("From zero: climb along x_bar*(1-exp(-k t)) (dashed)")
ax.legend(fontsize=8)

# (bottom-left) std vs mean with sqrt(x_bar) reference
ax = axes[1, 0]
ax.loglog(means, stds, "o", ms=9, label="measured (std vs mean)")
ref = np.logspace(np.log10(means.min()) - 0.2, np.log10(means.max()) + 0.2, 50)
ax.loglog(ref, np.sqrt(ref), "-", label=r"$\sqrt{\bar{x}}$ (Poisson)")
ax.set_xlabel("mean")
ax.set_ylabel("standard deviation")
ax.set_title("std vs mean tracks sqrt(x_bar)")
ax.legend(fontsize=9)

# (bottom-right) relative noise vs 1/sqrt(x_bar)
ax = axes[1, 1]
ax.loglog(x_bars_arr, stds / means, "s", ms=9, label="measured std/mean")
ax.loglog(ref, 1.0 / np.sqrt(ref), "-", label=r"$1/\sqrt{\bar{x}}$")
ax.set_xlabel("x_bar")
ax.set_ylabel("relative noise (std/mean)")
ax.set_title("Relative noise falls as 1/sqrt(x_bar)")
ax.legend(fontsize=9)

# One-sentence explanation of why the check confirms the result:
# Because a birth-death process with constant production and linear
# degradation has a Poisson steady state, its variance equals its mean,
# so std = sqrt(x_bar); confirming std ~ sqrt(x_bar) therefore confirms
# the intrinsic (Poissonian) noise and the 1/sqrt(x_bar) scaling of
# relative noise.
print("\nWhy the check confirms the result:")
print("The birth-death process has a Poisson steady state (variance = mean), "
      "so std = sqrt(x_bar); observing std track sqrt(x_bar) confirms the intrinsic "
      "Poissonian noise and hence relative noise ~ 1/sqrt(x_bar).")

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/8D.2.1_s5.png", dpi=130)
