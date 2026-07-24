import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Gillespie SSA for a birth-death process:
#   0 -> x   at rate g        (constitutive production)
#   x -> 0   at rate k*x      (first-order degradation)
# Deterministic steady state: x_bar = g/k
# ---------------------------------------------------------------

def gillespie_birth_death(g, k, x0, tmax, rng):
    """Explicit Gillespie SSA. Returns times[], counts[], and
    time-weighted mean and std of the molecule count."""
    t = 0.0
    x = int(x0)
    times = [t]
    counts = [x]

    # running accumulators for time-weighted (dwell-time) statistics
    # each state x is held for a duration dt, so we weight by dt
    T_total = 0.0          # total accumulated time
    sum_x = 0.0            # sum of x*dt   -> for mean
    sum_x2 = 0.0           # sum of x^2*dt -> for variance

    while t < tmax:
        # 1) compute reaction propensities for current state
        a_birth = g            # production propensity (constant)
        a_death = k * x        # degradation propensity (proportional to x)
        a_total = a_birth + a_death

        if a_total <= 0.0:
            # no reactions possible (x==0 and g==0); jump to tmax
            dt = tmax - t
            T_total += dt
            sum_x += x * dt
            sum_x2 += x * x * dt
            t = tmax
            break

        # 2) draw waiting time to next reaction (exponential)
        dt = rng.exponential(1.0 / a_total)

        # accumulate time-weighted stats for the state we are LEAVING,
        # clipping the last interval so we never weight beyond tmax
        dt_eff = min(dt, tmax - t)
        T_total += dt_eff
        sum_x += x * dt_eff
        sum_x2 += x * x * dt_eff

        t += dt
        if t > tmax:
            break

        # 3) choose which reaction fired, proportional to its propensity
        if rng.random() * a_total < a_birth:
            x += 1             # birth
        else:
            x -= 1             # death

        times.append(t)
        counts.append(x)

    # 4) finalize time-weighted mean and standard deviation
    mean = sum_x / T_total
    var = sum_x2 / T_total - mean * mean
    std = np.sqrt(max(var, 0.0))
    return np.array(times), np.array(counts), mean, std


# ---------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------
k = 0.1
g_values = [100.0, 10.0, 1.0]
x_bars = [g / k for g in g_values]     # 1000, 100, 10
tmax = 100.0
seed = 101

# store measured statistics for the std-vs-mean check
means_ss, stds_ss = [], []   # started at steady state (used for noise check)

# ---------------------------------------------------------------
# Run simulations and plot
# ---------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for i, (g, x_bar) in enumerate(zip(g_values, x_bars)):
    ax = axes.flat[i]

    # --- start AT the deterministic steady state ---
    rng = np.random.default_rng(seed)
    t_ss, c_ss, m_ss, s_ss = gillespie_birth_death(g, k, round(x_bar), tmax, rng)
    means_ss.append(m_ss)
    stds_ss.append(s_ss)

    # --- start FROM zero (relaxation toward steady state) ---
    rng = np.random.default_rng(seed)
    t_z, c_z, m_z, s_z = gillespie_birth_death(g, k, 0, tmax, rng)

    print(f"g={g:g}, k={k:g}, x_bar={x_bar:g}")
    print(f"  from steady state: time-weighted mean = {m_ss:.4f}, std = {s_ss:.4f}, sqrt(x_bar) = {np.sqrt(x_bar):.4f}")
    print(f"  from zero        : time-weighted mean = {m_z:.4f}, std = {s_z:.4f}")
    print(f"  relative noise (std/mean) from steady state = {s_ss/m_ss:.4f}, 1/sqrt(x_bar) = {1/np.sqrt(x_bar):.4f}")

    # trajectories (step plots, since counts are piecewise-constant)
    ax.step(t_ss, c_ss, where='post', lw=0.8, label='from steady state')
    ax.step(t_z, c_z, where='post', lw=0.8, label='from zero')
    # deterministic relaxation curve x_bar*(1 - exp(-k*t))
    tt = np.linspace(0, tmax, 500)
    ax.plot(tt, x_bar * (1 - np.exp(-k * tt)), 'k--', lw=1.5,
            label=r'$\bar x(1-e^{-kt})$')
    ax.axhline(x_bar, color='gray', ls=':', lw=1)
    ax.set_title(f"g={g:g}  (x_bar={x_bar:g})")
    ax.set_xlabel("time")
    ax.set_ylabel("molecule count x")
    ax.legend(fontsize=8)

# ---------------------------------------------------------------
# std vs mean, compared against sqrt(x_bar) (Poisson prediction)
# ---------------------------------------------------------------
ax = axes.flat[3]
means_ss = np.array(means_ss)
stds_ss = np.array(stds_ss)
order = np.argsort(means_ss)
ax.loglog(means_ss[order], stds_ss[order], 'o-', label='measured std')
mgrid = np.logspace(np.log10(means_ss.min())-0.2, np.log10(means_ss.max())+0.2, 100)
ax.loglog(mgrid, np.sqrt(mgrid), 'k--', label=r'$\sqrt{\bar x}$ (Poisson)')
ax.set_xlabel("time-weighted mean")
ax.set_ylabel("time-weighted std")
ax.set_title("std vs mean  (Poisson: std = sqrt(mean))")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/8D.2.1_s4.png", dpi=120)

# ---------------------------------------------------------------
# Numerical summary of the Poisson check
# ---------------------------------------------------------------
print("\n--- Poisson / intrinsic-noise check (started at steady state) ---")
for g, x_bar, m, s in zip(g_values, x_bars, means_ss, stds_ss):
    print(f"g={g:g}: std/sqrt(mean) = {s/np.sqrt(m):.4f}  (should be ~1 for Poisson); "
          f"relative noise {s/m:.4f} vs 1/sqrt(x_bar) {1/np.sqrt(x_bar):.4f}")

# Explanation:
print("\nWhy this confirms the result: for constant birth and first-order death the "
      "stationary distribution is Poisson, whose variance equals its mean, so observing "
      "std ~ sqrt(mean) ~ sqrt(x_bar) confirms that intrinsic noise is Poissonian and "
      "hence relative noise std/mean = 1/sqrt(x_bar) shrinks as the copy number grows.")
