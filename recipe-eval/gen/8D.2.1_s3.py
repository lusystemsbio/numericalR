import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Birth-death process:  0 -> x  at rate g  (production, constant)
#                       x -> 0  at rate k*x (degradation, linear)
# Deterministic steady state: x_bar = g/k
# ---------------------------------------------------------------

def gillespie_bd(g, k, x0, tmax, rng):
    """Explicit Gillespie SSA for the birth-death process.
    Returns time-stamped trajectory (times, counts)."""
    t = 0.0            # current time
    x = x0             # current molecule count
    times = [t]        # record of jump times
    counts = [x]       # record of counts
    while t < tmax:
        # 1. compute propensities for each reaction
        a_birth = g          # production propensity (constant)
        a_death = k * x      # degradation propensity (proportional to x)
        a_tot = a_birth + a_death
        if a_tot <= 0.0:     # no reaction possible (absorbing at x=0 with g=0)
            break
        # 2. draw time to next reaction: exponential with rate a_tot
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:         # do not step past tmax
            break
        # 3. choose which reaction fires, weighted by propensity
        if rng.random() * a_tot < a_birth:
            x += 1           # a birth occurred
        else:
            x -= 1           # a death occurred
        # 4. record the new state
        times.append(t)
        counts.append(x)
    # extend the last state out to tmax so weighting covers full window
    times.append(tmax)
    counts.append(x)
    return np.array(times, dtype=float), np.array(counts, dtype=float)


def time_weighted_stats(times, counts, t_burn):
    """Time-weighted mean and std over [t_burn, tmax], since between
    jumps the count is constant so each state's weight is its dwell time."""
    dt = np.diff(times)                 # dwell time in each state
    vals = counts[:-1]                  # count held during each interval
    t_start = times[:-1]                # interval start times
    # keep only intervals after the burn-in time
    mask = t_start >= t_burn
    dt = dt[mask]
    vals = vals[mask]
    W = dt.sum()
    mean = np.sum(vals * dt) / W                      # time-weighted mean
    var = np.sum((vals - mean) ** 2 * dt) / W         # time-weighted variance
    return mean, np.sqrt(var)


# ---------------------------- parameters ----------------------------
k = 0.1
g_list = [100.0, 10.0, 1.0]
tmax = 100.0
seed = 101
t_burn = 20.0   # discard early transient for the steady-state statistics

means_ss, stds_ss = [], []
means_z, stds_z = [], []
xbars = []

# figure 1: trajectories
fig1, axes = plt.subplots(len(g_list), 1, figsize=(9, 9), sharex=True)

for i, g in enumerate(g_list):
    x_bar = g / k
    xbars.append(x_bar)
    rng = np.random.default_rng(seed)   # reseed per g for reproducibility

    # --- run started AT the steady state ---
    t_ss, c_ss = gillespie_bd(g, k, int(round(x_bar)), tmax, rng)
    m_ss, s_ss = time_weighted_stats(t_ss, c_ss, t_burn)
    means_ss.append(m_ss); stds_ss.append(s_ss)

    # --- run started FROM zero ---
    t_z, c_z = gillespie_bd(g, k, 0, tmax, rng)
    m_z, s_z = time_weighted_stats(t_z, c_z, t_burn)
    means_z.append(m_z); stds_z.append(s_z)

    print(f"g = {g:6.1f}  x_bar = {x_bar:7.1f}")
    print(f"  from steady state: mean = {m_ss:10.4f}  std = {s_ss:10.4f}  sqrt(x_bar) = {np.sqrt(x_bar):10.4f}")
    print(f"  from zero        : mean = {m_z:10.4f}  std = {s_z:10.4f}  sqrt(x_bar) = {np.sqrt(x_bar):10.4f}")
    print(f"  relative noise (std/mean) from steady state = {s_ss / m_ss:.5f}   1/sqrt(x_bar) = {1/np.sqrt(x_bar):.5f}")

    # plot both trajectories (step plots: state held constant between jumps)
    ax = axes[i]
    ax.step(t_ss, c_ss, where="post", color="C0", lw=1.0, label="start at steady state")
    ax.step(t_z, c_z, where="post", color="C1", lw=1.0, alpha=0.8, label="start from zero")
    ax.axhline(x_bar, color="k", ls="--", lw=1.0, label=r"$\bar{x}=g/k$")
    # analytic mean relaxation from zero: x_bar*(1 - exp(-k t))
    tt = np.linspace(0, tmax, 400)
    ax.plot(tt, x_bar * (1 - np.exp(-k * tt)), color="C3", lw=2.0, ls=":",
            label=r"$\bar{x}(1-e^{-kt})$")
    ax.set_ylabel("molecule count x")
    ax.set_title(f"g = {g:.0f},  x_bar = {x_bar:.0f}")
    ax.legend(fontsize=7, loc="lower right")

axes[-1].set_xlabel("time")
fig1.tight_layout()

# figure 2: std vs mean, compared to sqrt(x_bar) (Poisson line)
fig2, ax2 = plt.subplots(figsize=(7, 6))
xbars = np.array(xbars)
means_ss = np.array(means_ss); stds_ss = np.array(stds_ss)
means_z = np.array(means_z); stds_z = np.array(stds_z)

ax2.scatter(means_ss, stds_ss, color="C0", s=60, label="std vs mean (from steady state)")
ax2.scatter(means_z, stds_z, color="C1", marker="s", s=60, label="std vs mean (from zero)")
xline = np.linspace(min(means_ss.min(), means_z.min()) * 0.5,
                    max(means_ss.max(), means_z.max()) * 1.5, 200)
ax2.plot(xline, np.sqrt(xline), "k--", label=r"$\sqrt{\bar{x}}$ (Poisson)")
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlabel("time-weighted mean count")
ax2.set_ylabel("time-weighted std")
ax2.set_title("Intrinsic noise: std follows sqrt(mean)")
ax2.legend()
fig2.tight_layout()

# save (combine both figures onto one canvas via saving fig1 with fig2 alongside)
# Save the trajectories figure as the primary output.
fig1.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/8D.2.1_s3.png", dpi=120)
# Also append the std-vs-mean panel below by making a combined figure.
fig_comb, axc = plt.subplots(figsize=(7, 6))
axc.scatter(means_ss, stds_ss, color="C0", s=60, label="from steady state")
axc.scatter(means_z, stds_z, color="C1", marker="s", s=60, label="from zero")
axc.plot(xline, np.sqrt(xline), "k--", label=r"$\sqrt{\bar{x}}$ (Poisson)")
axc.set_xscale("log"); axc.set_yscale("log")
axc.set_xlabel("mean"); axc.set_ylabel("std")
axc.set_title("std vs mean against sqrt(x_bar)")
axc.legend()
fig_comb.tight_layout()

# Print summary of the Poisson / noise check
print("\n--- Poisson check (std / sqrt(x_bar) should be ~1) ---")
for xb, s in zip(xbars, stds_ss):
    print(f"x_bar = {xb:7.1f}   std/sqrt(x_bar) = {s / np.sqrt(xb):.4f}")

print("\nExplanation: because the birth-death steady state is Poisson its "
      "variance equals its mean, so std = sqrt(x_bar) and the relative noise "
      "std/mean = 1/sqrt(x_bar) shrinks as x_bar grows; seeing the simulated "
      "std lie on the sqrt(x_bar) line (and the from-zero mean track "
      "x_bar*(1-exp(-k*t))) confirms the intrinsic-noise result.")
