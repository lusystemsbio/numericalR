import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Gillespie SSA for birth-death: 0 -> x at rate g, x -> 0 at rate k*x ---
def gillespie_bd(g, k, tmax, x0, rng):
    t = 0.0
    x = x0
    ts = [t]          # event times
    xs = [x]          # molecule count right after each event
    while t < tmax:
        a_prod = g          # propensity of production (constant)
        a_deg = k * x       # propensity of degradation (proportional to x)
        a_tot = a_prod + a_deg
        if a_tot <= 0:
            break
        # time to next reaction: exponential with rate a_tot
        tau = rng.exponential(1.0 / a_tot)
        t += tau
        if t > tmax:
            break
        # choose which reaction fires, weighted by propensities
        if rng.random() < a_prod / a_tot:
            x += 1          # production
        else:
            x -= 1          # degradation
        ts.append(t)
        xs.append(x)
    return np.array(ts), np.array(xs)

# --- time-weighted mean and std over the trajectory (piecewise-constant x) ---
def time_weighted_stats(ts, xs, tmax, burn_frac=0.0):
    # append tmax as final boundary so the last state has a dwell time
    t_edges = np.append(ts, tmax)
    dwell = np.diff(t_edges)          # time spent in each state xs[i]
    x_vals = xs                        # state held during each dwell interval
    # optional burn-in: drop early transient by start time
    t_start = burn_frac * tmax
    mask = t_edges[:-1] >= t_start
    dwell = dwell[mask]
    x_vals = x_vals[mask]
    W = dwell.sum()
    mean = np.sum(dwell * x_vals) / W
    var = np.sum(dwell * (x_vals - mean) ** 2) / W
    return mean, np.sqrt(var)

k = 0.1
g_list = [100.0, 10.0, 1.0]
tmax = 100.0
seed = 101

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

means_ss, stds_ss, xbars = [], [], []

# ---- Trajectories started at steady state (top-left) + collect stats ----
ax = axes[0, 0]
for g in g_list:
    x_bar = g / k
    xbars.append(x_bar)
    rng = np.random.default_rng(seed)          # reproducible per condition
    ts, xs = gillespie_bd(g, k, tmax, int(round(x_bar)), rng)
    ax.step(ts, xs, where="post", label=f"g={g:g}, x_bar={x_bar:g}")
    # time-weighted stats (burn-in small since we start at steady state)
    m, s = time_weighted_stats(ts, xs, tmax, burn_frac=0.1)
    means_ss.append(m)
    stds_ss.append(s)
    print(f"[start=steady state] g={g:g}  x_bar={x_bar:g}  "
          f"time-weighted mean={m:.3f}  std={s:.3f}  "
          f"sqrt(x_bar)={np.sqrt(x_bar):.3f}  rel_noise(std/mean)={s/m:.4f}")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Trajectories started at steady state")
ax.set_yscale("log")
ax.legend()

# ---- Trajectories started from zero (top-right) ----
ax = axes[0, 1]
for g in g_list:
    x_bar = g / k
    rng = np.random.default_rng(seed)
    ts, xs = gillespie_bd(g, k, tmax, 0, rng)
    ax.step(ts, xs, where="post", label=f"g={g:g} (from 0)")
    # deterministic relaxation curve x_bar*(1 - exp(-k t)) for comparison
    tt = np.linspace(0, tmax, 400)
    ax.plot(tt, x_bar * (1 - np.exp(-k * tt)), "k--", lw=1)
    # report how close the trajectory end is to the analytic climb
    x_end_theory = x_bar * (1 - np.exp(-k * tmax))
    print(f"[start=zero]        g={g:g}  x_bar={x_bar:g}  "
          f"analytic x(tmax)=x_bar*(1-exp(-k*tmax))={x_end_theory:.3f}  "
          f"sim x at last event={xs[-1]:d}")
ax.set_xlabel("time")
ax.set_ylabel("molecule count x")
ax.set_title("Started from zero (dashed = x_bar*(1-exp(-k t)))")
ax.legend()

# ---- std vs mean, compared to sqrt(x_bar) (bottom-left) ----
means_ss = np.array(means_ss)
stds_ss = np.array(stds_ss)
xbars = np.array(xbars)
ax = axes[1, 0]
ax.loglog(means_ss, stds_ss, "o", ms=9, label="measured std vs mean")
xline = np.logspace(np.log10(xbars.min()) - 0.2, np.log10(xbars.max()) + 0.2, 50)
ax.loglog(xline, np.sqrt(xline), "-", label="sqrt(x_bar) (Poisson)")
ax.set_xlabel("mean x")
ax.set_ylabel("std of x")
ax.set_title("Noise check: std vs mean tracks sqrt(x_bar)")
ax.legend()

# ---- relative noise vs 1/sqrt(x_bar) (bottom-right) ----
ax = axes[1, 1]
rel_measured = stds_ss / means_ss
ax.loglog(xbars, rel_measured, "o", ms=9, label="measured std/mean")
ax.loglog(xline, 1.0 / np.sqrt(xline), "-", label="1/sqrt(x_bar)")
ax.set_xlabel("x_bar")
ax.set_ylabel("relative noise (std/mean)")
ax.set_title("Relative noise falls as 1/sqrt(x_bar)")
ax.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/8D.2.1_s1.png")

# ---- summary of the Poisson check ----
print("\n--- Poisson check summary (steady-state runs) ---")
for xb, m, s in zip(xbars, means_ss, stds_ss):
    print(f"x_bar={xb:g}: std={s:.3f}, sqrt(x_bar)={np.sqrt(xb):.3f}, "
          f"ratio std/sqrt(x_bar)={s/np.sqrt(xb):.3f}, "
          f"rel_noise={s/m:.4f}, 1/sqrt(x_bar)={1/np.sqrt(xb):.4f}")

# One-sentence explanation:
print("\nExplanation: Because the birth-death steady state is Poisson (variance = mean = "
      "x_bar), the standard deviation equals sqrt(x_bar), so relative noise std/mean = "
      "1/sqrt(x_bar) shrinks as x_bar grows; matching this scaling confirms the simulation "
      "reproduces the correct intrinsic noise.")
