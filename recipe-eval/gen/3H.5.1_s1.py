import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Right-hand side of the repressilator ODEs (three-gene repression ring) ---
# X repressed by Z, Y repressed by X, Z repressed by Y; Hill coefficient 3, unit degradation.
def rhs(state, g, h, l):
    x, y, z = state
    dx = g / (1.0 + z**3) - x   # X repressed by Z
    dy = h / (1.0 + x**3) - y   # Y repressed by X
    dz = l / (1.0 + y**3) - z   # Z repressed by Y
    return np.array([dx, dy, dz])

# --- Generic explicit RK4 integrator (implemented step by step, not a library call) ---
def rk4(rhs, state0, t0, t1, dt, *params):
    n = int(round((t1 - t0) / dt))          # number of steps
    ts = np.empty(n + 1)                     # storage for time
    ys = np.empty((n + 1, len(state0)))      # storage for state
    ts[0] = t0
    ys[0] = state0
    y = np.array(state0, dtype=float)
    t = t0
    for i in range(n):
        k1 = rhs(y, *params)                 # slope at start
        k2 = rhs(y + 0.5 * dt * k1, *params) # slope at midpoint using k1
        k3 = rhs(y + 0.5 * dt * k2, *params) # slope at midpoint using k2
        k4 = rhs(y + dt * k3, *params)       # slope at end using k3
        y = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)  # weighted average
        t = t + dt
        ts[i + 1] = t
        ys[i + 1] = y
    return ts, ys

# --- Parameters and integration ---
g = h = l = 5.0
state0 = [1.0, 1.5, 2.0]   # asymmetric start so the ring is off-equilibrium
dt = 0.01
t0, t1 = 0.0, 60.0
ts, ys = rk4(rhs, state0, t0, t1, dt, g, h, l)
x, y, z = ys[:, 0], ys[:, 1], ys[:, 2]

# --- Time-series plot of the three genes oscillating in turn ---
plt.figure(figsize=(10, 5))
plt.plot(ts, x, label="x (repressed by z)")
plt.plot(ts, y, label="y (repressed by x)")
plt.plot(ts, z, label="z (repressed by y)")
plt.xlabel("time (nondimensional)")
plt.ylabel("protein concentration")
plt.title("Repressilator: three-gene repression ring (g=h=l=5, Hill n=3)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.5.1_s1.png")

# --- Check: sustained limit cycle with a fixed peak order ---
# Find interior local maxima (peaks) of each species in the latter half of the run,
# where transients have died out, by simple three-point comparison.
def find_peaks(signal, tarr, tmin):
    peak_times = []
    peak_vals = []
    for i in range(1, len(signal) - 1):
        if tarr[i] >= tmin and signal[i] > signal[i-1] and signal[i] > signal[i+1]:
            peak_times.append(tarr[i])
            peak_vals.append(signal[i])
    return np.array(peak_times), np.array(peak_vals)

t_settle = 30.0  # ignore early transient; look only at the settled regime
xt, xv = find_peaks(x, ts, t_settle)
yt, yv = find_peaks(y, ts, t_settle)
zt, zv = find_peaks(z, ts, t_settle)

# Sustained oscillation: amplitude of the last few peaks stays roughly constant (not decaying).
def amp_stats(vals, tag):
    if len(vals) >= 2:
        recent = vals[-min(4, len(vals)):]
        print(f"{tag}: number of settled peaks = {len(vals)}")
        print(f"{tag}: last peak heights = {np.round(recent, 4).tolist()}")
        print(f"{tag}: peak-height spread (max-min of last peaks) = {recent.max() - recent.min():.6f}")
    else:
        print(f"{tag}: too few peaks detected")

amp_stats(xv, "x")
amp_stats(yv, "y")
amp_stats(zv, "z")

# Period from spacing between consecutive x-peaks (should be nearly constant for a limit cycle).
if len(xt) >= 2:
    periods = np.diff(xt)
    print(f"x inter-peak periods = {np.round(periods, 4).tolist()}")
    print(f"mean period = {periods.mean():.6f}")
    print(f"period standard deviation = {periods.std():.6f}")

# Peak order: within one representative cycle, x should peak, then y, then z (fixed cyclic order).
if len(xt) >= 1 and len(yt) >= 1 and len(zt) >= 1:
    tx = xt[0]                                   # reference: first settled x-peak
    ty = yt[yt > tx][0] if np.any(yt > tx) else np.nan   # next y-peak after it
    tz = zt[zt > tx][0] if np.any(zt > tx) else np.nan   # next z-peak after it
    print(f"first settled x-peak time = {tx:.4f}")
    print(f"next y-peak time = {ty:.4f}")
    print(f"next z-peak time = {tz:.4f}")
    ordered = tx < ty < tz
    print(f"peaks occur in fixed order x -> y -> z: {ordered}")

# Explanation of why this check confirms the result:
print("Explanation: constant peak heights and a constant inter-peak period (rather than "
      "decaying or growing) show the trajectory has settled onto a self-sustaining closed "
      "limit cycle, and the fixed x->y->z peak ordering confirms the phase-shifted "
      "repression ring rather than a decay to steady state.")
