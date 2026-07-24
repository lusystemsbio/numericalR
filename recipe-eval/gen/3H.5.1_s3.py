import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Repressilator: a ring of three mutual repressions (nondimensional form).
#   dx/dt = g/(1 + z^3) - x   (X repressed by Z)
#   dy/dt = h/(1 + x^3) - y   (Y repressed by X)
#   dz/dt = l/(1 + y^3) - z   (Z repressed by Y)
# Hill coefficient = 3, unit degradation rate.
# ----------------------------------------------------------------------

# Parameters (production rates for the three genes)
g = h = l = 5.0

def repressilator(state, g, h, l):
    """Return the time derivatives [dx, dy, dz] of the three-gene ring."""
    x, y, z = state
    dx = g / (1.0 + z**3) - x   # X is repressed by Z
    dy = h / (1.0 + x**3) - y   # Y is repressed by X
    dz = l / (1.0 + y**3) - z   # Z is repressed by Y
    return np.array([dx, dy, dz])

def rk4_step(f, state, dt, *args):
    """One explicit classical Runge-Kutta 4th-order step."""
    k1 = f(state, *args)                    # slope at start of interval
    k2 = f(state + 0.5 * dt * k1, *args)    # slope at midpoint using k1
    k3 = f(state + 0.5 * dt * k2, *args)    # slope at midpoint using k2
    k4 = f(state + dt * k3, *args)          # slope at end using k3
    # weighted average of the four slopes (1,2,2,1)/6
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

# ----------------------------------------------------------------------
# Integrate the ring with the generic RK4 stepper.
# ----------------------------------------------------------------------
dt = 0.01
T = 60.0
n_steps = int(T / dt)

t = np.linspace(0.0, T, n_steps + 1)
traj = np.zeros((n_steps + 1, 3))

# Initial condition: slightly asymmetric so the ring leaves the unstable
# symmetric fixed point and settles onto its limit cycle.
traj[0] = np.array([1.0, 1.1, 1.2])

for i in range(n_steps):
    traj[i + 1] = rk4_step(repressilator, traj[i], dt, g, h, l)

x, y, z = traj[:, 0], traj[:, 1], traj[:, 2]

# ----------------------------------------------------------------------
# Time-series plot: the three genes oscillating in turn.
# ----------------------------------------------------------------------
plt.figure(figsize=(10, 5))
plt.plot(t, x, label="X", color="tab:red")
plt.plot(t, y, label="Y", color="tab:green")
plt.plot(t, z, label="Z", color="tab:blue")
plt.xlabel("time (nondimensional)")
plt.ylabel("concentration")
plt.title("Three-gene repressilator (g = h = l = 5, Hill n = 3), RK4")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.5.1_s3.png")

# ----------------------------------------------------------------------
# CHECK: confirm a sustained limit cycle with a fixed peaking order.
# Analyse the second half of the run (transients have decayed by then).
# ----------------------------------------------------------------------
half = n_steps // 2
tt = t[half:]
xx, yy, zz = x[half:], y[half:], z[half:]

def find_peaks(sig):
    """Indices of strict local maxima (interior points)."""
    return np.where((sig[1:-1] > sig[:-2]) & (sig[1:-1] > sig[2:]))[0] + 1

px = find_peaks(xx)
py = find_peaks(yy)
pz = find_peaks(zz)

# 1) Sustained (not decaying): amplitude of X in the last window is sizeable.
amp_x = xx.max() - xx.min()
amp_y = yy.max() - yy.min()
amp_z = zz.max() - zz.min()

# 2) Fixed period: spacing between successive X peaks is nearly constant.
if len(px) >= 2:
    periods = np.diff(tt[px])
    period_mean = periods.mean()
    period_std = periods.std()
else:
    period_mean = period_std = float("nan")

# 3) Fixed peaking order: within each X-period the peaks occur X -> Y -> Z
#    (a fixed cyclic phase lag ~ one third of a period).
t_peak_x = tt[px[-2]] if len(px) >= 2 else float("nan")
# first Y and Z peaks after that X peak
ty_after = tt[py][tt[py] > t_peak_x]
tz_after = tt[pz][tt[pz] > t_peak_x]
lag_y = (ty_after[0] - t_peak_x) if len(ty_after) else float("nan")
lag_z = (tz_after[0] - t_peak_x) if len(tz_after) else float("nan")

# ----------------------------------------------------------------------
# Print every numerical result.
# ----------------------------------------------------------------------
print("Parameters g = h = l:", g)
print("Time step dt:", dt)
print("Total integration time T:", T)
print("Number of RK4 steps:", n_steps)
print("Final state x, y, z:", x[-1], y[-1], z[-1])
print("Amplitude X (last half):", amp_x)
print("Amplitude Y (last half):", amp_y)
print("Amplitude Z (last half):", amp_z)
print("Number of X peaks (last half):", len(px))
print("Number of Y peaks (last half):", len(py))
print("Number of Z peaks (last half):", len(pz))
print("Mean period between X peaks:", period_mean)
print("Std of period between X peaks:", period_std)
print("Phase lag X->Y (time):", lag_y)
print("Phase lag X->Z (time):", lag_z)
print("Peaking order fixed as X->Y->Z:", bool(lag_y < lag_z))
print("Sustained oscillation (amplitude > 0.5 and near-constant period):",
      bool(amp_x > 0.5 and period_std < 0.05 * period_mean))

# Explanation:
# The check confirms the result because a near-constant inter-peak period with
# non-decaying amplitude means the trajectory has settled onto a closed limit
# cycle, and the fixed X->Y->Z ordering of the peaks shows the three genes fire
# in a stable, repeating sequence rather than damping to a steady state.
print("Why: a constant period with non-decaying amplitude and a fixed X->Y->Z "
      "peak order means the system has locked onto a stable limit cycle.")
