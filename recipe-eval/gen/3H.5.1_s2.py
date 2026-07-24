import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Repressilator: ring of three mutual repressions (nondimensional)
#   dx/dt = g/(1 + z^3) - x     (X repressed by Z)
#   dy/dt = h/(1 + x^3) - y     (Y repressed by X)
#   dz/dt = l/(1 + y^3) - z     (Z repressed by Y)
# Hill coefficient n = 3, unit degradation rate.
# ---------------------------------------------------------------

# Parameters (equal maximal production rates)
g = h = l = 5.0

def deriv(state):
    # unpack the three protein concentrations
    x, y, z = state
    dx = g / (1.0 + z**3) - x   # X production repressed by Z, minus decay
    dy = h / (1.0 + x**3) - y   # Y production repressed by X, minus decay
    dz = l / (1.0 + y**3) - z   # Z production repressed by Y, minus decay
    return np.array([dx, dy, dz])

# ---------------------------------------------------------------
# Generic classical Runge-Kutta 4th order (RK4), written out explicitly
# ---------------------------------------------------------------
def rk4_step(state, dt):
    k1 = deriv(state)                    # slope at start of interval
    k2 = deriv(state + 0.5 * dt * k1)    # slope at midpoint using k1
    k3 = deriv(state + 0.5 * dt * k2)    # slope at midpoint using k2
    k4 = deriv(state + dt * k3)          # slope at end using k3
    # weighted average of the four slopes advances the state
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

# ---------------------------------------------------------------
# Integrate the system
# ---------------------------------------------------------------
dt = 0.01                # time step
T = 60.0                 # total integration time
nsteps = int(T / dt)     # number of RK4 steps

t = np.zeros(nsteps + 1)
traj = np.zeros((nsteps + 1, 3))

# Asymmetric initial condition so the ring is off the unstable steady state
state = np.array([1.0, 1.5, 3.0])
traj[0] = state

for i in range(nsteps):
    state = rk4_step(state, dt)   # advance one RK4 step
    traj[i + 1] = state
    t[i + 1] = t[i] + dt

x, y, z = traj[:, 0], traj[:, 1], traj[:, 2]

# ---------------------------------------------------------------
# Check: confirm a sustained limit cycle with a fixed peak order.
# Compare successive periods late in the run: (a) peak amplitudes are
# steady (not decaying to a fixed point), and (b) X, Y, Z peak in a
# repeating fixed cyclic order.
# ---------------------------------------------------------------

def find_peaks(signal):
    # local maxima: interior points larger than both neighbors
    idx = np.where((signal[1:-1] > signal[:-2]) & (signal[1:-1] > signal[2:]))[0] + 1
    return idx

# use the second half of the trajectory (after transients die out)
half = nsteps // 2
px = find_peaks(x); px = px[px >= half]
py = find_peaks(y); py = py[py >= half]
pz = find_peaks(z); pz = pz[pz >= half]

# amplitude stability: spread of peak heights in the steady portion
x_peak_amps = x[px]
amp_mean = np.mean(x_peak_amps)
amp_spread = np.max(x_peak_amps) - np.min(x_peak_amps)

# period from spacing between consecutive X peaks
periods = np.diff(t[px])
period_mean = np.mean(periods)
period_spread = np.max(periods) - np.min(periods)

# peak order: merge all peaks by time, label by gene, read the cyclic pattern
events = sorted(
    [(t[i], 'X') for i in px] +
    [(t[i], 'Y') for i in py] +
    [(t[i], 'Z') for i in pz]
)
order = [g_ for _, g_ in events]

# check that consecutive triples repeat a single fixed cyclic ordering
first_triple = order[:3]
is_fixed_order = all(
    order[k:k+3] == first_triple
    for k in range(0, len(order) - 2, 3)
)

print("Parameters g = h = l:", g)
print("Time step dt:", dt)
print("Total integration time T:", T)
print("Final state (x, y, z):", x[-1], y[-1], z[-1])
print("Number of X peaks in steady portion:", len(px))
print("Mean X peak amplitude:", amp_mean)
print("X peak amplitude spread (max-min):", amp_spread)
print("Mean oscillation period:", period_mean)
print("Period spread (max-min):", period_spread)
print("Observed peak order sequence:", order)
print("Repeating peak-order triple:", first_triple)
print("Fixed cyclic peak order confirmed:", is_fixed_order)
sustained = (amp_spread < 1e-2) and (period_spread < 1e-2)
print("Sustained limit cycle confirmed:", sustained)

# ---------------------------------------------------------------
# Time-series plot
# ---------------------------------------------------------------
plt.figure(figsize=(10, 5))
plt.plot(t, x, label="x (Gene X)")
plt.plot(t, y, label="y (Gene Y)")
plt.plot(t, z, label="z (Gene Z)")
plt.xlabel("time")
plt.ylabel("concentration")
plt.title("Repressilator: three genes oscillating in turn (g=h=l=5, RK4)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.5.1_s2.png")

# One-sentence explanation:
# The check confirms the result because constant peak amplitudes and periods
# across late cycles rule out decay to a steady state, while the unchanging
# cyclic X->Z->Y peak order shows a genuine, phase-locked sustained limit cycle.
print("Why the check works: steady peak amplitudes/periods rule out settling to a fixed point, and the unchanging cyclic peak order proves a phase-locked sustained limit cycle.")
