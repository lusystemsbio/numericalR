import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---- Generic multi-variable RK4 from Part 3A ----
# Integrates dy/dt = f(t, y) where y is a list/vector of state variables.
def rk4_step(f, t, y, h):
    # each k is a vector (one entry per state variable)
    k1 = f(t, y)
    k2 = f(t + h / 2, [yi + h / 2 * k1i for yi, k1i in zip(y, k1)])
    k3 = f(t + h / 2, [yi + h / 2 * k2i for yi, k2i in zip(y, k2)])
    k4 = f(t + h, [yi + h * k3i for yi, k3i in zip(y, k3)])
    # weighted average of the four slope estimates
    return [yi + h / 6 * (a + 2 * b + 2 * c + d)
            for yi, a, b, c, d in zip(y, k1, k2, k3, k4)]


def integrate(f, y0, t0, t_end, h):
    # explicit time stepping loop (not a one-call routine)
    t = t0
    y = list(y0)
    ts = [t]
    ys = [y]
    n_steps = int(round((t_end - t0) / h))
    for _ in range(n_steps):
        y = rk4_step(f, t, y, h)
        t += h
        ts.append(t)
        ys.append(y)
    return ts, ys


# ---- Chemostat model ----
# y = [N, C]
# dN/dt = a1*(C/(C+1))*N - N        growth by Michaelis-Menten uptake minus washout
# dC/dt = -(C/(C+1))*N - C + a2     consumption minus washout plus scaled feed
def make_chemostat(a1, a2):
    def f(t, y):
        N, C = y
        mm = C / (C + 1.0)          # Michaelis-Menten saturation term
        dN = a1 * mm * N - N
        dC = -mm * N - C + a2
        return [dN, dC]
    return f


# ---- Parameters ----
a1, a2 = 2.0, 5.0
f = make_chemostat(a1, a2)

t0, t_end, h = 0.0, 40.0, 0.01

# Run 1: N(0) = 0 (no seed) -> should reach washout (0, 5)
ts1, ys1 = integrate(f, [0.0, 5.0], t0, t_end, h)
N1_final, C1_final = ys1[-1]

# Run 2: tiny seed N(0) = 0.01 -> should reach coexistence (8, 1)
ts2, ys2 = integrate(f, [0.01, 5.0], t0, t_end, h)
N2_final, C2_final = ys2[-1]

# ---- Report numerical results ----
print("Parameters: a1 = {:.1f}, a2 = {:.1f}".format(a1, a2))
print("Washout equilibrium (analytic):        N = 0.0, C = 5.0")
print("Coexistence equilibrium (analytic):     N = 8.0, C = 1.0")
print("Run 1 start N(0) = 0.0  -> final N = {:.6f}".format(N1_final))
print("Run 1 start N(0) = 0.0  -> final C = {:.6f}".format(C1_final))
print("Run 2 start N(0) = 0.01 -> final N = {:.6f}".format(N2_final))
print("Run 2 start N(0) = 0.01 -> final C = {:.6f}".format(C2_final))

# Minimum N along run 2 shows the initial drift toward washout before peeling away
N2_series = [y[0] for y in ys2]
min_N2 = min(N2_series)
min_idx = N2_series.index(min_N2)
print("Run 2 minimum N along path = {:.6f} at t = {:.3f} (drift toward washout)".format(
    min_N2, ts2[min_idx]))
print("Run 2 then peels away to coexistence, final N = {:.6f}".format(N2_final))

# ---- Phase-plane plot ----
N1 = [y[0] for y in ys1]
C1 = [y[1] for y in ys1]
N2 = [y[0] for y in ys2]
C2 = [y[1] for y in ys2]

plt.figure(figsize=(8, 6))
plt.plot(N1, C1, 'b-', lw=2, label='N(0)=0 -> washout')
plt.plot(N2, C2, 'r-', lw=2, label='N(0)=0.01 -> coexistence')
plt.plot(0, 5, 'ks', ms=10, label='washout (0, 5)')
plt.plot(8, 1, 'k^', ms=10, label='coexistence (8, 1)')
plt.plot(N1[0], C1[0], 'bo', ms=6)
plt.plot(N2[0], C2[0], 'ro', ms=6)
plt.xlabel('Population N')
plt.ylabel('Substrate C')
plt.title('Chemostat phase plane: washout vs coexistence')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3C.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because N=0 makes dN/dt=0 exactly, no population can ever "
      "appear and the system must relax to washout (0,5), whereas the tiny seed "
      "N=0.01 is enough to let growth take over once it nears the nutrient-rich "
      "washout point, so it peels away to the stable coexistence state (8,1) -- "
      "confirming washout is only reachable with literally zero population.")
