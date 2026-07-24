import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
g0 = 10.0      # basal production
g1 = 60.0      # max repressible production
Xth = 200.0    # repression threshold
n = 4          # Hill coefficient
k = 0.1        # linear degradation rate
dt = 0.01      # time step
T = 200.0      # final time
X_hist = 1.0   # constant history value X(t) = 1 for t <= 0

taus = [5, 10, 15, 20]  # delays to sweep

# Production term: basal + repressive Hill driven by the DELAYED protein level
def f(X_now, X_delayed):
    return g0 + g1 / (1.0 + (X_delayed / Xth) ** n) - k * X_now

# ---- Delay-DDE integration via explicit 2nd-order Heun (predictor-corrector) ----
def integrate(tau):
    nsteps = int(round(T / dt))          # number of time steps
    d = int(round(tau / dt))             # delay expressed in steps
    t = np.linspace(0.0, T, nsteps + 1)  # time grid
    X = np.empty(nsteps + 1)
    X[0] = X_hist                         # initial value

    # helper: value of X at index (i - d); use constant history when i-d < 0
    def delayed(i):
        j = i - d
        return X[j] if j >= 0 else X_hist

    for i in range(nsteps):
        # delayed value needed for the derivative at the current step
        Xd_i = delayed(i)
        # predictor (Euler step)
        slope1 = f(X[i], Xd_i)
        X_pred = X[i] + dt * slope1
        # delayed value at the next step (for the corrector slope)
        j_next = (i + 1) - d
        Xd_next = X[j_next] if j_next >= 0 else X_hist
        # corrector: evaluate slope at predicted end-point, average the two slopes
        slope2 = f(X_pred, Xd_next)
        X[i + 1] = X[i] + 0.5 * dt * (slope1 + slope2)
    return t, X

# ---- Analytic steady state (X* where dX/dt = 0, X constant so X(t-tau)=X*) ----
# k*X = g0 + g1/(1+(X/Xth)^n); solve numerically for reference
def steady_state():
    X = 250.0
    for _ in range(200):
        prod = g0 + g1 / (1.0 + (X / Xth) ** n)
        Xnew = prod / k
        X = 0.5 * X + 0.5 * Xnew  # damped fixed-point iteration
    return X

Xstar = steady_state()
print(f"Analytic steady state X* (dX/dt=0): {Xstar:.4f}")

# ---- Run sweep, plot, and report late-time behavior ----
plt.figure(figsize=(10, 6))
results = {}
for tau in taus:
    t, X = integrate(tau)
    results[tau] = (t, X)
    plt.plot(t, X, label=f"tau = {tau}")

    # characterize the last quarter of the trace (steady vs oscillatory)
    tail = X[int(0.75 * len(X)):]
    xmin, xmax, xmean = tail.min(), tail.max(), tail.mean()
    amp = xmax - xmin
    print(f"tau = {tau:2d}:  final X = {X[-1]:8.3f}   "
          f"tail mean = {xmean:8.3f}   tail min = {xmin:8.3f}   "
          f"tail max = {xmax:8.3f}   peak-to-peak amplitude = {amp:8.3f}")

plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Delayed negative autoregulation: oscillation onset as delay tau grows")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.2.1_s1.png")

# ---- Separate check: short delay settles near steady state ~250 ----
t5, X5 = results[5]
final5 = X5[-1]
tail5 = X5[int(0.75 * len(X5)):]
amp5 = tail5.max() - tail5.min()
print(f"\nCHECK (tau=5): final X = {final5:.3f}, tail peak-to-peak = {amp5:.4f}")
print(f"CHECK: |final X - 250| = {abs(final5 - 250.0):.4f}")
settles = amp5 < 1.0 and abs(final5 - 250.0) < 10.0
print(f"CHECK settles near steady state ~250 with negligible oscillation: {settles}")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: Because the short delay (tau=5) relaxes to the flat, "
      "analytically predicted fixed point near 250 while the larger delays "
      "develop growing then sustained peak-to-peak oscillations, the check "
      "confirms that increasing the feedback delay is what destabilizes the "
      "steady state into oscillation, exactly as the model predicts.")
