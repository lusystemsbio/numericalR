import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters (toggle switch with mutual repression) ---
gX0 = gY0 = 10.0     # basal production rate
gX1 = gY1 = 40.0     # max regulated production rate
X0 = Y0 = 100.0      # repression thresholds
nX = nY = 4.0        # Hill coefficients
kX = kY = 0.1        # degradation rates
b = 20.0             # additive noise strength

# --- Simulation controls ---
T = 1000.0
dt = 0.01
N = int(T / dt)
sqrt_dt = np.sqrt(dt)     # dW ~ N(0, dt) => sqrt(dt) * N(0,1)
np.random.seed(3)

# --- Storage ---
t = np.linspace(0.0, T, N + 1)
X = np.empty(N + 1)
Y = np.empty(N + 1)
X[0], Y[0] = 50.0, 200.0   # initial condition (low-X / high-Y state)

# --- Drift functions (deterministic part of each SDE) ---
def driftX(x, y):
    return gX0 + gX1 / (1.0 + (y / Y0) ** nY) - kX * x

def driftY(x, y):
    return gY0 + gY1 / (1.0 + (x / X0) ** nX) - kY * y

# --- Euler-Maruyama integration for the two-variable SDE ---
for i in range(N):
    x, y = X[i], Y[i]
    dWx = sqrt_dt * np.random.randn()   # Wiener increment for X
    dWy = sqrt_dt * np.random.randn()   # Wiener increment for Y
    xn = x + driftX(x, y) * dt + b * dWx   # deterministic step + noise kick
    yn = y + driftY(x, y) * dt + b * dWy
    X[i + 1] = max(xn, 0.0)   # keep concentrations non-negative
    Y[i + 1] = max(yn, 0.0)

# --- Time-series and phase-plane plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(t, X, lw=0.6, label="X")
ax1.plot(t, Y, lw=0.6, label="Y")
ax1.set_xlabel("t"); ax1.set_ylabel("concentration")
ax1.set_title("Toggle switch time series (hops between states)")
ax1.legend()
ax2.plot(X, Y, lw=0.3, color="purple")
ax2.set_xlabel("X"); ax2.set_ylabel("Y")
ax2.set_title("Phase plane")
fig.tight_layout()
fig.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.1.1_s2.png", dpi=120)

# --- Bistability / transition check ---
# Classify each timepoint by which gene wins (X high vs Y high).
X_high = X > Y            # boolean mask for high-X/low-Y state
Y_high = ~X_high          # high-Y/low-X state
frac_X_high = X_high.mean()
frac_Y_high = Y_high.mean()

# Count hops = sign changes in (X - Y), i.e. state swaps along the trajectory.
sign = np.sign(X - Y)
sign[sign == 0] = 1
n_transitions = int(np.sum(np.abs(np.diff(sign)) > 0))

# Mean levels conditioned on each state (should show two well-separated modes).
meanX_in_Xhigh = X[X_high].mean() if X_high.any() else float("nan")
meanY_in_Xhigh = Y[X_high].mean() if X_high.any() else float("nan")
meanX_in_Yhigh = X[Y_high].mean() if Y_high.any() else float("nan")
meanY_in_Yhigh = Y[Y_high].mean() if Y_high.any() else float("nan")

print("Fraction of time in high-X/low-Y state:", frac_X_high)
print("Fraction of time in high-Y/low-X state:", frac_Y_high)
print("Number of state transitions (hops):", n_transitions)
print("Mean X in high-X state:", meanX_in_Xhigh)
print("Mean Y in high-X state:", meanY_in_Xhigh)
print("Mean X in high-Y state:", meanX_in_Yhigh)
print("Mean Y in high-Y state:", meanY_in_Yhigh)
print("Final (X, Y):", X[-1], Y[-1])
# One-sentence explanation:
print("Explanation: Finding that the trajectory spends substantial time in BOTH "
      "a high-X/low-Y state and a high-Y/low-X state while repeatedly switching "
      "between them confirms bistability with noise-driven transitions, because a "
      "monostable system would settle into a single state and never hop.")
