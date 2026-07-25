import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters (toggle switch with mutual repression) ----
gX0 = gY0 = 10.0      # basal production
gX1 = gY1 = 40.0      # max regulated production
X0 = Y0 = 100.0       # repression thresholds
nX = nY = 4.0         # Hill coefficients
kX = kY = 0.1         # degradation rates
b = 20.0              # additive noise strength

# ---- Simulation settings ----
T = 1000.0
dt = 0.01
N = int(T / dt)
np.random.seed(3)

# Deterministic drift terms (repressive Hill functions)
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Y0) ** nY) - kX * X

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / X0) ** nX) - kY * Y

# ---- Euler-Maruyama integration (explicit) ----
X = np.empty(N + 1)
Y = np.empty(N + 1)
X[0], Y[0] = 50.0, 200.0
sqrt_dt = np.sqrt(dt)

for i in range(N):
    # drift * dt
    dX_det = fX(X[i], Y[i]) * dt
    dY_det = fY(X[i], Y[i]) * dt
    # diffusion: b * dW, with dW ~ N(0, dt) = sqrt(dt) * N(0,1)
    dWX = sqrt_dt * np.random.randn()
    dWY = sqrt_dt * np.random.randn()
    # update
    Xn = X[i] + dX_det + b * dWX
    Yn = Y[i] + dY_det + b * dWY
    # keep concentrations non-negative
    X[i + 1] = max(Xn, 0.0)
    Y[i + 1] = max(Yn, 0.0)

t = np.linspace(0.0, T, N + 1)

# ---- Report basic trajectory statistics ----
print("Final X:", X[-1])
print("Final Y:", Y[-1])
print("Mean X:", np.mean(X))
print("Mean Y:", np.mean(Y))

# ---- Bistability / hopping check ----
# Deterministic fixed points near each pure state: gene i high => X_high = (g0+g1)/k, X_low = g0/k
X_high = (gX0 + gX1) / kX
X_low = gX0 / kX
print("Predicted high state level:", X_high)
print("Predicted low state level:", X_low)

# Classify each time point by which gene dominates (X - Y sign)
state = np.sign(X - Y)                 # +1 : high-X/low-Y ; -1 : low-X/high-Y
frac_highX = np.mean(state > 0)
frac_highY = np.mean(state < 0)
print("Fraction of time in high-X/low-Y state:", frac_highX)
print("Fraction of time in low-X/high-Y state:", frac_highY)

# Count transitions (sign changes of X - Y) => hops between states
transitions = np.sum(np.abs(np.diff(np.sign(X - Y))) > 0)
print("Number of state transitions (hops):", int(transitions))

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

ax1.plot(t, X, label="X", color="tab:blue", lw=0.7)
ax1.plot(t, Y, label="Y", color="tab:red", lw=0.7)
ax1.set_xlabel("time")
ax1.set_ylabel("concentration")
ax1.set_title("Toggle switch time series (hops between states)")
ax1.legend()

ax2.plot(X, Y, color="tab:purple", lw=0.4, alpha=0.6)
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_title("Phase plane (two clusters = two states)")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.1.1_s5.png")

# One-sentence explanation of the check:
print("Explanation: Finding that the trajectory spends substantial time in BOTH "
      "a high-X/low-Y and a low-X/high-Y state and repeatedly transitions between "
      "them confirms bistability with noise-driven hops, since a monostable system "
      "would settle into a single state and never swap which gene is highly expressed.")
