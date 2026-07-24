import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Stiff model parameters -------------------------------------------------
# Constitutively expressed gene: dX/dt = g - k*X
# Large degradation rate k makes the problem stiff.
g = 50.0      # constitutive production rate
k = 10.0      # degradation rate (large -> stiff)
X0 = 3.0      # initial amount
t0, tf = 0.0, 4.0
dt = 0.2

# Steady state (exact solution decays toward this value)
Xss = g / k
print(f"Steady state g/k = {Xss}")

# Time grid
t = np.arange(t0, tf + dt/2, dt)
n = len(t)

# --- Exact solution ---------------------------------------------------------
# X(t) = g/k + (X0 - g/k) * exp(-k t)
X_exact = Xss + (X0 - Xss) * np.exp(-k * t)

# --- Forward (explicit) Euler ----------------------------------------------
# X_{n+1} = X_n + (g - k*X_n)*dt
# Stability requires dt < 2/k = 0.2; at dt = 0.2 it sits on the edge/unstable.
X_fe = np.empty(n)
X_fe[0] = X0
for i in range(n - 1):
    rhs = g - k * X_fe[i]          # evaluate RHS at the CURRENT (known) point
    X_fe[i + 1] = X_fe[i] + rhs * dt

# --- Backward (implicit) Euler ---------------------------------------------
# X_{n+1} = X_n + (g - k*X_{n+1})*dt
# Because the RHS is linear, solve algebraically for X_{n+1}:
#   X_{n+1}(1 + dt*k) = X_n + dt*g
#   X_{n+1} = (X_n + dt*g) / (1 + dt*k)
X_be = np.empty(n)
X_be[0] = X0
for i in range(n - 1):
    X_be[i + 1] = (X_be[i] + dt * g) / (1.0 + dt * k)  # closed-form implicit update

# --- Print results ----------------------------------------------------------
print(f"Amplification factor forward Euler |1 - dt*k| = {abs(1 - dt*k)}")
print(f"Amplification factor backward Euler 1/(1 + dt*k) = {1.0/(1.0 + dt*k)}")
print(f"Final time t = {t[-1]}")
print(f"Forward Euler  X(final) = {X_fe[-1]}")
print(f"Backward Euler X(final) = {X_be[-1]}")
print(f"Exact          X(final) = {X_exact[-1]}")
print(f"Forward Euler  max |X| over run = {np.max(np.abs(X_fe))}")
print(f"Backward Euler max |X| over run = {np.max(np.abs(X_be))}")
print(f"Forward Euler diverges (final magnitude > 1e3)? {abs(X_fe[-1]) > 1e3}")
print(f"Backward Euler stable (final within 0.1 of g/k)? {abs(X_be[-1] - Xss) < 0.1}")

# --- Plot -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t, X_exact, 'k-', lw=2, label='Exact')
ax.plot(t, X_be, 'bo-', label='Backward (implicit) Euler')
ax.plot(t, X_fe, 'rs--', label='Forward (explicit) Euler')
ax.axhline(Xss, color='gray', ls=':', label='g/k = 5')
ax.set_xlabel('t')
ax.set_ylabel('X(t)')
ax.set_title(f'Stiff gene expression: g={g}, k={k}, dt={dt}')
ax.legend()
# Clip y-axis so the stable solutions remain visible despite FE blow-up
ax.set_ylim(-5, 20)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.8.1_s1.png")

# The check confirms the result because at dt=0.2 the forward-Euler amplification
# factor |1 - dt*k| = 1 is not below 1 (so errors fail to decay and the numerical
# solution oscillates/diverges), while the backward-Euler factor 1/(1 + dt*k) < 1
# is unconditionally < 1, so its solution decays smoothly to g/k regardless of dt.
print("Check: forward Euler oscillates/diverges (|1-dt*k|>=1) while backward Euler stays stable (1/(1+dt*k)<1) and tracks the exact decay to g/k.")
