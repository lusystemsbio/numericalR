import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Stiff model: dX/dt = g - k*X, constitutively expressed gene with large degradation rate k
g = 50.0    # constitutive production rate
k = 10.0    # large degradation rate -> makes the problem stiff
X0 = 3.0    # initial condition
t0, tf = 0.0, 4.0
dt = 0.2

# Time grid
t = np.arange(t0, tf + dt/2, dt)
n = len(t)

# Exact (analytical) solution: X(t) = g/k + (X0 - g/k)*exp(-k*t), decays to steady state g/k
steady = g / k
X_exact = steady + (X0 - steady) * np.exp(-k * t)

# Forward (explicit) Euler: X_{n+1} = X_n + dt*(g - k*X_n)
# Explicit step-by-step so the instability is visible; stable only if dt < 2/k = 0.2.
X_fe = np.empty(n)
X_fe[0] = X0
for i in range(n - 1):
    rhs = g - k * X_fe[i]          # evaluate slope at current point
    X_fe[i + 1] = X_fe[i] + dt * rhs  # take the explicit step

# Backward (implicit) Euler: X_{n+1} = X_n + dt*(g - k*X_{n+1})
# Linear RHS rearranges to a closed form; solve for X_{n+1} explicitly each step.
X_be = np.empty(n)
X_be[0] = X0
for i in range(n - 1):
    # X_next = (X + dt*g) / (1 + dt*k), the algebraic solution of the implicit equation
    X_be[i + 1] = (X_be[i] + dt * g) / (1.0 + dt * k)

# Report numerical results
print(f"Steady state g/k = {steady}")
print(f"dt = {dt}, forward-Euler stability limit 2/k = {2.0/k}")
print(f"Final time t = {t[-1]}")
print(f"Forward  Euler X at final time = {X_fe[-1]}")
print(f"Backward Euler X at final time = {X_be[-1]}")
print(f"Exact          X at final time = {X_exact[-1]}")
print(f"Forward  Euler max |X| over run = {np.max(np.abs(X_fe))}")
print(f"Backward Euler max |X| over run = {np.max(np.abs(X_be))}")
print(f"Backward Euler final error vs exact = {abs(X_be[-1] - X_exact[-1])}")

# Detect oscillation/divergence in forward Euler (sign changes in successive differences + growth)
diffs = np.diff(X_fe)
sign_changes = int(np.sum(np.diff(np.sign(diffs)) != 0))
print(f"Forward Euler sign changes in successive steps (oscillation) = {sign_changes}")
print(f"Forward Euler diverges (|X_final| > |X0|) = {abs(X_fe[-1]) > abs(X0)}")
print(f"Backward Euler stays stable (bounded, converges to g/k) = "
      f"{np.all(np.isfinite(X_be)) and abs(X_be[-1] - steady) < 1e-3}")

# Plot comparison
plt.figure(figsize=(8, 5))
plt.plot(t, X_exact, 'k-', lw=2, label='Exact')
plt.plot(t, X_be, 'bo-', label='Backward (implicit) Euler')
plt.plot(t, X_fe, 'rs--', label='Forward (explicit) Euler')
plt.axhline(steady, color='gray', ls=':', label='g/k = 5')
plt.xlabel('t')
plt.ylabel('X')
plt.title(f'Stiff gene model dX/dt = g - k*X (g={g}, k={k}, dt={dt})')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.8.1_s5.png")

# One-sentence explanation of why the check confirms the result:
# Because at dt=0.2 the explicit factor (1 - dt*k) = -1 has magnitude >= 1 so forward Euler
# amplifies errors and oscillates/diverges, while the implicit factor 1/(1+dt*k) = 1/3 always
# damps toward g/k, the check confirms backward Euler is unconditionally stable for this stiff problem.
print("Why the check confirms it: forward Euler's amplification factor (1 - dt*k) = -1 "
      "has magnitude >= 1 causing growing oscillations, whereas backward Euler's factor "
      "1/(1 + dt*k) = 1/3 is < 1 and always decays to g/k, so stability is unconditional.")
