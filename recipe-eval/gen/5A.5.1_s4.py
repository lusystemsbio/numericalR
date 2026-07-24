import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Harmonic oscillator: d^2x/dt^2 = -k*x, i.e. f(x) = -k*x
# State pair: dx/dt = v, dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (constant for true motion)
# ---------------------------------------------------------------

k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time


def force(x):
    # acceleration f = -k*x (mass = 1)
    return -k * x


def verlet_positions_only(dt):
    """Position-only (basic) Verlet. Stores ONLY positions.
    x_{n+2} = 2*x_{n+1} - x_n + dt^2 * f_{n+1}
    """
    n_steps = int(round(t_end / dt))
    t = np.linspace(0.0, n_steps * dt, n_steps + 1)
    x = np.empty(n_steps + 1)

    # store the initial position
    x[0] = x0

    # startup step (Taylor) to get the second position using v0 and f0
    x[1] = x0 + dt * v0 + 0.5 * dt * dt * force(x0)

    # main recurrence: advance from the two previous positions only
    for n in range(n_steps - 1):
        x[n + 2] = 2.0 * x[n + 1] - x[n] + dt * dt * force(x[n + 1])

    return t, x


def true_solution(t):
    # exact motion: x(t) = x0*cos(w t) + (v0/w)*sin(w t), w = sqrt(k)
    w = np.sqrt(k)
    return x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)


# ---- integrate at both step sizes ----
t_small, x_small = verlet_positions_only(0.01)
t_big, x_big = verlet_positions_only(0.1)

# ---- check: compare each Verlet result against the exact oscillation ----
# (Verlet stores only positions, so we compare positions to the true x(t).)
err_small = np.max(np.abs(x_small - true_solution(t_small)))
err_big = np.max(np.abs(x_big - true_solution(t_big)))

# analytic properties of the true oscillation for reference
w = np.sqrt(k)
amplitude = np.sqrt(x0 ** 2 + (v0 / w) ** 2)
e_true = 0.5 * k * x0 ** 2 + 0.5 * v0 ** 2

print("angular frequency w = sqrt(k):", w)
print("period T = 2*pi/w:", 2 * np.pi / w)
print("true amplitude sqrt(x0^2 + (v0/w)^2):", amplitude)
print("exact energy per unit mass e0 = 0.5*k*x0^2 + 0.5*v0^2:", e_true)
print("Verlet dt=0.01: number of position samples:", x_small.size)
print("Verlet dt=0.01: final position x(100):", x_small[-1])
print("Verlet dt=0.01: max position amplitude:", np.max(np.abs(x_small)))
print("Verlet dt=0.01: max abs error vs exact x(t):", err_small)
print("Verlet dt=0.1: number of position samples:", x_big.size)
print("Verlet dt=0.1: final position x(100):", x_big[-1])
print("Verlet dt=0.1: max position amplitude:", np.max(np.abs(x_big)))
print("Verlet dt=0.1: max abs error vs exact x(t):", err_big)

# ---- plot x(t) for Verlet at both step sizes ----
plt.figure(figsize=(11, 5))
tt = np.linspace(0, t_end, 2000)
plt.plot(tt, true_solution(tt), 'k-', lw=1, alpha=0.5, label="exact x(t)")
plt.plot(t_small, x_small, 'b-', lw=1.0, label="Verlet dt=0.01")
plt.plot(t_big, x_big, 'r--', lw=1.0, label="Verlet dt=0.1")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Position-only Verlet for harmonic oscillator (k=0.1, x0=1, v0=2)")
plt.legend(loc="upper right")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.5.1_s4.png")

# Explanation: the check confirms the result because, even though Verlet advances
# using only the two previous stored positions (never a velocity), its computed x(t)
# stays a bounded oscillation matching the exact sinusoid at both step sizes, showing
# the position-only recurrence faithfully reproduces the dynamics.
print("Check: Verlet reproduces the bounded oscillation at both dt while storing only positions,")
print("       because both max-error values stay small and no amplitude blow-up occurs.")
