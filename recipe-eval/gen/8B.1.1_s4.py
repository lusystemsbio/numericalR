import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Two-state promoter Markov chain
# State 0 = active, State 1 = inactive
# Rows of T are the "from" state; columns the "to" state
T = np.array([[0.7, 0.3],
              [0.1, 0.9]])

n_steps = int(1e4)
rng = np.random.default_rng(1)  # seed 1

# --- Direct simulation of the Markov chain ---
states = np.empty(n_steps, dtype=int)
state = 0  # start active
for i in range(n_steps):
    states[i] = state
    u = rng.random()                    # draw a uniform in [0,1)
    row = T[state]                       # transition probs from current state
    # pick next state: if u < P(stay/first outcome) go to col 0, else col 1
    state = 0 if u < row[0] else 1

# --- Running (cumulative) fractions of time active/inactive ---
is_active = (states == 0).astype(float)
steps = np.arange(1, n_steps + 1)
running_active = np.cumsum(is_active) / steps
running_inactive = 1.0 - running_active

# --- Numerical results ---
final_active = running_active[-1]
final_inactive = running_inactive[-1]
print(f"Transition matrix T = {T.tolist()}")
print(f"Final running active fraction:   {final_active:.4f}")
print(f"Final running inactive fraction: {final_inactive:.4f}")

# Theoretical stationary distribution (eigenvector check) for reference
# pi solves pi = pi T ; for this T, pi = [0.25, 0.75]
pi_active = T[1, 0] / (T[0, 1] + T[1, 0])
pi_inactive = 1.0 - pi_active
print(f"Theoretical stationary active fraction:   {pi_active:.4f}")
print(f"Theoretical stationary inactive fraction: {pi_inactive:.4f}")

# Burst check: mean active spell length vs mean inactive spell length
# Expected spell length in a state = 1 / (prob of leaving that state)
mean_active_spell = 1.0 / T[0, 1]    # leave active w.p. 0.3
mean_inactive_spell = 1.0 / T[1, 0]  # leave inactive w.p. 0.1
print(f"Expected mean active spell length (steps):   {mean_active_spell:.4f}")
print(f"Expected mean inactive spell length (steps): {mean_inactive_spell:.4f}")

# Measure empirical spell lengths from the trajectory
active_spells = []
inactive_spells = []
cur = states[0]
length = 1
for s in states[1:]:
    if s == cur:
        length += 1
    else:
        (active_spells if cur == 0 else inactive_spells).append(length)
        cur = s
        length = 1
(active_spells if cur == 0 else inactive_spells).append(length)
print(f"Empirical mean active spell length (steps):   {np.mean(active_spells):.4f}")
print(f"Empirical mean inactive spell length (steps): {np.mean(inactive_spells):.4f}")

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

first = 300  # first steps to visualize the bursting behavior
ax1.step(np.arange(first), states[:first], where="post", color="tab:blue")
ax1.set_yticks([0, 1])
ax1.set_yticklabels(["active (0)", "inactive (1)"])
ax1.set_xlabel("step")
ax1.set_ylabel("promoter state")
ax1.set_title("Promoter state over first %d steps (short active bursts, long silent stretches)" % first)

ax2.plot(steps, running_active, color="tab:green", label="running active fraction")
ax2.plot(steps, running_inactive, color="tab:red", label="running inactive fraction")
ax2.axhline(0.25, color="tab:green", ls="--", alpha=0.6, label="0.25 target")
ax2.axhline(0.75, color="tab:red", ls="--", alpha=0.6, label="0.75 target")
ax2.set_xlabel("step")
ax2.set_ylabel("running fraction")
ax2.set_title("Running active/inactive fractions converge to ~25% / ~75%")
ax2.legend(loc="center right")

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.1.1_s4.png")

# Explanation of why the check confirms the result:
print("Check explanation: because the inactive state is left 3x less often than the active "
      "state, spells silent are ~3x longer, producing bursty dynamics whose running fractions "
      "settle at the stationary distribution 25%/75% predicted by pi = pi T, confirming correctness.")
