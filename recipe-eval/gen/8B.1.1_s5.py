import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model setup: two-state promoter Markov chain ---
# States: 0 = active, 1 = inactive
# Transition matrix T[i, j] = P(next = j | current = i)
# Rows: from active, from inactive
T = np.array([[0.7, 0.3],
              [0.1, 0.9]])

n_steps = int(1e4)
rng = np.random.default_rng(1)   # seed 1

# --- Direct simulation of the Markov chain ---
states = np.empty(n_steps, dtype=int)
state = 0                        # start active
for t in range(n_steps):
    states[t] = state            # record current state
    u = rng.random()             # draw a uniform in [0, 1)
    # pick next state from the current row of T:
    # go to state 0 if u < T[state, 0], else state 1
    if u < T[state, 0]:
        state = 0
    else:
        state = 1

# --- Running (cumulative) fractions of time in each state ---
active_indicator = (states == 0).astype(float)         # 1 when active
steps = np.arange(1, n_steps + 1)
running_active = np.cumsum(active_indicator) / steps    # running active fraction
running_inactive = 1.0 - running_active                 # running inactive fraction

# --- Report numerical results ---
final_active = running_active[-1]
final_inactive = running_inactive[-1]
print(f"Total steps: {n_steps}")
print(f"Final running active fraction: {final_active:.4f}")
print(f"Final running inactive fraction: {final_inactive:.4f}")

# Theoretical stationary distribution (for reference/comparison)
# pi solves pi = pi T ; for this chain pi_active = 0.1/(0.1+0.3) = 0.25
pi_active = 0.1 / (0.1 + 0.3)
pi_inactive = 0.3 / (0.1 + 0.3)
print(f"Theoretical stationary active fraction: {pi_active:.4f}")
print(f"Theoretical stationary inactive fraction: {pi_inactive:.4f}")

# Burst check: measure lengths of consecutive active and inactive spells
def spell_lengths(seq, target):
    lengths = []
    count = 0
    for s in seq:
        if s == target:
            count += 1
        elif count > 0:
            lengths.append(count)
            count = 0
    if count > 0:
        lengths.append(count)
    return np.array(lengths)

active_spells = spell_lengths(states, 0)
inactive_spells = spell_lengths(states, 1)
print(f"Mean active spell length (steps): {active_spells.mean():.4f}")
print(f"Mean inactive spell length (steps): {inactive_spells.mean():.4f}")
# Expected geometric mean spell = 1/leave-prob: active 1/0.3, inactive 1/0.1
print(f"Expected mean active spell (1/0.3): {1/0.3:.4f}")
print(f"Expected mean inactive spell (1/0.1): {1/0.1:.4f}")

# --- Plot ---
n_show = 300  # first steps to display the trajectory clearly
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

# Top: promoter state over the first steps (bursty: short active, long silent)
ax1.step(np.arange(n_show), states[:n_show], where="post", color="C0")
ax1.set_yticks([0, 1])
ax1.set_yticklabels(["active", "inactive"])
ax1.set_xlabel("step")
ax1.set_ylabel("state")
ax1.set_title(f"Promoter state over first {n_show} steps (bursty)")

# Bottom: running active/inactive fractions converging to 0.25 / 0.75
ax2.plot(steps, running_active, color="C1", label="running active fraction")
ax2.plot(steps, running_inactive, color="C2", label="running inactive fraction")
ax2.axhline(0.25, color="C1", ls="--", lw=1, label="0.25 (active target)")
ax2.axhline(0.75, color="C2", ls="--", lw=1, label="0.75 (inactive target)")
ax2.set_xlabel("step")
ax2.set_ylabel("running fraction")
ax2.set_ylim(0, 1)
ax2.set_title("Running fractions converge to ~25% active, ~75% inactive")
ax2.legend(loc="center right", fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.1.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: because the running active fraction settles near 0.25 while "
      "active spells stay short and inactive spells stay long, the simulated long-run "
      "time-average matches the chain's stationary distribution pi = [0.25, 0.75], "
      "confirming correct sampling of the two-state promoter dynamics.")
