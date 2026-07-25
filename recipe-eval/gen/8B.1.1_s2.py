import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Transition matrix: rows are current state, cols are next state
# state 0 = active, state 1 = inactive
T = np.array([[0.7, 0.3],
              [0.1, 0.9]])

n_steps = int(1e4)
rng = np.random.default_rng(1)  # seed 1

# --- Direct simulation of the two-state Markov chain ---
states = np.empty(n_steps + 1, dtype=int)
states[0] = 0  # start active

for i in range(n_steps):
    cur = states[i]
    u = rng.random()               # draw a uniform in [0,1)
    # pick next state from the current row of T:
    # if u < T[cur,0] go to state 0 (active), else state 1 (inactive)
    states[i + 1] = 0 if u < T[cur, 0] else 1

# --- Running fractions of active / inactive over time ---
is_active = (states == 0).astype(float)
steps = np.arange(n_steps + 1)
cum_active = np.cumsum(is_active)
running_active = cum_active / (steps + 1)      # fraction active up to each step
running_inactive = 1.0 - running_active        # fraction inactive up to each step

# --- Report numerical results ---
final_active = running_active[-1]
final_inactive = running_inactive[-1]
print(f"Final running active fraction: {final_active:.4f}")
print(f"Final running inactive fraction: {final_inactive:.4f}")

# Theoretical stationary distribution (for reference/check)
# solve pi = pi T, pi sums to 1
# for a 2-state chain: pi_active = q10/(q01+q10) with q01=0.3, q10=0.1
pi_active = T[1, 0] / (T[0, 1] + T[1, 0])
pi_inactive = 1.0 - pi_active
print(f"Theoretical stationary active fraction: {pi_active:.4f}")
print(f"Theoretical stationary inactive fraction: {pi_inactive:.4f}")

# --- Burst check: measure lengths of consecutive active vs inactive spells ---
def spell_lengths(seq, target):
    lengths = []
    run = 0
    for s in seq:
        if s == target:
            run += 1
        elif run > 0:
            lengths.append(run)
            run = 0
    if run > 0:
        lengths.append(run)
    return np.array(lengths)

active_spells = spell_lengths(states, 0)
inactive_spells = spell_lengths(states, 1)
print(f"Mean active spell length (steps): {active_spells.mean():.4f}")
print(f"Mean inactive spell length (steps): {inactive_spells.mean():.4f}")
print(f"Number of active bursts: {len(active_spells)}")
print(f"Number of inactive stretches: {len(inactive_spells)}")

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

# Top: promoter state over the first steps (show first 300 for visibility)
n_show = 300
ax1.step(steps[:n_show], states[:n_show], where="post", color="tab:blue")
ax1.set_yticks([0, 1])
ax1.set_yticklabels(["active", "inactive"])
ax1.set_xlabel("step")
ax1.set_ylabel("promoter state")
ax1.set_title(f"Promoter state over first {n_show} steps (bursty: short active, long silent)")

# Bottom: running fractions
ax2.plot(steps, running_active, label="running active fraction", color="tab:green")
ax2.plot(steps, running_inactive, label="running inactive fraction", color="tab:red")
ax2.axhline(0.25, color="tab:green", ls="--", lw=0.8, label="0.25 target")
ax2.axhline(0.75, color="tab:red", ls="--", lw=0.8, label="0.75 target")
ax2.set_xlabel("step")
ax2.set_ylabel("fraction")
ax2.set_ylim(0, 1)
ax2.set_title("Running active/inactive fractions converge to 25% / 75%")
ax2.legend(loc="center right", fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.1.1_s2.png")

# --- One-sentence explanation of why the check confirms the result ---
print("Explanation: The check confirms the result because the running fractions "
      "settling near 0.25 active / 0.75 inactive match the chain's exact stationary "
      "distribution (leave-rate ratio 0.1:0.3), and the short active spells with long "
      "silent stretches are exactly the bursting expected from the higher active-to-inactive "
      "transition probability.")
