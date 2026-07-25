import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Target distribution (stationary distribution of the 8B.1 bursting-gene promoter) ----
# States: 0 = inactive, 1 = active
p_ss = np.array([0.25, 0.75])   # target P(x)

# ---- Metropolis-Hastings with symmetric two-state proposal ----
rng = np.random.default_rng(1)  # seed 1
n_steps = int(1e4)

x = 1                # start active (state 1)
chain = np.empty(n_steps + 1, dtype=int)
chain[0] = x

for t in range(n_steps):
    x_prop = 1 - x                         # symmetric proposal: always propose the other state
    a = min(1.0, p_ss[x_prop] / p_ss[x])   # acceptance prob (symmetric => ratio of targets)
    if rng.random() < a:                   # accept/reject
        x = x_prop
    chain[t + 1] = x                        # record current state (accepted or held)

# ---- Sampled state fractions ----
frac0 = np.mean(chain == 0)
frac1 = np.mean(chain == 1)
print(f"Sampled fraction state 0 (inactive): {frac0:.4f}")
print(f"Sampled fraction state 1 (active):   {frac1:.4f}")
print(f"Target p_ss:                         [{p_ss[0]:.4f}, {p_ss[1]:.4f}]")

# ---- Recover transition matrix from the sampled chain (empirical MC estimate) ----
counts = np.zeros((2, 2))
for t in range(n_steps):
    counts[chain[t], chain[t + 1]] += 1
T_sampler = counts / counts.sum(axis=1, keepdims=True)
print("Recovered (empirical) transition matrix of MH sampler:")
print(f"  T[0,0]={T_sampler[0,0]:.4f}  T[0,1]={T_sampler[0,1]:.4f}")
print(f"  T[1,0]={T_sampler[1,0]:.4f}  T[1,1]={T_sampler[1,1]:.4f}")

# ---- Theoretical MH transition matrix (for comparison) ----
# From state 0: propose 1, a = min(1, .75/.25)=1  -> always move  => T[0,0]=0, T[0,1]=1
# From state 1: propose 0, a = min(1, .25/.75)=1/3 -> move w.p. 1/3 => T[1,1]=2/3, T[1,0]=1/3
T_mh = np.array([[0.0, 1.0],
                 [1.0/3.0, 2.0/3.0]])
print("Theoretical MH transition matrix:")
print(f"  T[0,0]={T_mh[0,0]:.4f}  T[0,1]={T_mh[0,1]:.4f}")
print(f"  T[1,0]={T_mh[1,0]:.4f}  T[1,1]={T_mh[1,1]:.4f}")

# ---- Separate check: a DIFFERENT original chain (8B.1) with the SAME stationary distribution ----
# Any chain satisfying detailed balance k01*p0 = k10*p1 shares p_ss. Pick a different one.
# Choose T_orig[0,1]=0.6 -> T_orig[0,0]=0.4 ; stationarity: 0.25*0.6 = 0.75*T_orig[1,0]
t10 = (0.25 * 0.6) / 0.75           # = 0.2
T_orig = np.array([[0.4, 0.6],
                   [t10, 1.0 - t10]])
print("Original 8B.1-style chain transition matrix (different from MH):")
print(f"  T[0,0]={T_orig[0,0]:.4f}  T[0,1]={T_orig[0,1]:.4f}")
print(f"  T[1,0]={T_orig[1,0]:.4f}  T[1,1]={T_orig[1,1]:.4f}")

# ---- Confirm both matrices are different but yield the same stationary distribution ----
def stationary(T):
    vals, vecs = np.linalg.eig(T.T)          # left eigenvector for eigenvalue 1
    v = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    return v / v.sum()

pi_mh = stationary(T_mh)
pi_orig = stationary(T_orig)
matrices_differ = not np.allclose(T_mh, T_orig)
same_stationary = np.allclose(pi_mh, pi_orig) and np.allclose(pi_mh, p_ss)
print(f"Stationary dist of MH matrix:       [{pi_mh[0]:.4f}, {pi_mh[1]:.4f}]")
print(f"Stationary dist of original matrix: [{pi_orig[0]:.4f}, {pi_orig[1]:.4f}]")
print(f"Matrices differ:                    {matrices_differ}")
print(f"Same stationary distribution:       {same_stationary}")

# Explanation: two demonstrably different transition matrices both leave p_ss invariant,
# so the stationary distribution does not uniquely determine the dynamics -- many chains
# (the original 8B.1 chain and the MH sampler) share one stationary distribution.
print("Why the check confirms the result: two different transition matrices both fix the "
      "same stationary p_ss, proving many distinct chains can share one stationary distribution.")

# ---- Plot ----
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
labels = ["inactive (0)", "active (1)"]
xpos = np.arange(2)
ax[0].bar(xpos - 0.2, [frac0, frac1], width=0.4, label="MH sampled")
ax[0].bar(xpos + 0.2, p_ss, width=0.4, label="target p_ss")
ax[0].set_xticks(xpos); ax[0].set_xticklabels(labels)
ax[0].set_ylabel("probability"); ax[0].set_title("Sampled vs target distribution")
ax[0].legend()

ax[1].plot(chain[:200], drawstyle="steps-post")
ax[1].set_xlabel("step"); ax[1].set_ylabel("state")
ax[1].set_yticks([0, 1]); ax[1].set_title("First 200 MH states")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.3.1_s4.png")
