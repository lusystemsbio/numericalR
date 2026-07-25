import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Target: two-state bursting-gene promoter of 8B.1, now given ONLY by its
# stationary distribution.  State 0 = inactive, State 1 = active.
# ----------------------------------------------------------------------
P = np.array([0.25, 0.75])   # target P(x): P(inactive)=0.25, P(active)=0.75

# ----------------------------------------------------------------------
# Metropolis-Hastings with a symmetric two-state proposal.
# Proposal is deterministic: from any state we always propose the OTHER
# state, so the proposal is symmetric and acceptance is a = min(1,P(x')/P(x)).
# ----------------------------------------------------------------------
np.random.seed(1)            # reproducible run, seed 1
n_steps = int(1e4)           # 1e4 steps
x = 1                        # start active (state 1)

chain = np.empty(n_steps + 1, dtype=int)
chain[0] = x
for t in range(n_steps):
    x_prop = 1 - x                       # symmetric proposal: always the other state
    a = min(1.0, P[x_prop] / P[x])       # MH acceptance ratio
    if np.random.rand() < a:             # accept with probability a
        x = x_prop                       # move
    # else: reject and stay at x (implicitly)
    chain[t + 1] = x

# ----------------------------------------------------------------------
# Sampled state fractions (drop the initial state; report the visited states)
# ----------------------------------------------------------------------
samples = chain[1:]
frac = np.array([np.mean(samples == 0), np.mean(samples == 1)])
print(f"Sampled fraction inactive (state 0): {frac[0]:.4f}")
print(f"Sampled fraction   active (state 1): {frac[1]:.4f}")
print(f"Target stationary distribution p_ss: {P[0]:.4f}, {P[1]:.4f}")

# ----------------------------------------------------------------------
# Recover the sampler's transition matrix by counting empirical transitions.
# ----------------------------------------------------------------------
counts = np.zeros((2, 2))
for i in range(len(chain) - 1):
    counts[chain[i], chain[i + 1]] += 1
T_recovered = counts / counts.sum(axis=1, keepdims=True)
print("Recovered (empirical) MH transition matrix:")
print(f"  T[0,0]={T_recovered[0,0]:.4f}  T[0,1]={T_recovered[0,1]:.4f}")
print(f"  T[1,0]={T_recovered[1,0]:.4f}  T[1,1]={T_recovered[1,1]:.4f}")

# ----------------------------------------------------------------------
# Analytic MH transition matrix implied by the deterministic proposal + a=min(1,ratio):
#   from 0: propose 1, ratio P1/P0 = 3 >= 1 -> accept surely: T[0,1]=1
#   from 1: propose 0, ratio P0/P1 = 1/3    -> accept w.p. 1/3: T[1,0]=1/3
# ----------------------------------------------------------------------
T_MH = np.array([[0.0,      1.0],
                 [1.0/3.0,  2.0/3.0]])
print("Analytic MH transition matrix:")
print(f"  T[0,0]={T_MH[0,0]:.4f}  T[0,1]={T_MH[0,1]:.4f}")
print(f"  T[1,0]={T_MH[1,0]:.4f}  T[1,1]={T_MH[1,1]:.4f}")

# ----------------------------------------------------------------------
# The ORIGINAL bursting-gene chain of 8B.1 (an example dynamics that ALSO
# has p_ss=(0.25,0.75) via detailed balance 0.25*k_on = 0.75*k_off).
# ----------------------------------------------------------------------
T_orig = np.array([[0.4, 0.6],
                   [0.2, 0.8]])
print("Original 8B.1 chain transition matrix:")
print(f"  T[0,0]={T_orig[0,0]:.4f}  T[0,1]={T_orig[0,1]:.4f}")
print(f"  T[1,0]={T_orig[1,0]:.4f}  T[1,1]={T_orig[1,1]:.4f}")

# ----------------------------------------------------------------------
# Check: the two chains are DIFFERENT but share the SAME stationary vector.
# Stationary pi solves pi = pi T ; obtain it as the left eigenvector (eig=1).
# ----------------------------------------------------------------------
def stationary(T):
    vals, vecs = np.linalg.eig(T.T)          # left eigenvectors of T
    v = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    return v / v.sum()

pi_MH = stationary(T_MH)
pi_orig = stationary(T_orig)
matrices_differ = not np.allclose(T_MH, T_orig)
same_stationary = np.allclose(pi_MH, pi_orig, atol=1e-8)
print(f"Stationary of MH chain:       {pi_MH[0]:.4f}, {pi_MH[1]:.4f}")
print(f"Stationary of original chain: {pi_orig[0]:.4f}, {pi_orig[1]:.4f}")
print(f"Transition matrices differ:   {matrices_differ}")
print(f"Same stationary distribution: {same_stationary}")

# One-sentence explanation of why this check confirms the result:
# Because the MH chain and the original chain have distinct transition matrices
# yet identical stationary vectors, the stationary distribution alone does not
# fix the dynamics -- many chains share one stationary distribution.

# ----------------------------------------------------------------------
# Figure: compare sampled fractions with the target.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
xpos = np.arange(2)
ax.bar(xpos - 0.2, P, width=0.4, label="target p_ss", color="steelblue")
ax.bar(xpos + 0.2, frac, width=0.4, label="MH samples", color="darkorange")
ax.set_xticks(xpos)
ax.set_xticklabels(["inactive (0)", "active (1)"])
ax.set_ylabel("probability / fraction")
ax.set_title("Metropolis-Hastings reproduces the target p_ss")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.3.1_s5.png")
