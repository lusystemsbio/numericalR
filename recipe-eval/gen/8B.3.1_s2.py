import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Target: stationary distribution of the 8B.1 bursting-gene promoter.
# State 0 = "active", State 1 = "inactive".
# We are given ONLY the target P(x); the dynamics are up to us.
# ---------------------------------------------------------------
P = np.array([0.25, 0.75])   # target P(x) = p_ss

# ---------------------------------------------------------------
# Metropolis-Hastings with a symmetric two-state proposal.
# Proposal always suggests the *other* state, so it is symmetric
# (q(x'|x) = q(x|x') = 1) and the acceptance reduces to
#   a = min(1, P(x')/P(x)).
# ---------------------------------------------------------------
rng = np.random.default_rng(1)   # seed 1
n_steps = int(1e4)

x = 0                    # start "active" (state 0)
chain = np.empty(n_steps + 1, dtype=int)
chain[0] = x

for t in range(n_steps):
    x_prop = 1 - x                       # symmetric proposal: the other state
    a = min(1.0, P[x_prop] / P[x])       # MH acceptance ratio
    if rng.random() < a:                 # accept with probability a
        x = x_prop
    # else: reject, stay at x (implicit)
    chain[t + 1] = x

# ---------------------------------------------------------------
# Sampled state fractions (should approach the target P(x)).
# ---------------------------------------------------------------
frac0 = np.mean(chain == 0)
frac1 = np.mean(chain == 1)
print("Target P(x)                : active=%.4f  inactive=%.4f" % (P[0], P[1]))
print("Sampled fraction active    : %.4f" % frac0)
print("Sampled fraction inactive  : %.4f" % frac1)

# ---------------------------------------------------------------
# Recover the transition matrix empirically by counting the
# transitions x_t -> x_{t+1} realized by the sampler.
# ---------------------------------------------------------------
T_hat = np.zeros((2, 2))
for i in range(n_steps):
    T_hat[chain[i], chain[i + 1]] += 1
T_hat = T_hat / T_hat.sum(axis=1, keepdims=True)   # normalize each row

print("Recovered MH transition matrix (empirical):")
print("  T_hat[0]: %.4f  %.4f" % (T_hat[0, 0], T_hat[0, 1]))
print("  T_hat[1]: %.4f  %.4f" % (T_hat[1, 0], T_hat[1, 1]))

# Analytic MH transition matrix for comparison:
#   from 0: propose 1, accept min(1,0.75/0.25)=1        -> [0, 1]
#   from 1: propose 0, accept min(1,0.25/0.75)=1/3      -> [1/3, 2/3]
T_mh = np.array([[0.0, 1.0],
                 [1.0/3.0, 2.0/3.0]])
print("Analytic MH transition matrix:")
print("  T_mh[0]: %.4f  %.4f" % (T_mh[0, 0], T_mh[0, 1]))
print("  T_mh[1]: %.4f  %.4f" % (T_mh[1, 0], T_mh[1, 1]))

# ---------------------------------------------------------------
# Separate check: the ORIGINAL 8B.1 chain has a DIFFERENT transition
# matrix, yet shares the SAME stationary distribution p_ss.
# (For a 2-state chain, p_ss is fixed once T[1,0]/T[0,1] = p0/p1 = 1/3.)
# ---------------------------------------------------------------
T_orig = np.array([[0.7, 0.3],
                   [0.1, 0.9]])   # a representative original chain, p0/p1 = 0.1/0.3 = 1/3

# Stationary distribution of a 2-state chain: proportional to (T[1,0], T[0,1]).
def stationary_2state(T):
    v = np.array([T[1, 0], T[0, 1]])
    return v / v.sum()

ss_mh = stationary_2state(T_mh)
ss_orig = stationary_2state(T_orig)
matrices_differ = not np.allclose(T_mh, T_orig)
same_stationary = np.allclose(ss_mh, ss_orig) and np.allclose(ss_mh, P)

print("Original 8B.1 transition matrix:")
print("  T_orig[0]: %.4f  %.4f" % (T_orig[0, 0], T_orig[0, 1]))
print("  T_orig[1]: %.4f  %.4f" % (T_orig[1, 0], T_orig[1, 1]))
print("Stationary of MH chain     : active=%.4f  inactive=%.4f" % (ss_mh[0], ss_mh[1]))
print("Stationary of original     : active=%.4f  inactive=%.4f" % (ss_orig[0], ss_orig[1]))
print("Transition matrices differ : %s" % matrices_differ)
print("Same stationary as target  : %s" % same_stationary)
# This check confirms the result because two chains with different transition
# matrices but the identical stationary distribution show that MH only needs to
# preserve p_ss (not the original dynamics) to sample the correct target.

# ---------------------------------------------------------------
# Plot: sampled fractions vs target, and running estimate of P(active).
# ---------------------------------------------------------------
running = np.cumsum(chain == 0) / np.arange(1, len(chain) + 1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

ax1.bar([0, 1], [frac0, frac1], width=0.5, color="steelblue", label="MH sampled")
ax1.plot([0, 1], P, "ro", markersize=10, label="target p_ss")
ax1.set_xticks([0, 1]); ax1.set_xticklabels(["active", "inactive"])
ax1.set_ylabel("fraction"); ax1.set_title("Sampled vs target distribution")
ax1.legend()

ax2.plot(running, color="steelblue", label="running P(active)")
ax2.axhline(P[0], color="red", linestyle="--", label="target 0.25")
ax2.set_xlabel("step"); ax2.set_ylabel("P(active)")
ax2.set_title("Convergence of MH estimate"); ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.3.1_s2.png")
