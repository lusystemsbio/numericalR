import os
import re
import sys
import shutil
import tempfile
import subprocess
import ctypes

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_PNG = "/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.3.1_s3.png"

# ----------------------------------------------------------------------
# Model / method parameters (fixed by the request)
# ----------------------------------------------------------------------
SEED    = 11        # integer seed -> drives the RNG that lives inside Fortran
TTOT    = 1.0e5     # total integration time for the reported rates
DT      = 0.01      # Euler-Maruyama step
TRELAX  = 100.0     # relaxation time after each crossing
G0, G1  = 10.0, 40.0
THR     = 100.0     # Hill threshold
K       = 0.1       # degradation rate
N       = 4.0       # Hill coefficient per gene
B       = 20.0      # additive noise amplitude

# ----------------------------------------------------------------------
# Fortran source: ONE compiled routine that both R (.Fortran) and Python
# (ctypes) will call.  It implements the mean-first-passage-time method
# EXPLICITLY (Euler-Maruyama + separatrix-crossing detection), and it owns
# its own RNG so that a given seed reproduces the same stream in ANY caller.
# ----------------------------------------------------------------------
FORTRAN_SRC = r"""
      subroutine mfpt(seed, ttot, dt, trelax, g0, g1, thr, kk, en, b,
     &                rate, tau, npass)
      implicit none
      integer seed, npass
      double precision ttot, dt, trelax, g0, g1, thr, kk, en, b
      double precision rate, tau
c     --- host-local state shared with the internal procedures ---
      integer*8 rng_state
      logical have_spare
      double precision spare
      double precision x, y, t, tstart, sumt, sdt
      integer s

c     seed the internal linear-congruential RNG (so both callers agree)
      rng_state = int(seed, 8)
      have_spare = .false.
      sdt = dsqrt(dt)

c     initial condition in the "x-high" basin, then let it settle
      x = 300.0d0
      y = 100.0d0
      t = 0.0d0
      call relax(trelax)

c     record the current basin (sign of x-y) and start the clock
      s      = sign_of(x - y)
      tstart = t
      sumt   = 0.0d0
      npass  = 0

c     main loop: integrate until we cross the separatrix x = y
      do while (t .lt. ttot)
         call step()
         if (sign_of(x - y) .ne. s) then
c           a first passage just completed -> store its duration
            sumt  = sumt + (t - tstart)
            npass = npass + 1
c           relax into the new basin so we start fresh, then re-arm timer
            call relax(trelax)
            s      = sign_of(x - y)
            tstart = t
         end if
      end do

c     mean first-passage time and rate kappa = 1/(2*tau)
      if (npass .gt. 0) then
         tau  = sumt / dble(npass)
         rate = 1.0d0 / (2.0d0 * tau)
      else
         tau  = ttot
         rate = 0.0d0
      end if

      contains

c        one Euler-Maruyama step of the two-gene toggle switch
         subroutine step()
         double precision fx, fy
         fx = g0 + g1*thr**en/(thr**en + y**en) - kk*x
         fy = g0 + g1*thr**en/(thr**en + x**en) - kk*y
         x = x + fx*dt + b*sdt*ran_n()
         y = y + fy*dt + b*sdt*ran_n()
         if (x .lt. 0.0d0) x = 0.0d0
         if (y .lt. 0.0d0) y = 0.0d0
         t = t + dt
         end subroutine step

c        advance by time tr without counting (relaxation)
         subroutine relax(tr)
         double precision tr, tend
         tend = t + tr
         do while (t .lt. tend)
            call step()
         end do
         end subroutine relax

c        which side of the separatrix are we on?
         integer function sign_of(v)
         double precision v
         if (v .ge. 0.0d0) then
            sign_of = 1
         else
            sign_of = -1
         end if
         end function sign_of

c        uniform(0,1) via a 64-bit-safe LCG (glibc constants)
         double precision function ran_u()
         rng_state = mod(rng_state*1103515245_8 + 12345_8, 2147483648_8)
         ran_u = dble(rng_state) / 2147483648.0d0
         end function ran_u

c        standard normal via Box-Muller (cached spare)
         double precision function ran_n()
         double precision u1, u2, r, th
         if (have_spare) then
            have_spare = .false.
            ran_n = spare
         else
            u1 = ran_u()
            u2 = ran_u()
            if (u1 .lt. 1.0d-12) u1 = 1.0d-12
            r  = dsqrt(-2.0d0*dlog(u1))
            th = 2.0d0*3.14159265358979323846d0*u2
            ran_n = r*dcos(th)
            spare = r*dsin(th)
            have_spare = .true.
         end if
         end function ran_n

      end subroutine mfpt
"""

# ----------------------------------------------------------------------
# Build the shared library
# ----------------------------------------------------------------------
if shutil.which("gfortran") is None:
    print("ERROR: gfortran not found; cannot build the shared library.")
    sys.exit(1)

workdir = tempfile.mkdtemp(prefix="toggle_")
src_path = os.path.join(workdir, "toggle.f")
so_path  = os.path.join(workdir, "toggle.so")
with open(src_path, "w") as f:
    f.write(FORTRAN_SRC)

build = subprocess.run(
    ["gfortran", "-O2", "-fPIC", "-shared", "-fno-range-check",
     "-o", so_path, src_path],
    capture_output=True, text=True)
if build.returncode != 0:
    print("ERROR: Fortran compilation failed.")
    print(build.stderr)
    sys.exit(1)
print("Shared library built at:", so_path)

# ----------------------------------------------------------------------
# Python side: call the compiled routine via ctypes
# ----------------------------------------------------------------------
lib = ctypes.CDLL(so_path)
_mfpt = lib.mfpt_          # gfortran mangles "mfpt" -> "mfpt_"
_mfpt.restype = None
_mfpt.argtypes = [ctypes.POINTER(ctypes.c_int)] + \
                 [ctypes.POINTER(ctypes.c_double)]*9 + \
                 [ctypes.POINTER(ctypes.c_double)]*2 + \
                 [ctypes.POINTER(ctypes.c_int)]

def py_mfpt(ttot):
    seed   = ctypes.c_int(SEED)
    a_ttot = ctypes.c_double(ttot)
    dt     = ctypes.c_double(DT)
    trelax = ctypes.c_double(TRELAX)
    g0, g1, thr = ctypes.c_double(G0), ctypes.c_double(G1), ctypes.c_double(THR)
    kk, en, b   = ctypes.c_double(K), ctypes.c_double(N), ctypes.c_double(B)
    rate = ctypes.c_double(0.0)
    tau  = ctypes.c_double(0.0)
    npass = ctypes.c_int(0)
    _mfpt(ctypes.byref(seed), ctypes.byref(a_ttot), ctypes.byref(dt),
          ctypes.byref(trelax), ctypes.byref(g0), ctypes.byref(g1),
          ctypes.byref(thr), ctypes.byref(kk), ctypes.byref(en),
          ctypes.byref(b), ctypes.byref(rate), ctypes.byref(tau),
          ctypes.byref(npass))
    return rate.value, tau.value, npass.value

# main Python result at the requested total time
py_rate, py_tau, py_npass = py_mfpt(TTOT)

# ----------------------------------------------------------------------
# R side: call the SAME compiled routine via .Fortran, same seed
# ----------------------------------------------------------------------
r_rate = r_tau = None
r_npass = None
if shutil.which("Rscript") is None:
    print("NOTE: Rscript not found; skipping the R (.Fortran) call.")
else:
    r_script = os.path.join(workdir, "call.R")
    with open(r_script, "w") as f:
        f.write(f"""
dyn.load("{so_path}")
res <- .Fortran("mfpt",
    seed   = as.integer({SEED}),
    ttot   = as.double({TTOT}),
    dt     = as.double({DT}),
    trelax = as.double({TRELAX}),
    g0     = as.double({G0}),
    g1     = as.double({G1}),
    thr    = as.double({THR}),
    kk     = as.double({K}),
    en     = as.double({N}),
    b      = as.double({B}),
    rate   = as.double(0),
    tau    = as.double(0),
    npass  = as.integer(0))
cat(sprintf("R_RATE %.17e\\n",  res$rate))
cat(sprintf("R_TAU %.17e\\n",   res$tau))
cat(sprintf("R_NPASS %d\\n",    res$npass))
""")
    rr = subprocess.run(["Rscript", "--vanilla", r_script],
                        capture_output=True, text=True)
    if rr.returncode != 0:
        print("NOTE: R call failed; showing Python results only.")
        print(rr.stderr)
    else:
        m_rate  = re.search(r"R_RATE\s+([-\d.eE+]+)",  rr.stdout)
        m_tau   = re.search(r"R_TAU\s+([-\d.eE+]+)",   rr.stdout)
        m_npass = re.search(r"R_NPASS\s+(\d+)",        rr.stdout)
        if m_rate:  r_rate  = float(m_rate.group(1))
        if m_tau:   r_tau   = float(m_tau.group(1))
        if m_npass: r_npass = int(m_npass.group(1))

# ----------------------------------------------------------------------
# Convergence check: longer total time -> more accurate rate (same seed)
# ----------------------------------------------------------------------
ttot_sweep = [1.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e5]
sweep_rates = []
sweep_npass = []
for tt in ttot_sweep:
    rr_rate, rr_tau, rr_np = py_mfpt(tt)
    sweep_rates.append(rr_rate)
    sweep_npass.append(rr_np)

# ----------------------------------------------------------------------
# Print every numerical result, one per line
# ----------------------------------------------------------------------
print("seed:", SEED)
print("total_time:", TTOT)
print("dt:", DT)
print("relaxation_time:", TRELAX)
print("g0:", G0)
print("g1:", G1)
print("threshold:", THR)
print("k:", K)
print("n_per_gene:", N)
print("noise_b:", B)
print("")
print("Python mean_first_passage_time tau: %.17e" % py_tau)
print("Python number_of_passages npass:", py_npass)
print("Python transition_rate kappa = 1/(2*tau): %.17e" % py_rate)
print("")
if r_rate is not None:
    print("R mean_first_passage_time tau: %.17e" % r_tau)
    print("R number_of_passages npass:", r_npass)
    print("R transition_rate kappa = 1/(2*tau): %.17e" % r_rate)
    print("")
    diff = abs(py_rate - r_rate)
    print("abs_difference_python_minus_R_rate: %.3e" % diff)
    print("rates_identical (R == Python for same seed):", bool(py_rate == r_rate))
else:
    print("R transition_rate: unavailable (R not run)")
print("")
print("Convergence sweep (same seed, growing total time):")
for tt, rt, npv in zip(ttot_sweep, sweep_rates, sweep_npass):
    print("  total_time = %.1e  npass = %5d  kappa = %.6e" % (tt, npv, rt))

print("")
print("Explanation: Because R (.Fortran) and Python (ctypes) invoke the exact "
      "same compiled symbol with the same seed, the RNG stream and every "
      "arithmetic step happen inside the shared library, so byte-identical "
      "rates prove the numerics live in the .so and not in either language.")

# ----------------------------------------------------------------------
# Plot: rate vs total time (convergence) with the R value overlaid
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5.0))
ax.semilogx(ttot_sweep, sweep_rates, "o-", color="tab:blue",
            label="Python (ctypes) MFPT rate")
ax.axhline(py_rate, color="tab:blue", ls=":", lw=1,
           label="Python rate @ T=1e5")
if r_rate is not None:
    ax.plot([TTOT], [r_rate], "x", color="tab:red", ms=14, mew=3,
            label="R (.Fortran) rate @ T=1e5")
ax.set_xlabel("total simulation time T")
ax.set_ylabel(r"transition rate $\kappa = 1/(2\tau)$")
ax.set_title("Noisy toggle switch: MFPT transition rate\n"
             "(shared Fortran routine, seed=11; longer T -> more accurate)")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
plt.savefig(OUT_PNG, dpi=130)
print("")
print("Figure saved to:", OUT_PNG)
