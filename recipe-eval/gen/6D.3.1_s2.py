import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
import re
import shutil
import subprocess
import tempfile
import ctypes
from ctypes import byref, c_int, c_double

# ---------------------------------------------------------------------------
# Model / measurement parameters (shared by R and Python)
# ---------------------------------------------------------------------------
SEED    = 11        # RNG seed -> lives inside the Fortran routine
TTOT    = 1.0e5     # total integration time for the reported rate
DT      = 0.01      # Euler-Maruyama time step
TRELAX  = 100.0     # relaxation time after each crossing
G0, G1  = 10.0, 40.0
THR     = 100.0     # Hill threshold
KDEG    = 0.1       # degradation rate
NEXP    = 4.0       # Hill coefficient per gene
BNOISE  = 20.0      # noise amplitude

workdir = tempfile.mkdtemp(prefix="toggle_mfpt_")
f90path = os.path.join(workdir, "mfpt.f90")
libpath = os.path.join(workdir, "libmfpt.so")
rpath   = os.path.join(workdir, "run_mfpt.R")

# ---------------------------------------------------------------------------
# The compiled Fortran routine.  All randomness is generated *inside* the
# shared library from the integer seed (a portable Park-Miller LCG feeding a
# Box-Muller Gaussian), so R (.Fortran) and Python (ctypes) that call the same
# .so with the same seed must produce byte-identical trajectories and rates.
# The mean-first-passage method is written out explicitly, step by step.
# ---------------------------------------------------------------------------
FORTRAN_SRC = r"""
      subroutine punif(state, u)
      implicit none
      integer*8, intent(inout) :: state
      double precision, intent(out) :: u
      ! minimal-standard (Park-Miller) LCG, 16807 * state mod 2^31-1
      state = mod(16807_8 * state, 2147483647_8)
      u = dble(state) / 2147483647.0d0
      end subroutine

      subroutine pgauss2(state, g1, g2)
      implicit none
      integer*8, intent(inout) :: state
      double precision, intent(out) :: g1, g2
      double precision :: u1, u2, r, th
      ! Box-Muller: two uniforms -> two standard normals
      call punif(state, u1)
      call punif(state, u2)
      if (u1 < 1.0d-12) u1 = 1.0d-12
      r  = sqrt(-2.0d0*log(u1))
      th = 6.283185307179586d0*u2
      g1 = r*cos(th)
      g2 = r*sin(th)
      end subroutine

      subroutine estep(x1,x2,dt,sq,g0,g1,thr,kdeg,nexp,bnoise,state)
      implicit none
      double precision, intent(inout) :: x1, x2
      double precision, intent(in)    :: dt,sq,g0,g1,thr,kdeg,nexp,bnoise
      integer*8, intent(inout)        :: state
      double precision :: h1,h2,n1,n2,tn
      ! one Euler-Maruyama step of the mutually repressing toggle switch
      tn = thr**nexp
      h1 = g0 + g1*tn/(tn + x2**nexp)    ! production of x1, repressed by x2
      h2 = g0 + g1*tn/(tn + x1**nexp)    ! production of x2, repressed by x1
      call pgauss2(state, n1, n2)
      x1 = x1 + (h1 - kdeg*x1)*dt + bnoise*sq*n1
      x2 = x2 + (h2 - kdeg*x2)*dt + bnoise*sq*n2
      if (x1 < 0.0d0) x1 = 0.0d0
      if (x2 < 0.0d0) x2 = 0.0d0
      end subroutine

      subroutine mfpt(seed, ttot, dt, trelax, g0, g1, thr, kdeg,
     &                nexp, bnoise, tau, ncross, rate)
      implicit none
      integer, intent(in)  :: seed
      integer, intent(out) :: ncross
      double precision, intent(in)  :: ttot,dt,trelax,g0,g1,thr,kdeg,
     &                                 nexp,bnoise
      double precision, intent(out) :: tau, rate
      integer*8 :: state
      double precision :: x1, x2, t, tstart, sq, tsum, side, trend

      state  = int(seed, kind=8)   ! seed the shared-library RNG
      x1     = 200.0d0             ! start in the (x1 high, x2 low) basin
      x2     = 20.0d0
      sq     = sqrt(dt)
      t      = 0.0d0
      tstart = 0.0d0
      tsum   = 0.0d0
      ncross = 0
      side   = sign(1.0d0, x1 - x2)   ! which side of separatrix x1=x2

      do while (t < ttot)
         ! advance one step and check for a separatrix crossing
         call estep(x1,x2,dt,sq,g0,g1,thr,kdeg,nexp,bnoise,state)
         t = t + dt
         if (sign(1.0d0, x1 - x2) /= side) then
            ! crossed x1=x2: record this first-passage time
            tsum   = tsum + (t - tstart)
            ncross = ncross + 1
            side   = sign(1.0d0, x1 - x2)
            ! relax so the trajectory settles into the new basin,
            ! then start timing the next passage from t
            trend = t + trelax
            do while (t < trend .and. t < ttot)
               call estep(x1,x2,dt,sq,g0,g1,thr,kdeg,nexp,bnoise,state)
               t    = t + dt
               side = sign(1.0d0, x1 - x2)
            end do
            tstart = t
         end if
      end do

      if (ncross > 0) then
         tau  = tsum / dble(ncross)   ! mean first-passage time
         rate = 1.0d0 / (2.0d0*tau)   ! kappa = 1/(2*tau)
      else
         tau  = 0.0d0
         rate = 0.0d0
      end if
      end subroutine
"""

with open(f90path, "w") as f:
    f.write(FORTRAN_SRC)

# --- compile the shared library used by BOTH languages ---------------------
gfortran = shutil.which("gfortran")
if gfortran is None:
    raise RuntimeError("gfortran not found; required to build the shared library.")
subprocess.run([gfortran, "-O2", "-shared", "-fPIC", "-o", libpath, f90path],
               check=True)

# ---------------------------------------------------------------------------
# Python side: call the compiled routine through ctypes (pass everything by
# reference, as Fortran expects).
# ---------------------------------------------------------------------------
lib = ctypes.CDLL(libpath)
fn = lib.mfpt_   # gfortran name mangling appends an underscore

def run_python(seed, ttot, dt, trelax):
    s   = c_int(seed);     tt = c_double(ttot); h = c_double(dt)
    tr  = c_double(trelax)
    g0  = c_double(G0);    g1 = c_double(G1);   th = c_double(THR)
    kd  = c_double(KDEG);  ne = c_double(NEXP); bn = c_double(BNOISE)
    tau = c_double(0.0);   nc = c_int(0);       rt = c_double(0.0)
    fn(byref(s), byref(tt), byref(h), byref(tr),
       byref(g0), byref(g1), byref(th), byref(kd), byref(ne), byref(bn),
       byref(tau), byref(nc), byref(rt))
    return tau.value, nc.value, rt.value

tau_py, ncross_py, rate_py = run_python(SEED, TTOT, DT, TRELAX)

# ---------------------------------------------------------------------------
# R side: call the *same* .so through .Fortran with the same seed.
# ---------------------------------------------------------------------------
R_SRC = r'''
args <- commandArgs(TRUE)
dyn.load(args[1])
res <- .Fortran("mfpt",
  seed   = as.integer(%d),
  ttot   = as.double(%.10g),
  dt     = as.double(%.10g),
  trelax = as.double(%.10g),
  g0     = as.double(%.10g),
  g1     = as.double(%.10g),
  thr    = as.double(%.10g),
  kdeg   = as.double(%.10g),
  nexp   = as.double(%.10g),
  bnoise = as.double(%.10g),
  tau    = as.double(0),
  ncross = as.integer(0),
  rate   = as.double(0))
cat(sprintf("TAU_R %%.15e\n",  res$tau))
cat(sprintf("NCROSS_R %%d\n",  res$ncross))
cat(sprintf("RATE_R %%.15e\n", res$rate))
''' % (SEED, TTOT, DT, TRELAX, G0, G1, THR, KDEG, NEXP, BNOISE)

with open(rpath, "w") as f:
    f.write(R_SRC)

rate_R = tau_R = None
ncross_R = None
Rscript = shutil.which("Rscript")
if Rscript is None:
    print("NOTE: Rscript not found; R result unavailable in this environment.")
else:
    out = subprocess.run([Rscript, rpath, libpath],
                         capture_output=True, text=True)
    txt = out.stdout + out.stderr
    m = re.search(r"RATE_R\s+([-\d.eE+]+)", txt)
    if m:
        rate_R = float(m.group(1))
    m = re.search(r"TAU_R\s+([-\d.eE+]+)", txt)
    if m:
        tau_R = float(m.group(1))
    m = re.search(r"NCROSS_R\s+(\d+)", txt)
    if m:
        ncross_R = int(m.group(1))
    if rate_R is None:
        print("NOTE: could not parse R output:")
        print(txt)

# ---------------------------------------------------------------------------
# Convergence check: longer total time -> more crossings -> more accurate rate.
# ---------------------------------------------------------------------------
ttots = [1.0e4, 3.0e4, 1.0e5, 3.0e5]
conv_rates = []
conv_ncross = []
for tt in ttots:
    _, nc, rt = run_python(SEED, tt, DT, TRELAX)
    conv_rates.append(rt)
    conv_ncross.append(nc)

# ---------------------------------------------------------------------------
# Report every numerical result, one per line.
# ---------------------------------------------------------------------------
print("Python tau (mean first-passage time):        %.15e" % tau_py)
print("Python number of crossings:                  %d" % ncross_py)
print("Python transition rate kappa = 1/(2*tau):    %.15e" % rate_py)
if rate_R is not None:
    print("R tau (mean first-passage time):             %.15e" % tau_R)
    print("R number of crossings:                       %d" % ncross_R)
    print("R transition rate kappa = 1/(2*tau):         %.15e" % rate_R)
    print("Absolute rate difference |R - Python|:       %.3e" % abs(rate_R - rate_py))
    print("Rates identical (exact bit match):           %s"
          % str(rate_R == rate_py))
else:
    print("R transition rate kappa:                     UNAVAILABLE")

print("Convergence of Python rate with total time:")
for tt, rt, nc in zip(ttots, conv_rates, conv_ncross):
    print("  ttot = %.3e : ncross = %6d : kappa = %.15e" % (tt, nc, rt))
print("Most accurate (largest ttot) kappa:          %.15e" % conv_rates[-1])

# ---------------------------------------------------------------------------
# Plot: rate vs total time, with R and Python points at ttot = 1e5.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(ttots, conv_rates, "o-", color="steelblue",
            label="Python (ctypes) vs total time")
ax.semilogx([TTOT], [rate_py], "s", color="darkorange", markersize=11,
            label="Python @ ttot=1e5")
if rate_R is not None:
    ax.semilogx([TTOT], [rate_R], "x", color="black", markersize=13,
                markeredgewidth=3, label="R (.Fortran) @ ttot=1e5")
ax.set_xlabel("total simulation time")
ax.set_ylabel(r"transition rate $\kappa = 1/(2\tau)$")
ax.set_title("Toggle-switch transition rate from mean first-passage time")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.3.1_s2.png")

# One-sentence explanation of why the cross-language check confirms the result:
print("Explanation: because R and Python both invoke the same compiled .so "
      "with the same seed and the RNG lives inside that library, identical "
      "rates prove the numerics come entirely from the shared Fortran routine, "
      "not from either language's own code.")
