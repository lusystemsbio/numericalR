import os
import sys
import shutil
import subprocess
import ctypes
from ctypes import c_int, c_double, byref, POINTER

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Output locations
# ----------------------------------------------------------------------
PNG = "/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.3.1_s4.png"
WORK = os.path.dirname(PNG)
os.makedirs(WORK, exist_ok=True)
F90 = os.path.join(WORK, "toggle.f90")
LIB = os.path.join(WORK, "libtoggle.so")
RSCRIPT = os.path.join(WORK, "toggle_call.R")

# ----------------------------------------------------------------------
# Fortran source: the numerics (RNG + SDE integrator) live ENTIRELY in
# this compiled shared library, so R (.Fortran) and Python (ctypes) call
# byte-for-byte the same code path.  A self-written Park-Miller LCG keyed
# only by an integer seed guarantees identical random streams in both.
# ----------------------------------------------------------------------
FORTRAN = r"""
module rng_mod
  implicit none
  integer(kind=8) :: rng_state
contains
  subroutine rng_seed(s)
    integer, intent(in) :: s
    rng_state = int(s, kind=8)
    if (rng_state <= 0_8) rng_state = 1_8
  end subroutine rng_seed

  double precision function rng_uniform()          ! Park-Miller minimal standard
    rng_state = mod(16807_8 * rng_state, 2147483647_8)
    rng_uniform = dble(rng_state) / 2147483647.0d0
  end function rng_uniform

  double precision function rng_gauss()            ! Box-Muller standard normal
    double precision :: u1, u2
    double precision, parameter :: twopi = 6.283185307179586d0
    u1 = rng_uniform()
    u2 = rng_uniform()
    if (u1 < 1.0d-12) u1 = 1.0d-12
    rng_gauss = sqrt(-2.0d0*log(u1)) * cos(twopi*u2)
  end function rng_gauss
end module rng_mod

! Mean first-passage-time estimate of the toggle-switch transition rate.
subroutine toggle_mfpt(seed, ttot, dt, trelax, g0, g1, Kthr, kdeg, nexp, b, &
                       tau, rate, ncross)
  use rng_mod
  implicit none
  integer :: seed, nexp, ncross
  double precision :: ttot, dt, trelax, g0, g1, Kthr, kdeg, b, tau, rate
  double precision :: x1, x2, h1, h2, sq, Kn, t, tstart, sumfpt
  double precision :: prevdiff, curdiff
  integer :: nrelax, i

  call rng_seed(seed)
  sq = b * sqrt(dt)                ! noise increment scale  b*sqrt(dt)
  Kn = Kthr ** nexp               ! threshold^n (Hill numerator)
  nrelax = int(trelax/dt)         ! number of steps to relax after a crossing

  ! start in the high-low stable state (basin 1); separatrix is x1 = x2
  x1 = 500.0d0
  x2 = 100.0d0
  t = 0.0d0
  tstart = 0.0d0
  sumfpt = 0.0d0
  ncross = 0
  prevdiff = x1 - x2

  do while (t < ttot)
     ! --- one Euler-Maruyama step of the two-gene toggle SDE ---
     h1 = Kn / (Kn + x2**nexp)     ! repression of gene 1 by product 2
     h2 = Kn / (Kn + x1**nexp)     ! repression of gene 2 by product 1
     x1 = x1 + (g0 + g1*h1 - kdeg*x1)*dt + sq*rng_gauss()
     x2 = x2 + (g0 + g1*h2 - kdeg*x2)*dt + sq*rng_gauss()
     if (x1 < 0.0d0) x1 = -x1      ! reflect at 0 (concentrations >= 0)
     if (x2 < 0.0d0) x2 = -x2
     t = t + dt
     curdiff = x1 - x2

     ! --- detect a separatrix crossing (sign change of x1 - x2) ---
     if (curdiff * prevdiff < 0.0d0) then
        sumfpt = sumfpt + (t - tstart)   ! record this first-passage time
        ncross = ncross + 1
        ! relax into the new basin so we do not re-count the same event
        do i = 1, nrelax
           h1 = Kn / (Kn + x2**nexp)
           h2 = Kn / (Kn + x1**nexp)
           x1 = x1 + (g0 + g1*h1 - kdeg*x1)*dt + sq*rng_gauss()
           x2 = x2 + (g0 + g1*h2 - kdeg*x2)*dt + sq*rng_gauss()
           if (x1 < 0.0d0) x1 = -x1
           if (x2 < 0.0d0) x2 = -x2
           t = t + dt
        end do
        tstart = t                        ! restart the clock after relaxation
        prevdiff = x1 - x2
     else
        prevdiff = curdiff
     end if
  end do

  if (ncross > 0) then
     tau  = sumfpt / dble(ncross)         ! mean first-passage time
     rate = 1.0d0 / (2.0d0 * tau)         ! transition rate  kappa = 1/(2 tau)
  else
     tau  = -1.0d0
     rate = 0.0d0
  end if
end subroutine toggle_mfpt

! Store a (strided) sample trajectory for plotting only.
subroutine toggle_traj(seed, nsteps, stride, dt, g0, g1, Kthr, kdeg, nexp, b, &
                       nout, o1, o2)
  use rng_mod
  implicit none
  integer :: seed, nsteps, stride, nexp, nout
  double precision :: dt, g0, g1, Kthr, kdeg, b
  double precision :: o1(nout), o2(nout)
  double precision :: x1, x2, h1, h2, sq, Kn
  integer :: i, j
  call rng_seed(seed)
  sq = b * sqrt(dt)
  Kn = Kthr ** nexp
  x1 = 500.0d0
  x2 = 100.0d0
  j = 0
  do i = 1, nsteps
     h1 = Kn / (Kn + x2**nexp)
     h2 = Kn / (Kn + x1**nexp)
     x1 = x1 + (g0 + g1*h1 - kdeg*x1)*dt + sq*rng_gauss()
     x2 = x2 + (g0 + g1*h2 - kdeg*x2)*dt + sq*rng_gauss()
     if (x1 < 0.0d0) x1 = -x1
     if (x2 < 0.0d0) x2 = -x2
     if (mod(i, stride) == 0) then
        j = j + 1
        if (j <= nout) then
           o1(j) = x1
           o2(j) = x2
        end if
     end if
  end do
end subroutine toggle_traj
"""

# ----------------------------------------------------------------------
# Model parameters (shared by R and Python; passed into the library)
# ----------------------------------------------------------------------
SEED    = 11
G0      = 10.0
G1      = 40.0
KTHR    = 100.0     # threshold
KDEG    = 0.1
NEXP    = 4         # Hill coefficient per gene
BNOISE  = 20.0
DT      = 0.01
TRELAX  = 100.0
TTOT    = 1.0e5
TTOT_LONG = 3.0e5

print("=== Noisy toggle switch: transition rate from mean first-passage time ===")
print(f"seed                 = {SEED}")
print(f"g0                   = {G0}")
print(f"g1                   = {G1}")
print(f"threshold K          = {KTHR}")
print(f"k (degradation)      = {KDEG}")
print(f"n (Hill exponent)    = {NEXP}")
print(f"noise b              = {BNOISE}")
print(f"dt                   = {DT}")
print(f"relaxation time      = {TRELAX}")
print(f"total time           = {TTOT}")

# ----------------------------------------------------------------------
# Compile the Fortran into a shared library
# ----------------------------------------------------------------------
gfortran = shutil.which("gfortran")
if gfortran is None:
    print("ERROR: gfortran not found; cannot build the shared library.")
    sys.exit(1)

with open(F90, "w") as fh:
    fh.write(FORTRAN)

subprocess.run([gfortran, "-O2", "-shared", "-fPIC", "-o", LIB, F90], check=True)
print(f"compiled shared library: {LIB}")

# ----------------------------------------------------------------------
# Python side: call the compiled routine via ctypes
# ----------------------------------------------------------------------
lib = ctypes.CDLL(LIB)

lib.toggle_mfpt_.argtypes = [
    POINTER(c_int),                                     # seed
    POINTER(c_double), POINTER(c_double), POINTER(c_double),  # ttot, dt, trelax
    POINTER(c_double), POINTER(c_double), POINTER(c_double),  # g0, g1, Kthr
    POINTER(c_double), POINTER(c_int), POINTER(c_double),     # kdeg, nexp, b
    POINTER(c_double), POINTER(c_double), POINTER(c_int),     # tau, rate, ncross
]
lib.toggle_mfpt_.restype = None

nd = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
lib.toggle_traj_.argtypes = [
    POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double),
    POINTER(c_int), nd, nd,
]
lib.toggle_traj_.restype = None


def mfpt_python(seed, ttot):
    """Explicit ctypes call to the shared library's MFPT routine."""
    tau = c_double(0.0)
    rate = c_double(0.0)
    ncross = c_int(0)
    lib.toggle_mfpt_(
        byref(c_int(seed)),
        byref(c_double(ttot)), byref(c_double(DT)), byref(c_double(TRELAX)),
        byref(c_double(G0)), byref(c_double(G1)), byref(c_double(KTHR)),
        byref(c_double(KDEG)), byref(c_int(NEXP)), byref(c_double(BNOISE)),
        byref(tau), byref(rate), byref(ncross),
    )
    return tau.value, rate.value, ncross.value


tau_py, rate_py, nc_py = mfpt_python(SEED, TTOT)
print(f"Python ncross        = {nc_py}")
print(f"Python mean FPT tau  = {tau_py:.15e}")
print(f"Python rate kappa    = {rate_py:.15e}")

# ----------------------------------------------------------------------
# R side: call the IDENTICAL compiled routine via .Fortran
# ----------------------------------------------------------------------
r_src = f'''
dyn.load("{LIB}")
res <- .Fortran("toggle_mfpt",
   seed   = as.integer({SEED}),
   ttot   = as.double({TTOT}),
   dt     = as.double({DT}),
   trelax = as.double({TRELAX}),
   g0     = as.double({G0}),
   g1     = as.double({G1}),
   Kthr   = as.double({KTHR}),
   kdeg   = as.double({KDEG}),
   nexp   = as.integer({NEXP}),
   b      = as.double({BNOISE}),
   tau    = as.double(0),
   rate   = as.double(0),
   ncross = as.integer(0))
cat(sprintf("R_NCROSS %d\\n", res$ncross))
cat(sprintf("R_TAU %.15e\\n", res$tau))
cat(sprintf("R_RATE %.15e\\n", res$rate))
'''
with open(RSCRIPT, "w") as fh:
    fh.write(r_src)

rate_R = tau_R = None
nc_R = None
rbin = shutil.which("Rscript")
if rbin is None:
    print("R_STATUS: Rscript not found; skipping the R call.")
else:
    out = subprocess.run([rbin, "--vanilla", RSCRIPT],
                         capture_output=True, text=True)
    for line in out.stdout.splitlines():
        p = line.split()
        if len(p) == 2 and p[0] == "R_NCROSS":
            nc_R = int(p[1])
        elif len(p) == 2 and p[0] == "R_TAU":
            tau_R = float(p[1])
        elif len(p) == 2 and p[0] == "R_RATE":
            rate_R = float(p[1])
    if rate_R is None:
        print("R_STATUS: R call failed; stderr follows:")
        print(out.stderr.strip())

if rate_R is not None:
    print(f"R ncross             = {nc_R}")
    print(f"R mean FPT tau       = {tau_R:.15e}")
    print(f"R rate kappa         = {rate_R:.15e}")

# ----------------------------------------------------------------------
# Cross-check: same seed + same library => identical rates
# ----------------------------------------------------------------------
if rate_R is not None:
    diff = abs(rate_R - rate_py)
    print(f"|rate_R - rate_Python| = {diff:.3e}")
    print(f"identical rates      = {rate_R == rate_py}")

# ----------------------------------------------------------------------
# Accuracy check: a longer total time gives a more accurate rate estimate
# ----------------------------------------------------------------------
tau_long, rate_long, nc_long = mfpt_python(SEED, TTOT_LONG)
print(f"Longer-run total time = {TTOT_LONG}")
print(f"Longer-run ncross     = {nc_long}")
print(f"Longer-run mean FPT   = {tau_long:.15e}")
print(f"Longer-run rate kappa = {rate_long:.15e}")

# ----------------------------------------------------------------------
# Explanation of why the R/Python agreement confirms the result
# ----------------------------------------------------------------------
print("EXPLANATION: Because R and Python only pass the same seed and parameters "
      "into the one compiled routine and get bit-identical rates back, the RNG and "
      "SDE integration must both reside in the shared library, so the measured rate "
      "is a property of that shared numerics rather than of either language wrapper.")

# ----------------------------------------------------------------------
# Illustrative sample trajectory (from the same library) for the figure
# ----------------------------------------------------------------------
nsteps = 5_000_000
stride = 500
nout = nsteps // stride
o1 = np.zeros(nout, dtype=np.float64)
o2 = np.zeros(nout, dtype=np.float64)
lib.toggle_traj_(
    byref(c_int(SEED)), byref(c_int(nsteps)), byref(c_int(stride)),
    byref(c_double(DT)), byref(c_double(G0)), byref(c_double(G1)),
    byref(c_double(KTHR)), byref(c_double(KDEG)), byref(c_int(NEXP)),
    byref(c_double(BNOISE)), byref(c_int(nout)), o1, o2,
)
tgrid = np.arange(1, nout + 1) * stride * DT

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(tgrid, o1, lw=0.7, label="gene 1")
ax.plot(tgrid, o2, lw=0.7, label="gene 2")
ax.set_xlabel("time")
ax.set_ylabel("expression level")
title = f"Toggle switch (seed {SEED}):  Python kappa = {rate_py:.4e}"
if rate_R is not None:
    title += f",  R kappa = {rate_R:.4e}"
ax.set_title(title)
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(PNG)
print(f"figure saved         = {PNG}")
