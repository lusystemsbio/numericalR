
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
