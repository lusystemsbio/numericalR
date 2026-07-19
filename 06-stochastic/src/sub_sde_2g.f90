subroutine random(idum,prandom)
   implicit none
   INTEGER,intent(inout) :: idum
   double precision,intent(out) :: prandom
   integer,parameter :: IA=(16807),IM=(2147483647),IQ=(127773),  &
      IR=(2836),MASK=(123459876)
   double precision :: AM
   INTEGER :: k
   AM=1./IM
   idum=ieor(idum,MASK)
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   prandom=AM*idum
   idum=ieor(idum,MASK)
end subroutine random

subroutine init_rand(r)
   implicit none
   integer,intent(out) :: r
   double precision :: r1

   call random(r,r1)

end subroutine init_rand

subroutine gaussian(r,grandom)
   implicit none
   integer,intent(inout) :: r
   double precision,intent(out) :: grandom

   double precision :: u1, v1, u2, v2, s

   do
      call random(r,u1)
      call random(r,u2)
      v1 = 2 * u1 - 1.d0
      v2 = 2 * u2 - 1.d0
      s = v1**2 + v2**2
      if ((s.gt.0.d0).and.(s .lt. 1.d0))exit
   enddo

   grandom = dsqrt(-2.d0*dlog(s)/s)*v1

end subroutine gaussian

subroutine sde_euler(n,x,f_func,g_func,r,dt)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: dt
   double precision,intent(inout) :: x(n)
   integer,intent(inout) :: r
   external :: f_func, g_func

   double precision :: f(n), g(n), v(n), dt2
   integer :: i

   dt2 = dsqrt(dt)
   call f_func(n,x,f)
   call g_func(n,x,g)
   do i = 1, n
      call gaussian(r,v(i))
   enddo

   x = x + f * dt + g*v * dt2

end subroutine sde_euler

subroutine sde_milstein(n,x,f_func,g_func,dg_func,r,dt)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: dt
   double precision,intent(inout) :: x(n)
   integer,intent(inout) :: r
   external :: f_func, g_func, dg_func

   double precision :: f(n), g(n), dg(n), v(n), dt2
   integer :: i

   dt2 = dsqrt(dt)
   call f_func(n,x,f)
   call g_func(n,x,g)
   call dg_func(n,x,dg)
   do i = 1, n
      call gaussian(r,v(i))
   enddo

   x = x + f * dt + g*v * dt2 + &
      0.5d0* g * dg * (v*v-1.d0) * dt

end subroutine sde_milstein

subroutine sde_rk(n,x,f_func,g_func,r,dt)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: dt
   double precision,intent(inout) :: x(n)
   integer,intent(inout) :: r
   external :: f_func, g_func

   double precision :: f(n), g(n), v(n), dt2, x1(n), g1(n), x2(n)
   integer :: i

   dt2 = dsqrt(dt)
   call f_func(n,x,f)
   call g_func(n,x,g)
   do i = 1, n
      call gaussian(r,v(i))
   enddo

   x2 = x + f * dt !+ g*v * dt2
   x1 = x2 + g * dt2

   call g_func(n,x1,g1)

   x = x2 + g*v*dt2 + 0.5d0 *(g1 - g) * (v*v-1.d0) * dt2

end subroutine sde_rk

function hill_inh(x, x0, n)
   implicit none
   integer, intent(in) :: n
   double precision,intent(in) :: x, x0
   double precision :: hill_inh
   
   double precision :: a
   
   a = (x/x0)**n
   hill_inh = 1/(1+a)
   
end function hill_inh

subroutine f_ts(n, x, f)
   implicit none
   integer, intent(in) :: n
   double precision,intent(in) :: x(n)
   double precision,intent(out) :: f(n)
   double precision :: g0 = (10.d0), g1 = (40.d0), x_hill = (100.d0), k = (0.1d0)
   integer :: n_hill = (4)
   double precision :: hill_inh
   
   f(1) =  g0 + g1*hill_inh(x(2), x_hill, n_hill) - k * x(1)
   f(2) =  g0 + g1*hill_inh(x(1), x_hill, n_hill) - k * x(2)
  
end subroutine f_ts

subroutine g_ts(n, x, g)
   implicit none
   integer, intent(in) :: n
   double precision,intent(in) :: x(n)
   double precision,intent(out) :: g(n)
   
   g(1:2) = 20.d0
  
end subroutine g_ts

subroutine sde_simulation(iseed, t_total, kappa, nrun)
   implicit none
   integer,intent(in) :: iseed
   double precision,intent(in) :: t_total
   double precision,intent(out) :: kappa(2)
   integer, intent(out) :: nrun(2)
   
   double precision :: t, r1, r2, dt, t0, trun, x(2)
   integer :: n, state, r, t_relax
   double precision :: mfpt(2)
   logical :: start 
   external :: f_ts, g_ts
   
   n = 2
   dt = 0.01d0
   t_relax = 100.d0
   r = iseed
   call init_rand(r)  ! initiate random number generator
   
   call random(r, r1)
   call random(r, r2)
   x(:) = (/r1*300.d0, r2*300.d0/)   ! the initial condition
   
   t = 0
   t0 = 0
   start = .false.
   nrun(:) = 0
   mfpt(:) = 0.d0
   do
      call sde_euler(n,x,f_ts,g_ts,r,dt)
      t = t + dt
      if(t > t_total)exit
      if(.not.start) then
        if(t - t0 > t_relax)then
          if(x(1) < x(2))then
            state = 1
          else
            state = 2
          endif
          trun = 0.d0
          start = .true.
        endif
      else
        if((state.eq.1).and.(X(1) > X(2)))then
          nrun(1) = nrun(1) + 1
          mfpt(1) = mfpt(1) + trun
          t0 = t
          start = .false.
          trun = 0.d0
          state = 2
        elseif((state.eq.2).and.(X(1) <= X(2)))then
          nrun(2) = nrun(2) + 1
          mfpt(2) = mfpt(2) + trun
          t0 = t
          start = .false.
          trun = 0.d0
          state = 1
        endif
      endif
      
      trun = trun + dt
   enddo
   
   mfpt(:) = mfpt(:)/nrun(:)
   kappa(:) = 1.d0/2.d0/mfpt(:)
   
end subroutine sde_simulation
