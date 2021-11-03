function derivs_1g(t, x)
  implicit none
  double precision :: derivs_1g
  double precision,intent(in) :: t, x
  double precision ::  g = (50.d0), k = (0.1d0)
  
  derivs_1g = g - k * x
  
end function derivs_1g

subroutine ode_1g_euler(x0,dt,nsteps,results)
  implicit none
  integer, intent(in) :: nsteps
  double precision, intent(in) :: x0, dt
  double precision, dimension(2, nsteps), intent(out) :: results
  double precision :: derivs_1g
  integer :: i, n
  double precision ::  t, x
  
  t = 0.d0
  x = x0
  do i = 1, nsteps
    t = t + dt
    x = x + derivs_1g(t, x) * dt
    results(1,i) = t
    results(2,i) = x
  enddo
  
end subroutine ode_1g_euler

subroutine ode_1g_rk2(x0,dt,nsteps,results)
  implicit none
  integer, intent(in) :: nsteps
  double precision, intent(in) :: x0, dt
  double precision, dimension(2, nsteps), intent(out) :: results
  double precision :: derivs_1g
  integer :: i, n
  double precision ::  t, x, k1, k2
  
  t = 0.d0
  x = x0
  do i = 1, nsteps
    k1 = derivs_1g(t,x) * dt
    t = t + dt
    k2 = derivs_1g(t,x+k1) * dt
    x = x + (k1+k2)/2.d0
    results(1,i) = t
    results(2,i) = x
  enddo
end subroutine ode_1g_rk2

subroutine ode_1g_rk4(x0,dt,nsteps,results)
  implicit none
  integer, intent(in) :: nsteps
  double precision, intent(in) :: x0, dt
  double precision, dimension(2, nsteps), intent(out) :: results
  double precision :: derivs_1g
  integer :: i, n
  double precision ::  t, x, k1, k2, k3, k4, dt05
  
  dt05 = dt/2.d0
  t = 0.d0
  x = x0
  do i = 1, nsteps
    k1 = derivs_1g(t,x) * dt
    t = t + dt05
    k2 = derivs_1g(t,x+k1/2.d0) * dt
    k3 = derivs_1g(t,x+k2/2.d0) * dt
    t = t + dt05
    k4 = derivs_1g(t,x+k3) * dt
    x = x + (k1+2.d0*k2+2.d0*k3+k4)/6.d0
    results(1,i) = t
    results(2,i) = x
  enddo
end subroutine ode_1g_rk4