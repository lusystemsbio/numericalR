subroutine cal_sd(n,x,sd)
  implicit none
  integer,intent(in) :: n
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: sd
  integer :: i
  double precision :: mean, mean2
  mean = 0.d0
  mean2 = 0.d0
  do i = 1, n
    mean = mean + x(i)
    mean2 = mean2 + x(i)**2
  end do
  mean = mean/n
  mean2 = mean2/n
  sd = sqrt(mean2 - mean**2)
end subroutine cal_sd