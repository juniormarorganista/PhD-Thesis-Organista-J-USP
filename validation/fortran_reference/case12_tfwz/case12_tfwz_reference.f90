program case12_tfwz_reference
  !! Standalone Fortran reference exporter for MMS Case 12 tfwz.
  !!
  !! This utility does not modify or depend on the original thesis Fortran files.
  !! It mirrors the original expression from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_12.f
  !!
  !! Usage:
  !!   gfortran -O2 case12_tfwz_reference.f90 -o case12_tfwz_reference.exe
  !!   .\case12_tfwz_reference.exe case12_tfwz_5x5_default.csv 5 5
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer :: i, j, nx, ny
  real(dp) :: x0, x1, y0, y1, dx, dy, x, y, t
  real(dp) :: tfwz
  character(len=512) :: outfile
  character(len=512) :: arg

  nx = 3
  ny = 3
  x0 = 0.0_dp
  x1 = 1.0_dp
  y0 = 0.0_dp
  y1 = 1.0_dp
  t = 0.0_dp
  outfile = 'case12_tfwz_3x3_default.csv'

  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) outfile = trim(arg)

  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg, *) nx

  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg, *) ny

  if (nx <= 1 .or. ny <= 1) stop 'nx and ny must be greater than 1'

  dx = (x1 - x0) / real(nx - 1, dp)
  dy = (y1 - y0) / real(ny - 1, dp)

  open(unit=10, file=trim(outfile), status='replace', action='write')
  write(10,'(A)') 'x,y,tfwz'

  do j = 1, ny
    y = y0 + real(j - 1, dp) * dy
    do i = 1, nx
      x = x0 + real(i - 1, dp) * dx
      tfwz = case12_tfwz_value(x, y, t)
      write(10,'(ES26.17E3,A,ES26.17E3,A,ES26.17E3)') x, ',', y, ',', tfwz
    end do
  end do

  close(10)

contains

  pure real(dp) function case12_tfwz_value(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: Rey = 1.0_dp
    real(dp), parameter :: betann = 0.5_dp
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b2yy = 1.0_dp, b3yy = 1.0_dp

    value = (-(Rey*y) + 2.0_dp*Rey*x*y + &
      Rey*y*cos(2.0_dp*pii*x) - &
      2.0_dp*Rey*x*y*cos(2.0_dp*pii*x) + &
      2.0_dp*pii**2.0_dp*Rey*x*y*cos(2.0_dp*pii*x) - &
      6.0_dp*pii**2.0_dp*Rey*x**2.0_dp*y*cos(2.0_dp*pii*x) + &
      4.0_dp*pii**2.0_dp*Rey*x**3*y*cos(2.0_dp*pii*x) + &
      pii*Rey*y*sin(2.0_dp*pii*x) - &
      2.0_dp*pii*Rey*x*y*sin(2.0_dp*pii*x) + &
      2.0_dp*pii*Rey*x**2.0_dp*y*sin(2.0_dp*pii*x) - &
      2.0_dp*pii**3*Rey*x**2.0_dp*y*sin(2.0_dp*pii*x) + &
      4.0_dp*pii**3*Rey*x**3*y*sin(2.0_dp*pii*x) - &
      2.0_dp*pii**3*Rey*x**4*y*sin(2.0_dp*pii*x) + &
      Rey*y*cos(2.0_dp*pii*y)* &
      (-((-1.0_dp + 2.0_dp*x)* &
      (1.0_dp + 2.0_dp*pii**2.0_dp*(-1.0_dp + x)*x - cos(2.0_dp*pii*x))) + &
      pii*(-1.0_dp - 2.0_dp*(-1.0_dp + x)*x)*sin(2.0_dp*pii*x)) + &
      2.0_dp*exp(at*t)*pii*cos(pii*x)* &
      (-(a3xy*b1xy*(-1.0_dp + betann)*pii*Rey) + &
      pii*((a2xx*b2xx - a2yy*b2yy)*Rey + &
      betann*(-8.0_dp - a2xx*b2xx*Rey + &
      a2yy*b2yy*Rey + 16.0_dp*x))*cos(pii*y) + &
      ((a2xx*b3xx - a2yy*b3yy)*(-1.0_dp + betann)*pii* &
      Rey - 2.0_dp*(4.0_dp*betann*pii**2.0_dp - at*Rey)* &
      (-1.0_dp + 2.0_dp*x)*y)*sin(pii*y)) + &
      exp(at*t)*(-2.0_dp*a2xy*b1xy*(-1.0_dp + betann)*pii**2*Rey* &
      sin(pii*x) + &
      2.0_dp*pii*cos(pii*y)* &
      (a1xy*b3xy*(-1.0_dp + betann)*pii*Rey + &
      (Rey*(-(a3xx*b2xx*pii) + a3yy*b2yy*pii + &
      2.0_dp*at*(-1.0_dp + x)*x) + &
      betann* &
      (8.0_dp + a3xx*b2xx*pii*Rey - &
      pii*(a3yy*b2yy*Rey + 8.0_dp*pii*(-1.0_dp + x)*x) &
      ))*sin(pii*x)) + &
      2.0_dp*(a1xy*b2xy*(-1.0_dp + betann)*pii**2*Rey + &
      (-((a3xx*b3xx - a3yy*b3yy)*(-1.0_dp + betann)* &
      pii**2.0_dp*Rey) + &
      2.0_dp*(at*Rey + &
      pii**2* &
      (-(at*Rey*(-1.0_dp + x)*x) + &
      2.0_dp*betann*(-4.0_dp + pii**2*(-1.0_dp + x)*x))) &
      *y)*sin(pii*x))*sin(pii*y)) + &
      pii*Rey*((-1.0_dp + 2.0_dp*x)* &
      (y**2 + (-1.0_dp + x)*x*(-1.0_dp + 2.0_dp*pii**2*y**2)) + &
      (-1.0_dp + 2.0_dp*x)*((-1.0_dp + x)*x - y**2)*cos(2.0_dp*pii*x) - &
      pii*((-1.0_dp + x)**2*x**2 + &
      (-1.0_dp - 2.0_dp*(-1.0_dp + x)*x)*y**2)*sin(2.0_dp*pii*x))* &
      sin(2.0_dp*pii*y))/(2.0_dp*exp(2.0_dp*at*t)*Rey)
  end function case12_tfwz_value

end program case12_tfwz_reference
