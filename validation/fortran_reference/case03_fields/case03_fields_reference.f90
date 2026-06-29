program case03_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 03 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_03.f
  !!     ::exact_solutions_case_03(t)
  !!     ::source_term_case_03(t)
  !!
  !! Default parameters match the Case03Parameters/Case03SourceTermParameters
  !! Python defaults (Re=1, Wi=1, beta=0.5, epsilon=0, xi=0, alphaG=0, at=0,
  !! all trigonometric a/b coefficients = 1), not the transient Rey/Wi/epsylon
  !! values currently sitting in par.nn.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case03_fields_reference.f90 -o case03_fields_reference.exe
  !!   .\case03_fields_reference.exe case03_fields_5x5_default.csv 5 5
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer :: i, j, nx, ny
  real(dp) :: x0, x1, y0, y1, dx, dy, x, y, t
  character(len=512) :: outfile
  character(len=512) :: arg

  nx = 5
  ny = 5
  x0 = 0.0_dp
  x1 = 1.0_dp
  y0 = 0.0_dp
  y1 = 1.0_dp
  t = 0.0_dp
  outfile = 'case03_fields_5x5_default.csv'

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
  write(10,'(A)') 'x,y,u,v,wz,psi,txx,txy,tyy,tfwz,tfTxx,tfTxy,tfTyy'

  do j = 1, ny
    y = y0 + real(j - 1, dp) * dy
    do i = 1, nx
      x = x0 + real(i - 1, dp) * dx
      write(10,'(*(ES26.17E3,:,","))') x, y, &
        u_case03(x, y, t), v_case03(x, y, t), wz_case03(x, y, t), psi_case03(x, y, t), &
        txx_case03(x, y, t), txy_case03(x, y, t), tyy_case03(x, y, t), &
        tfwz_case03(x, y, t), tftxx_case03(x, y, t), tftxy_case03(x, y, t), tftyy_case03(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (8*(x**2 - 2*x**3 + x**4)* &
(-2*y + 4*y**3))/exp(at*t)
  end function u_case03

  pure real(dp) function v_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (-32*x*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.))/exp(at*t)
  end function v_case03

  pure real(dp) function wz_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (8*(x**2 - 2*x**3 + x**4)* &
(-2 + 12*y**2))/exp(at*t) + (32*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) + &
(32*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.))/exp(at*t)
  end function wz_case03

  pure real(dp) function psi_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (16*(x**2/2. - x**3 + x**4/2.)* &
y**2*(-1 + y**2))/exp(at*t)
  end function psi_case03

  pure real(dp) function txx_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    value = ((1 - betann)*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y)))/exp(at*t)
  end function txx_case03

  pure real(dp) function txy_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    value = ((1 - betann)*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y)))/exp(at*t)
  end function txy_case03

  pure real(dp) function tyy_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = ((1 - betann)*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b1yy + b3yy*cos(pii*y) + &
b2yy*sin(pii*y)))/exp(at*t)
  end function tyy_case03

  pure real(dp) function tfwz_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = (8*at*(x**2 - 2*x**3 + x**4)* &
(-2 + 12*y**2))/exp(at*t) + (32*at*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) + &
(32*at*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/ &
exp(at*t) + (32*x*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.)*((192*(x**2 - 2*x**3 + x**4)* &
y)/exp(at*t) + (32*x*(-3 + 4*x)*(-y + 2*y**3))/exp(at*t) &
+ (32*(1 - 3*x + 2*x**2)*(-y + 2*y**3))/exp(at*t))) &
/exp(at*t) + (betann*((192*(x**2 - 2*x**3 + x**4))/ &
exp(at*t) + (32*x*(-3 + 4*x)*(-1 + 6*y**2))/exp(at*t) &
+ (32*(1 - 3*x + 2*x**2)*(-1 + 6*y**2))/exp(at*t) &
+ (8*(2 - 12*x + 12*x**2)*(-2 + 12*y**2))/exp(at*t) &
+ (384*(-0.5*y**2 + y**4/2.))/exp(at*t)))/Rey - &
(8*(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)* &
((8*(2*x - 6*x**2 + 4*x**3)*(-2 + 12*y**2))/exp(at*t) &
+ (128*x*(-0.5*y**2 + y**4/2.))/exp(at*t) + &
(64*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/exp(at*t)))/exp(at*t) + &
(32*x*(1 - 3*x + 2*x**2)*(-y + 2*y**3)* &
((8*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2))/exp(at*t) &
+ (32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/exp(at*t) + &
(32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/ &
exp(at*t)))/exp(at*t) - (8*(2*x - 6*x**2 + 4*x**3)* &
(-2*y + 4*y**3)*((8*(x**2 - 2*x**3 + x**4)* &
(-2 + 12*y**2))/exp(at*t) + (32*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) + &
(32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/ &
exp(at*t)))/exp(at*t) - ((1 - betann)* &
(-(a3xy*pii**2*cos(pii*x)) - a2xy*pii**2*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/exp(at*t) &
+ ((1 - betann)*(a2xx*pii*cos(pii*x) - a3xx*pii*sin(pii*x))* &
(b2xx*pii*cos(pii*y) - b3xx*pii*sin(pii*y)))/exp(at*t) - &
((1 - betann)*(a2yy*pii*cos(pii*x) - a3yy*pii*sin(pii*x))* &
(b2yy*pii*cos(pii*y) - b3yy*pii*sin(pii*y)))/exp(at*t) + &
((1 - betann)*(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(-(b3xy*pii**2*cos(pii*y)) - b2xy*pii**2*sin(pii*y)))/ &
exp(at*t)
  end function tfwz_case03

  pure real(dp) function tftxx_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = (-16*(1 - betann)*(2*x - 6*x**2 + 4*x**3)* &
(-2*y + 4*y**3))/(exp(at*t)*Rey) + (alphaG*Rey*Wi* &
(((1 - betann)**2*(a1xx+a3xx*cos(pii*x)+a2xx*sin(pii*x))**2* &
(b1xx + b3xx*cos(pii*y) +b2xx*sin(pii*y))**2)/exp(2*at*t) + &
((1 - betann)**2*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))**2*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y))**2)/exp(2*at*t)))/(1 - betann) &
+ Wi*(-((at*(1 - betann)*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y)))/exp(at*t)) - (32*(1 - betann)*x* &
(1 - 3*x + 2*x**2)*(-y + 2*y**3)*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y)))/ &
exp(2*at*t) - (8*(1 - betann)*(2*x - 6*x**2 + 4*x**3)* &
(-2*y + 4*y**3)*(a1xx + a3xx*cos(pii*x) + a2xx*sin(pii*x))* &
(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y)))/exp(2*at*t) + &
(8*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)* &
(a2xx*pii*cos(pii*x) - a3xx*pii*sin(pii*x))* &
(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y)))/exp(2*at*t) - &
(16*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)* &
(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/ &
exp(2*at*t) - (32*(1 - betann)*x*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.)*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b2xx*pii*cos(pii*y) - &
b3xx*pii*sin(pii*y)))/exp(2*at*t) + xi*((16*(1 - betann)* &
(2*x - 6*x**2 + 4*x**3)*(-2*y + 4*y**3)* &
(a1xx + a3xx*cos(pii*x) + a2xx*sin(pii*x))* &
(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y)))/exp(2*at*t) + &
(8*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)* &
(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))*(b1xy + &
b3xy*cos(pii*y) + b2xy*sin(pii*y)))/exp(2*at*t) + &
((1 - betann)*((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.) &
)/exp(at*t) - (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.) &
)/exp(at*t))*(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/exp(at*t))) + &
((1 - betann)*(a1xx + a3xx*cos(pii*x) + a2xx*sin(pii*x))* &
(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y))* &
(1 + (epsylon*Rey*Wi*(((1 - betann)*(a1xx + &
a3xx*cos(pii*x) + a2xx*sin(pii*x))*(b1xx + &
b3xx*cos(pii*y) + b2xx*sin(pii*y)))/exp(at*t) + &
((1 - betann)*(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/exp(at*t)))/ &
(1 - betann)))/exp(at*t)
  end function tftxx_case03

  pure real(dp) function tftxy_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = -(((1 - betann)*((8*(x**2 - 2*x**3 + x**4)* &
(-2 + 12*y**2))/exp(at*t) - (32*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) - (32*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.))/exp(at*t)))/Rey) + (alphaG*Rey*Wi* &
(((1 - betann)**2*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y)))/exp(2*at*t) + &
((1 - betann)**2*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y))*(b1yy + b3yy*cos(pii*y) + &
b2yy*sin(pii*y)))/exp(2*at*t)))/(1 - betann) &
+ ((1 - betann)*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y))*(1 + (epsylon*Rey*Wi* &
(((1 - betann)*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y)))/exp(at*t) + ((1 - betann)* &
(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/ &
exp(at*t)))/(1 - betann)))/exp(at*t) + &
Wi*(-(((1 - betann)*((-32*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) - &
(32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.) &
)/exp(at*t))*(a1xx + a3xx*cos(pii*x) + &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y)))/exp(at*t)) - &
(at*(1 - betann)*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y)))/exp(at*t) + (8*(1 - betann)* &
(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)* &
(a2xy*pii*cos(pii*x) - a3xy*pii*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/ &
exp(2*at*t) - (8*(1 - betann)*(x**2 - 2*x**3 + x**4)* &
(-2 + 12*y**2)*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b1yy + b3yy*cos(pii*y) + &
b2yy*sin(pii*y)))/exp(2*at*t) - (32*(1 - betann)*x* &
(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.)* &
(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(b2xy*pii*cos(pii*y) - b3xy*pii*sin(pii*y)))/ &
exp(2*at*t) + xi*((4*(1 - betann)* &
(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)* &
(a1xx + a3xx*cos(pii*x) + a2xx*sin(pii*x))* &
(b1xx + b3xx*cos(pii*y) + b2xx*sin(pii*y)))/ &
exp(2*at*t) + ((1 - betann)*((-32*x*(-3 + 4*x)* &
(-0.5*y**2 + y**4/2.))/exp(at*t) - (32*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.))/exp(at*t))*(a1xx+a3xx*cos(pii*x)+ &
a2xx*sin(pii*x))*(b1xx + b3xx*cos(pii*y) + &
b2xx*sin(pii*y)))/(2.*exp(at*t)) - (32*(1 - betann)*x* &
(1 - 3*x + 2*x**2)*(-y + 2*y**3)*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y)))/exp(2*at*t) + (8*(1 - betann)* &
(2*x - 6*x**2 + 4*x**3)*(-2*y + 4*y**3)*(a1xy + &
a3xy*cos(pii*x) + a2xy*sin(pii*x))*(b1xy + &
b3xy*cos(pii*y) + b2xy*sin(pii*y)))/exp(2*at*t) + &
(4*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)* &
(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))*(b1yy + &
b3yy*cos(pii*y) + b2yy*sin(pii*y)))/exp(2*at*t) + &
((1 - betann)*((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.) &
)/exp(at*t) - (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.) &
)/exp(at*t))*(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/ &
(2.*exp(at*t))))
  end function tftxy_case03

  pure real(dp) function tftyy_case03(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = (64*(1 - betann)*x*(1 - 3*x + 2*x**2)* &
(-y + 2*y**3))/(exp(at*t)*Rey) + (alphaG*Rey*Wi* &
(((1 - betann)**2*(a1xy + a3xy*cos(pii*x) + &
a2xy*sin(pii*x))**2*(b1xy + b3xy*cos(pii*y) + &
b2xy*sin(pii*y))**2)/exp(2*at*t) + ((1 - betann)**2* &
(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))**2* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y))**2)/ &
exp(2*at*t)))/(1 - betann) + ((1 - betann)* &
(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y))* &
(1 + (epsylon*Rey*Wi*(((1 - betann)* &
(a1xx + a3xx*cos(pii*x) + a2xx*sin(pii*x))*(b1xx + &
b3xx*cos(pii*y) + b2xx*sin(pii*y)))/exp(at*t) + &
((1 - betann)*(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/exp(at*t)))/ &
(1 - betann)))/exp(at*t) + Wi*((-2*(1 - betann)* &
((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/exp(at*t) - &
(32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/exp(at*t))* &
(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/exp(at*t) - &
(at*(1 - betann)*(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/exp(at*t) + &
(32*(1 - betann)*x*(1 - 3*x + 2*x**2)*(-y + 2*y**3)* &
(a1yy + a3yy*cos(pii*x) + a2yy*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/ &
exp(2*at*t) + (8*(1 - betann)*(2*x - 6*x**2 + 4*x**3)* &
(-2*y + 4*y**3)*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b1yy + b3yy*cos(pii*y) + &
b2yy*sin(pii*y)))/exp(2*at*t) + (8*(1 - betann)* &
(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)* &
(a2yy*pii*cos(pii*x) - a3yy*pii*sin(pii*x))* &
(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/ &
exp(2*at*t) - (32*(1 - betann)*x*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.)*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b2yy*pii*cos(pii*y) - &
b3yy*pii*sin(pii*y)))/exp(2*at*t) + xi*((8*(1 - betann)* &
(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*(a1xy + &
a3xy*cos(pii*x) + a2xy*sin(pii*x))*(b1xy + &
b3xy*cos(pii*y) + b2xy*sin(pii*y)))/ &
exp(2*at*t) + ((1 - betann)* &
((-32*x*(-3 + 4*x)* (-0.5*y**2 + y**4/2.) &
)/exp(at*t) - (32*(1 - 3*x + 2*x**2)* &
(-0.5*y**2 + y**4/2.))/exp(at*t))* &
(a1xy + a3xy*cos(pii*x) + a2xy*sin(pii*x))* &
(b1xy + b3xy*cos(pii*y) + b2xy*sin(pii*y)))/ &
exp(at*t) - (64*(1 - betann)*x*(1 - 3*x + 2*x**2)* &
(-y + 2*y**3)*(a1yy + a3yy*cos(pii*x) + &
a2yy*sin(pii*x))*(b1yy + b3yy*cos(pii*y) + b2yy*sin(pii*y)))/ &
exp(2*at*t)))
  end function tftyy_case03

end program case03_fields_reference
