program case11_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 11 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_11.f
  !!     ::exact_solutions_case_11(t)
  !!     ::source_term_case_11(t)
  !!
  !! Default parameters match the Case11Parameters/Case11SourceTermParameters
  !! Python defaults (Re=1, Wi=1, beta=0.5, epsilon=0, xi=0, alphaG=0, at=0,
  !! all polynomial a/b coefficients = 1), not the transient Rey/Wi/epsylon
  !! values currently sitting in par.nn.
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
  outfile = 'case11_fields_5x5_default.csv'

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
        u_case11(x, y, t), v_case11(x, y, t), wz_case11(x, y, t), psi_case11(x, y, t), &
        txx_case11(x, y, t), txy_case11(x, y, t), tyy_case11(x, y, t), &
        tfwz_case11(x, y, t), tftxx_case11(x, y, t), tftxy_case11(x, y, t), tftyy_case11(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1.0d0 + x)*x*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)))/exp(at*t)
  end function u_case11

  pure real(dp) function v_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = -((y*(pii*(-1.0d0 + x)*x*cos(pii*x) + &
(-1.0d0 + 2.0d0*x)*sin(pii*x))*sin(pii*y))/exp(at*t))
  end function v_case11

  pure real(dp) function wz_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (2.0d0*(pii*(-1.0d0 + x)*x*cos(pii*y)* &
sin(pii*x) + &
y*(pii*(-1.0d0 + 2.0d0*x)*cos(pii*x) + &
(1.0d0 - pii**2*(-1.0d0 + x)*x)*sin(pii*x))* &
sin(pii*y)))/exp(at*t)
  end function wz_case11

  pure real(dp) function psi_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1.0d0 + x)*x*y*sin(pii*x)* &
sin(pii*y))/exp(at*t)
  end function psi_case11

  pure real(dp) function txx_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    value = -(((-1.0d0 + betann)* &
(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y)))/ &
exp(at*t))
  end function txx_case11

  pure real(dp) function txy_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    value = -(((-1.0d0 + betann)* &
(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y)))/ &
exp(at*t))
  end function txy_case11

  pure real(dp) function tyy_case11(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = -(((-1.0d0 + betann)* &
(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y)))/ &
exp(at*t))
  end function tyy_case11

  pure real(dp) function tfwz_case11(x, y, t) result(value)
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
    value = (Rey*(-1 + 2*x)*y + &
2*(-1 + betann)*exp(at*t)*Rey* &
(a2yy*b2yy - 2*a1xy*b3xy + 2*a2yy*b3yy*y - &
a2xx*(b2xx + 2*b3xx*y) + &
2*a3xy*(b1xy - b3xy*x**2 + b2xy*y + &
b3xy*y**2) - &
2*x*(a2xy*b3xy + a3xx*(b2xx + 2*b3xx*y) - &
a3yy*(b2yy + 2*b3yy*y))) + &
4*exp(at*t)*pii*(-1 + 2*x)*cos(pii*x)* &
(4*betann*pii*cos(pii*y) + &
(-4*betann*pii**2 + at*Rey)*y*sin(pii*y)) + &
4*exp(at*t)*sin(pii*x)* &
(pii*(at*Rey*(-1 + x)*x + &
betann*(4 - 4*pii**2*(-1 + x)*x))*cos(pii*y) &
+ (at*Rey + &
pii**2*(-(at*Rey*(-1 + x)*x) + &
2*betann*(-4 + pii**2*(-1 + x)*x)))*y* &
sin(pii*y)) + &
Rey*(-1 + 2*x)*cos(2*pii*x)* &
(y*(-1 + 2*pii**2*(-1 + x)*x + cos(2*pii*y)) - &
pii*(x - x**2 + y**2)*sin(2*pii*y)) + &
pii*Rey*sin(2*pii*x)* &
(y - 2*(-1 + x)*x*(-1 + pii**2*(-1 + x)*x)*y + &
(-1 - 2*(-1 + x)*x)*y*cos(2*pii*y) - &
pii*((-1 + x)**2*x**2 + &
(-1 - 2*(-1 + x)*x)*y**2)*sin(2*pii*y)) + &
Rey*(-1 + 2*x)* &
((-1 - 2*pii**2*(-1 + x)*x)*y*cos(2*pii*y) + &
pii*(y**2 + (-1 + x)*x*(-1 + 2*pii**2*y**2))* &
sin(2*pii*y)))/(2.*exp(2*at*t)*Rey)
  end function tfwz_case11

  pure real(dp) function tftxx_case11(x, y, t) result(value)
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
    value = (alphaG*(1 - betann)*Rey*Wi* &
((a1xx + x*(a2xx + a3xx*x))**2* &
(b1xx + y*(b2xx + b3xx*y))**2 + &
(a1xy + x*(a2xy + a3xy*x))**2* &
(b1xy + y*(b2xy + b3xy*y))**2) + &
(1 - betann)*(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))* &
(exp(at*t) + epsylon*Rey*Wi* &
(a1xx*(b1xx + y*(b2xx + b3xx*y)) + &
a1yy*(b1yy + y*(b2yy + b3yy*y)) + &
x*(a2yy*b1yy + a2yy*y*(b2yy + b3yy*y) + &
a2xx*(b1xx + y*(b2xx + b3xx*y)) + &
a3xx*x*(b1xx + y*(b2xx + b3xx*y)) + &
a3yy*x*(b1yy + y*(b2yy + b3yy*y))))) + &
(2*(-1 + betann)*exp(at*t)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)))/Rey + &
Wi*(at*(-1 + betann)*exp(at*t)* &
(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y)) + &
(-1 + betann)*pii*(a1xx + x*(a2xx + a3xx*x))*y* &
(b1xx + y*(b2xx + b3xx*y))*cos(pii*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x)) + &
(-1 + betann)*(a1xx + x*(a2xx + a3xx*x))*y* &
(b2xx + 2*b3xx*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) + &
(-1 + betann)*(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) - &
(-1 + betann)*pii*(-1 + x)*x* &
(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))*cos(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*(-1 + x)*x*(a2xx + 2*a3xx*x)* &
(b1xx + y*(b2xx + b3xx*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*(-1 + x)* &
(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*x*(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) + &
2*(-1 + betann)*(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
2*(-1 + betann)*pii*(-1 + x)*x* &
(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(-2*cos(pii*y) + pii*y*sin(pii*y)) + &
(1 - betann)*xi* &
((a1xy + x*(a2xy + a3xy*x))*y* &
(b1xy + y*(b2xy + b3xy*y))* &
(2*pii*(1 - 2*x)*cos(pii*x) + &
(-2 + pii**2*(-1 + x)*x)*sin(pii*x))* &
sin(pii*y) + &
2*(a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
pii*(-1 + x)*x*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(-2*cos(pii*y) + pii*y*sin(pii*y)))))/ &
exp(2*at*t)
  end function tftxx_case11

  pure real(dp) function tftxy_case11(x, y, t) result(value)
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
    value = (alphaG*(1 - betann)*Rey*Wi* &
(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))* &
((a1xx + x*(a2xx + a3xx*x))* &
(b1xx + y*(b2xx + b3xx*y)) + &
(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))) - &
(-1 + betann)*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))* &
(exp(at*t) + epsylon*Rey*Wi* &
(a1xx*(b1xx + y*(b2xx + b3xx*y)) + &
a1yy*(b1yy + y*(b2yy + b3yy*y)) + &
x*(a2yy*b1yy + a2yy*y*(b2yy + b3yy*y) + &
a2xx*(b1xx + y*(b2xx + b3xx*y)) + &
a3xx*x*(b1xx + y*(b2xx + b3xx*y)) + &
a3yy*x*(b1yy + y*(b2yy + b3yy*y))))) + &
(2*(-1 + betann)*exp(at*t)* &
(pii*(-1 + x)*x*cos(pii*y)*sin(pii*x) - &
y*(pii*(-1 + 2*x)*cos(pii*x) + sin(pii*x))* &
sin(pii*y)))/Rey + &
(-1 + betann)*Wi* &
(at*exp(at*t)*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y)) + &
pii*(a1xy + x*(a2xy + a3xy*x))*y* &
(b1xy + y*(b2xy + b3xy*y))*cos(pii*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x)) + &
(a1xy + x*(a2xy + a3xy*x))*y*(b2xy + 2*b3xy*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) + &
(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) + &
(a1xx + x*(a2xx + a3xx*x))*y* &
(b1xx + y*(b2xx + b3xx*y))* &
(2*pii*(1 - 2*x)*cos(pii*x) + &
(-2 + pii**2*(-1 + x)*x)*sin(pii*x))* &
sin(pii*y) - &
pii*(-1 + x)*x*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*cos(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + x)*x*(a2xy + 2*a3xy*x)* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + x)*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
x*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
pii*(-1 + x)*x*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))*sin(pii*x)* &
(-2*cos(pii*y) + pii*y*sin(pii*y)) + &
xi*(a1xx*(b1xx + y*(b2xx + b3xx*y)) + &
a1yy*(b1yy + y*(b2yy + b3yy*y)) + &
x*(a2yy*b1yy + a2yy*y*(b2yy + b3yy*y) + &
a2xx*(b1xx + y*(b2xx + b3xx*y)) + &
a3xx*x*(b1xx + y*(b2xx + b3xx*y)) + &
a3yy*x*(b1yy + y*(b2yy + b3yy*y))))* &
(-(pii*(-1 + x)*x*cos(pii*y)*sin(pii*x)) + &
y*(pii*(-1 + 2*x)*cos(pii*x) + sin(pii*x))* &
sin(pii*y))))/exp(2*at*t)
  end function tftxy_case11

  pure real(dp) function tftyy_case11(x, y, t) result(value)
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
    value = (alphaG*(1 - betann)*Rey*Wi* &
((a1xy + x*(a2xy + a3xy*x))**2* &
(b1xy + y*(b2xy + b3xy*y))**2 + &
(a1yy + x*(a2yy + a3yy*x))**2* &
(b1yy + y*(b2yy + b3yy*y))**2) + &
(1 - betann)*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))* &
(exp(at*t) + epsylon*Rey*Wi* &
(a1xx*(b1xx + y*(b2xx + b3xx*y)) + &
a1yy*(b1yy + y*(b2yy + b3yy*y)) + &
x*(a2yy*b1yy + a2yy*y*(b2yy + b3yy*y) + &
a2xx*(b1xx + y*(b2xx + b3xx*y)) + &
a3xx*x*(b1xx + y*(b2xx + b3xx*y)) + &
a3yy*x*(b1yy + y*(b2yy + b3yy*y))))) - &
(2*(-1 + betann)*exp(at*t)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)))/Rey + &
Wi*(at*(-1 + betann)*exp(at*t)* &
(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y)) + &
(-1 + betann)*pii*(a1yy + x*(a2yy + a3yy*x))*y* &
(b1yy + y*(b2yy + b3yy*y))*cos(pii*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x)) + &
(-1 + betann)*(a1yy + x*(a2yy + a3yy*x))*y* &
(b2yy + 2*b3yy*y)* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) + &
(-1 + betann)*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))*sin(pii*y) + &
2*(-1 + betann)*(a1xy + x*(a2xy + a3xy*x))*y* &
(b1xy + y*(b2xy + b3xy*y))* &
(2*pii*(1 - 2*x)*cos(pii*x) + &
(-2 + pii**2*(-1 + x)*x)*sin(pii*x))* &
sin(pii*y) - &
(-1 + betann)*pii*(-1 + x)*x* &
(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))*cos(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*(-1 + x)*x*(a2yy + 2*a3yy*x)* &
(b1yy + y*(b2yy + b3yy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*(-1 + x)* &
(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
(-1 + betann)*x*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
2*(-1 + betann)*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)) + &
(1 - betann)*xi* &
((a1xy + x*(a2xy + a3xy*x))*y* &
(b1xy + y*(b2xy + b3xy*y))* &
(2*pii*(1 - 2*x)*cos(pii*x) + &
(-2 + pii**2*(-1 + x)*x)*sin(pii*x))* &
sin(pii*y) - &
2*(a1yy + x*(a2yy + a3yy*x))* &
(b1yy + y*(b2yy + b3yy*y))* &
(pii*(-1 + x)*x*cos(pii*x) + &
(-1 + 2*x)*sin(pii*x))* &
(pii*y*cos(pii*y) + sin(pii*y)) - &
pii*(-1 + x)*x*(a1xy + x*(a2xy + a3xy*x))* &
(b1xy + y*(b2xy + b3xy*y))*sin(pii*x)* &
(-2*cos(pii*y) + pii*y*sin(pii*y)))))/ &
exp(2*at*t)
  end function tftyy_case11

end program case11_fields_reference
