program case05_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 05 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_05.f
  !!     ::exact_solutions_case_05(t)
  !!     ::source_term_case_05(t)
  !!
  !! Default parameters match the Case05Parameters/Case05SourceTermParameters
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
  outfile = 'case05_fields_5x5_default.csv'

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
        u_case05(x, y, t), v_case05(x, y, t), wz_case05(x, y, t), psi_case05(x, y, t), &
        txx_case05(x, y, t), txy_case05(x, y, t), tyy_case05(x, y, t), &
        tfwz_case05(x, y, t), tftxx_case05(x, y, t), tftxy_case05(x, y, t), tftyy_case05(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (8*cos(pii*x)*sin(pii*y))/exp(at*t)
  end function u_case05

  pure real(dp) function v_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (-8*cos(pii*y)*sin(pii*x))/exp(at*t)
  end function v_case05

  pure real(dp) function wz_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (16*pii*cos(pii*x)*cos(pii*y))/exp(at*t)
  end function wz_case05

  pure real(dp) function psi_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (-8*cos(pii*x)*cos(pii*y))/ &
(exp(at*t)*pii)
  end function psi_case05

  pure real(dp) function txx_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    value = ((1 - betann)* &
(a1xx + a2xx*x + a3xx*x**2)* &
(b1xx + b2xx*y + b3xx*y**2))/ &
exp(at*t)
  end function txx_case05

  pure real(dp) function txy_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    value = ((1 - betann)* &
(a1xy + a2xy*x + a3xy*x**2)* &
(b1xy + b2xy*y + b3xy*y**2))/ &
exp(at*t)
  end function txy_case05

  pure real(dp) function tyy_case05(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = ((1 - betann)* &
(a1yy + a2yy*x + a3yy*x**2)* &
(b1yy + b2yy*y + b3yy*y**2))/ &
exp(at*t)
  end function tyy_case05

  pure real(dp) function tfwz_case05(x, y, t) result(value)
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
    value = (2*b3xy*(1 - betann)* &
(a1xy + a2xy*x + a3xy*x**2))/ &
exp(at*t) + &
((1 - betann)*(a2xx + 2*a3xx*x)* &
(b2xx + 2*b3xx*y))/exp(at*t) - &
((1 - betann)*(a2yy + 2*a3yy*x)* &
(b2yy + 2*b3yy*y))/exp(at*t) - &
(2*a3xy*(1 - betann)* &
(b1xy + b2xy*y + b3xy*y**2))/ &
exp(at*t) + &
(16*at*pii*cos(pii*x)*cos(pii*y))/ &
exp(at*t) - &
(32*betann*pii**3*cos(pii*x)* &
cos(pii*y))/(exp(at*t)*Rey)
  end function tfwz_case05

  pure real(dp) function tftxx_case05(x, y, t) result(value)
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
    value = (alphaG*Rey*Wi* &
(((1 - betann)**2* &
(a1xx + a2xx*x + &
a3xx*x**2)**2* &
(b1xx + b2xx*y + &
b3xx*y**2)**2)/ &
exp(2*at*t) + &
((1 - betann)**2* &
(a1xy + a2xy*x + &
a3xy*x**2)**2* &
(b1xy + b2xy*y + &
b3xy*y**2)**2)/ &
exp(2*at*t)))/(1 - betann) &
+ ((1 - betann)* &
(a1xx + a2xx*x + a3xx*x**2)* &
(b1xx + b2xx*y + b3xx*y**2)* &
(1 + (epsylon*Rey*Wi* &
(((1 - betann)* &
(a1xx + a2xx*x + &
a3xx*x**2)* &
(b1xx + b2xx*y + &
b3xx*y**2))/exp(at*t) &
+ ((1 - betann)* &
(a1yy + a2yy*x + &
a3yy*x**2)* &
(b1yy + b2yy*y + &
b3yy*y**2))/exp(at*t) &
))/(1 - betann)))/ &
exp(at*t) + &
(16*(1 - betann)*pii*sin(pii*x)* &
sin(pii*y))/(exp(at*t)*Rey) + &
Wi*(-((at*(1 - betann)* &
(a1xx + a2xx*x + &
a3xx*x**2)* &
(b1xx + b2xx*y + &
b3xx*y**2))/exp(at*t)) &
- (16*(1 - betann)*pii* &
(a1xy + a2xy*x + a3xy*x**2)* &
(b1xy + b2xy*y + b3xy*y**2)* &
cos(pii*x)*cos(pii*y))/ &
exp(2*at*t) - &
(8*(1 - betann)* &
(a1xx + a2xx*x + a3xx*x**2)* &
(b2xx + 2*b3xx*y)* &
cos(pii*y)*sin(pii*x))/ &
exp(2*at*t) + &
(8*(1 - betann)* &
(a2xx + 2*a3xx*x)* &
(b1xx + b2xx*y + b3xx*y**2)* &
cos(pii*x)*sin(pii*y))/ &
exp(2*at*t) + &
(16*(1 - betann)*pii* &
(a1xx + a2xx*x + a3xx*x**2)* &
(b1xx + b2xx*y + b3xx*y**2)* &
sin(pii*x)*sin(pii*y))/ &
exp(2*at*t) - &
(16*(1 - betann)*pii* &
(a1xx + a2xx*x + a3xx*x**2)* &
xi* &
(b1xx + b2xx*y + b3xx*y**2)* &
sin(pii*x)*sin(pii*y))/ &
exp(2*at*t))
  end function tftxx_case05

  pure real(dp) function tftxy_case05(x, y, t) result(value)
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
    value = (alphaG*Rey*Wi* &
(((1 - betann)**2* &
(a1xx + a2xx*x + &
a3xx*x**2)* &
(a1xy + a2xy*x + &
a3xy*x**2)* &
(b1xx + b2xx*y + &
b3xx*y**2)* &
(b1xy + b2xy*y + &
b3xy*y**2))/exp(2*at*t) &
+ ((1 - betann)**2* &
(a1xy + a2xy*x + &
a3xy*x**2)* &
(a1yy + a2yy*x + &
a3yy*x**2)* &
(b1xy + b2xy*y + &
b3xy*y**2)* &
(b1yy + b2yy*y + &
b3yy*y**2))/exp(2*at*t)) &
)/(1 - betann) + &
((1 - betann)* &
(a1xy + a2xy*x + a3xy*x**2)* &
(b1xy + b2xy*y + b3xy*y**2)* &
(1 + (epsylon*Rey*Wi* &
(((1 - betann)* &
(a1xx + a2xx*x + &
a3xx*x**2)* &
(b1xx + b2xx*y + &
b3xx*y**2))/exp(at*t) &
+ ((1 - betann)* &
(a1yy + a2yy*x + &
a3yy*x**2)* &
(b1yy + b2yy*y + &
b3yy*y**2))/exp(at*t) &
))/(1 - betann)))/ &
exp(at*t) + &
Wi*(-((at*(1 - betann)* &
(a1xy + a2xy*x + &
a3xy*x**2)* &
(b1xy + b2xy*y + &
b3xy*y**2))/exp(at*t)) &
+ (8*(1 - betann)*pii* &
(a1xx + a2xx*x + a3xx*x**2)* &
(b1xx + b2xx*y + b3xx*y**2)* &
cos(pii*x)*cos(pii*y))/ &
exp(2*at*t) - &
(8*(1 - betann)*pii* &
(a1yy + a2yy*x + a3yy*x**2)* &
(b1yy + b2yy*y + b3yy*y**2)* &
cos(pii*x)*cos(pii*y))/ &
exp(2*at*t) - &
(8*(1 - betann)* &
(a1xy + a2xy*x + a3xy*x**2)* &
(b2xy + 2*b3xy*y)* &
cos(pii*y)*sin(pii*x))/ &
exp(2*at*t) + &
(8*(1 - betann)* &
(a2xy + 2*a3xy*x)* &
(b1xy + b2xy*y + b3xy*y**2)* &
cos(pii*x)*sin(pii*y))/ &
exp(2*at*t))
  end function tftxy_case05

  pure real(dp) function tftyy_case05(x, y, t) result(value)
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
    value = (alphaG*Rey*Wi*(((1 - betann)**2* &
(a1xy + a2xy*x + a3xy*x**2)**2*(b1xy + b2xy*y + b3xy*y**2)**2)/ &
exp(2*at*t) + ((1 - betann)**2*(a1yy + a2yy*x + a3yy*x**2)**2* &
(b1yy + b2yy*y + b3yy*y**2)**2)/exp(2*at*t)))/(1 - betann) &
+ ((1 - betann)*(a1yy + a2yy*x + a3yy*x**2)* &
(b1yy + b2yy*y + b3yy*y**2)*(1 + (epsylon*Rey*Wi*(((1 - betann)* &
(a1xx + a2xx*x + a3xx*x**2)*(b1xx + b2xx*y + &
b3xx*y**2))/exp(at*t) &
+ ((1 - betann)* &
(a1yy + a2yy*x + &
a3yy*x**2)* &
(b1yy + b2yy*y + &
b3yy*y**2))/exp(at*t) &
))/(1 - betann)))/ &
exp(at*t) - &
(16*(1 - betann)*pii*sin(pii*x)* &
sin(pii*y))/(exp(at*t)*Rey) + &
Wi*(-((at*(1 - betann)* &
(a1yy + a2yy*x + &
a3yy*x**2)* &
(b1yy + b2yy*y + &
b3yy*y**2))/exp(at*t)) &
+ (16*(1 - betann)*pii* &
(a1xy + a2xy*x + a3xy*x**2)* &
(b1xy + b2xy*y + b3xy*y**2)* &
cos(pii*x)*cos(pii*y))/ &
exp(2*at*t) - &
(8*(1 - betann)* &
(a1yy + a2yy*x + a3yy*x**2)* &
(b2yy + 2*b3yy*y)* &
cos(pii*y)*sin(pii*x))/ &
exp(2*at*t) + &
(8*(1 - betann)* &
(a2yy + 2*a3yy*x)* &
(b1yy + b2yy*y + b3yy*y**2)* &
cos(pii*x)*sin(pii*y))/ &
exp(2*at*t) - &
(16*(1 - betann)*pii* &
(a1yy + a2yy*x + a3yy*x**2)* &
(b1yy + b2yy*y + b3yy*y**2)* &
sin(pii*x)*sin(pii*y))/ &
exp(2*at*t) + &
(16*(1 - betann)*pii* &
(a1yy + a2yy*x + a3yy*x**2)* &
xi* &
(b1yy + b2yy*y + b3yy*y**2)* &
sin(pii*x)*sin(pii*y))/ &
exp(2*at*t))
  end function tftyy_case05

end program case05_fields_reference
