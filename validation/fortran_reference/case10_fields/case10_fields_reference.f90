program case10_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 10 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_10.f
  !!     ::exact_solutions_case_10(t)
  !!     ::source_term_case_10(t)
  !!
  !! Default parameters match the Case10Parameters/Case10SourceTermParameters
  !! Python defaults (Re=1, Wi=1, beta=0.5, epsilon=0, xi=0, alphaG=0, at=0,
  !! c1=1, c2=2, c3=3), not the transient Rey/Wi/epsylon values currently
  !! sitting in par.nn.
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
  outfile = 'case10_fields_5x5_default.csv'

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
        u_case10(x, y, t), v_case10(x, y, t), wz_case10(x, y, t), psi_case10(x, y, t), &
        txx_case10(x, y, t), txy_case10(x, y, t), tyy_case10(x, y, t), &
        tfwz_case10(x, y, t), tftxx_case10(x, y, t), tftxy_case10(x, y, t), tftyy_case10(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1.0d0 + x)*x*sin(pii*x)* &
(pii*y*cos(pii*y) + sin(pii*y)))/exp(at*t)
  end function u_case10

  pure real(dp) function v_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = -((y*(pii*(-1.0d0 + x)*x*cos(pii*x) + &
(-1.0d0 + 2.0d0*x)*sin(pii*x))*sin(pii*y))/exp(at*t))
  end function v_case10

  pure real(dp) function wz_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (2.0d0*(pii*(-1.0d0 + x)*x*cos(pii*y)* &
sin(pii*x) + &
y*(pii*(-1.0d0 + 2.0d0*x)*cos(pii*x) + &
(1.0d0 - pii**2*(-1.0d0 + x)*x)*sin(pii*x))* &
sin(pii*y)))/exp(at*t)
  end function wz_case10

  pure real(dp) function psi_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1.0d0 + x)*x*y*sin(pii*x)* &
sin(pii*y))/exp(at*t)
  end function psi_case10

  pure real(dp) function txx_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c1 = 1.0_dp
    value = -(((-1.0d0 + betann)*c1)/exp(at*t))
  end function txx_case10

  pure real(dp) function txy_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c2 = 2.0_dp
    value = -(((-1.0d0 + betann)*c2)/exp(at*t))
  end function txy_case10

  pure real(dp) function tyy_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c3 = 3.0_dp
    value = -(((-1.0d0 + betann)*c3)/exp(at*t))
  end function tyy_case10

  pure real(dp) function tfwz_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = (2*pii**2*Rey*(-1 + x)*x*(-1 + 2*x)* &
y*cos(pii*x)**2* &
sin(pii*y)*(pii*y*cos(pii*y) + sin(pii*y)) + &
2*pii*cos(pii*x)* &
(pii**2*Rey*x**2*y* &
(-1 + x - x**2 + x*cos(2*pii*y))*sin(pii*x) - &
exp(at*t)*(-1 + 2*x)* &
(-4*betann*pii*cos(pii*y) + &
(4*betann*pii**2 - at*Rey)*y*sin(pii*y))) + &
sin(pii*x)*(2*exp(at*t)* &
(pii*(at*Rey*(-1 + x)*x + &
betann*(4 - 4*pii**2*(-1 + x)*x))* &
cos(pii*y) + (at*Rey + pii**2* &
(-(at*Rey*(-1 + x)*x) + &
2*betann*(-4 + pii**2*(-1 + x)*x)))*y* &
sin(pii*y)) + &
pii*Rey*cos(pii*x)* &
(2*(1 + 2*x*(-1 + x + pii**2*x**2))*y* &
sin(pii*y)**2 - &
pii*((-1 + x)**2*x**2 + &
(-1 - 2*(-1 + x)*x)*y**2)*sin(2*pii*y)) &
+ Rey*(-1 + 2*x)*sin(pii*x)* &
(y - pii**2*(-1 + x)*x*y + &
(-1 - pii**2*(-1 + x)*x)*y*cos(2*pii*y) + &
pii*(y**2 + (-1 + x)*x*(-1 + pii**2*y**2))* &
sin(2*pii*y))))/(exp(2*at*t)*Rey)
  end function tfwz_case10

  pure real(dp) function tftxx_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1 + betann)*(Rey* &
(-((alphaG*(c1**2 + c2**2) + &
c1*(c1 + c3)*epsylon)*Rey*Wi) + &
c1*exp(at*t)*(-1 + at*Wi)) + &
2*pii*cos(pii*x)* &
(pii*(-1 + x)*x* &
(exp(at*t) - c1*Rey*Wi*(-1 + xi))*y* &
cos(pii*y) + &
(exp(at*t)*(-1 + x)*x + &
Rey*Wi* &
(-(c1*(-1 + x)*x*(-1 + xi)) + &
c2*(-1 + 2*x)*xi*y))*sin(pii*y)) + &
2*sin(pii*x)* &
(pii*(-(c2*Rey*Wi*(-1 + x)*x*(-2 + xi)) + &
(-1 + 2*x)* &
(exp(at*t) - c1*Rey*Wi*(-1 + xi))*y)* &
cos(pii*y) + &
(exp(at*t)*(-1 + 2*x) + &
Rey*Wi* &
(-(c1*(-1 + 2*x)*(-1 + xi)) + &
c2*(-(pii**2*(-1 + x)*x) + xi)*y))* &
sin(pii*y))))/(exp(2*at*t)*Rey)
  end function tftxx_case10

  pure real(dp) function tftxy_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1 + betann)*(c2*Rey* &
(-((c1 + c3)*(alphaG + epsylon)*Rey*Wi) + &
exp(at*t)*(-1 + at*Wi)) + &
pii*(-1 + x)*x* &
(2*exp(at*t) + Rey*Wi*(2*c3 - (c1 + c3)*xi))* &
cos(pii*y)*sin(pii*x) + &
y*(-(pii*(-1 + 2*x)* &
(2*exp(at*t) + &
Rey*Wi*(2*c1 - (c1 + c3)*xi))*cos(pii*x) &
) + (-2*exp(at*t) + &
Rey*Wi* &
(c3*(-(pii**2*(-1 + x)*x) + xi) + &
c1*(-2 + pii**2*(-1 + x)*x + xi)))* &
sin(pii*x))*sin(pii*y)))/(exp(2*at*t)*Rey)
  end function tftxy_case10

  pure real(dp) function tftyy_case10(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    real(dp), parameter :: pii = 3.141592653589793_dp
    value = ((-1 + betann)*(Rey* &
(-((alphaG*(c2**2 + c3**2) + &
c3*(c1 + c3)*epsylon)*Rey*Wi) + &
c3*exp(at*t)*(-1 + at*Wi)) + &
2*pii*cos(pii*x)* &
(-(pii*(-1 + x)*x* &
(exp(at*t) - c3*Rey*Wi*(-1 + xi))*y* &
cos(pii*y)) + &
(-(exp(at*t)*(-1 + x)*x) + &
Rey*Wi* &
(c3*(-1 + x)*x*(-1 + xi) + &
c2*(-1 + 2*x)*(-2 + xi)*y))*sin(pii*y)) &
+ 2*sin(pii*x)* &
(pii*(-(c2*Rey*Wi*(-1 + x)*x*xi) - &
(-1 + 2*x)* &
(exp(at*t) - c3*Rey*Wi*(-1 + xi))*y)* &
cos(pii*y) + &
(exp(at*t)*(1 - 2*x) + &
Rey*Wi* &
(c3*(-1 + 2*x)*(-1 + xi) + &
c2*(-2 + pii**2*(-1 + x)*x + xi)*y))* &
sin(pii*y))))/(exp(2*at*t)*Rey)
  end function tftyy_case10

end program case10_fields_reference
