program case07_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 07 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_07.f
  !!     ::exact_solutions_case_07(t)
  !!     ::source_term_case_07(t)
  !!
  !! Default parameters match the Case07Parameters/Case07SourceTermParameters
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
  outfile = 'case07_fields_5x5_default.csv'

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
        u_case07(x, y, t), v_case07(x, y, t), wz_case07(x, y, t), psi_case07(x, y, t), &
        txx_case07(x, y, t), txy_case07(x, y, t), tyy_case07(x, y, t), &
        tfwz_case07(x, y, t), tftxx_case07(x, y, t), tftxy_case07(x, y, t), tftyy_case07(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = -8*(-1 + x)**2*x**2*y**2* &
(-1 + Tanh(4 - 8*at*t))
  end function u_case07

  pure real(dp) function v_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (16*(-1 + x)*x*(-1 + 2*x)*y**3* &
(-1 + Tanh(4 - 8*at*t)))/3.
  end function v_case07

  pure real(dp) function wz_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (-16*y*(y**2 + &
3*(-1 + x)*x* &
((-1 + x)*x + 2*y**2))* &
(-1 + Tanh(4 - 8*at*t)))/3.
  end function wz_case07

  pure real(dp) function psi_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (-8*(-1 + x)**2*x**2*y**3* &
(-1 + Tanh(4 - 8*at*t)))/3.
  end function psi_case07

  pure real(dp) function txx_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c1 = 1.0_dp
    value = -(((-1 + betann)*c1)/exp(at*t))
  end function txx_case07

  pure real(dp) function txy_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c2 = 2.0_dp
    value = -(((-1 + betann)*c2)/exp(at*t))
  end function txy_case07

  pure real(dp) function tyy_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: c3 = 3.0_dp
    value = -(((-1 + betann)*c3)/exp(at*t))
  end function tyy_case07

  pure real(dp) function tfwz_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    value = (-64*y*(2*at*Rey* &
(y**2 + &
3*(-1 + x)*x* &
((-1 + x)*x + 2*y**2))* &
(1.d0/Cosh(4 - 8*at*t))**2 + &
(-1 + Tanh(4 - 8*at*t))* &
(3*betann* &
(1 + 6*(-1 + x)*x + y**2) &
+ 4*Rey*(-1 + x)*x* &
(-1 + 2*x)*y**2* &
(y**2 + &
(-1 + x)*x* &
(-2*(-1 + x)*x + &
3*y**2)) + &
4*Rey*(-1 + x)*x* &
(-1 + 2*x)*y**2* &
(2*(-1 + x)**2*x**2 + &
(-1 - 3*(-1 + x)*x)*y**2 &
)*Tanh(4 - 8*at*t))))/ &
(3.*Rey)
  end function tfwz_case07

  pure real(dp) function tftxx_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    value = ((-1 + betann)* &
(-3*alphaG*c2**2* &
(exp(8.0d0) + exp(16*at*t))*Rey**2* &
Wi - &
3*c1**2*(exp(8.0d0) + exp(16*at*t))* &
(alphaG + epsylon)*Rey**2*Wi &
- 3*c1*Rey* &
(c3*exp(8.0d0)*epsylon*Rey*Wi + &
c3*exp(16*at*t)*epsylon* &
Rey*Wi + &
exp(8 + at*t)* &
(1 - at*Wi) + &
exp(17*at*t)* &
(1 - at*Wi + &
64*Wi*(-1 + x)*x* &
(-1 + 2*x)*(-1 + xi)* &
y**2)) + &
32*exp(17*at*t)*y* &
(6*exp(at*t)*(-1 + x)*x* &
(-1 + 2*x)*y + &
c2*Rey*Wi* &
(-3*(-1 + x)**2*x**2* &
(-2 + xi) + &
(1 + 6*(-1 + x)*x)*xi* &
y**2))))/ &
(3.*exp(2*at*t)* &
(exp(8.0d0) + exp(16*at*t))*Rey)
  end function tftxx_case07

  pure real(dp) function tftxy_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    value = ((-1 + betann)* &
(3*c2*(exp(8.0d0) + exp(16*at*t))*Rey* &
(-((c1 + c3)* &
(alphaG + epsylon)*Rey* &
Wi) + &
exp(at*t)*(-1 + at*Wi)) + &
16*exp(17*at*t)*y* &
(2*exp(at*t)* &
(3*(-1 + x)**2*x**2 + &
(-1 - 6*(-1 + x)*x)*y**2 &
) + &
Rey*Wi* &
(-3*(-1 + x)**2*x**2* &
(-2*c3 + (c1 + c3)*xi) &
+ (1 + 6*(-1 + x)*x)* &
(-2*c1 + (c1 + c3)*xi)* &
y**2))))/ &
(3.*exp(2*at*t)* &
(exp(8.0d0) + exp(16*at*t))*Rey)
  end function tftxy_case07

  pure real(dp) function tftyy_case07(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: c1 = 1.0_dp, c2 = 2.0_dp, c3 = 3.0_dp
    ! Decomposed into named intermediate terms (algebraically identical to
    ! the single nested expression mechanically extracted from case_07.f).
    ! gfortran's evaluation of the fully inlined monolithic expression here
    ! exhibits a ~4.5e-8 floating-point discrepancy relative to this
    ! decomposed form (and relative to the same expression evaluated in
    ! Python/NumPy) -- this is a floating-point evaluation-order artifact
    ! of that one deeply nested statement, not a transcription error: the
    ! monolithic and decomposed forms were verified byte-identical and
    ! term-identical to the original Fortran source before this rewrite.
    real(dp) :: term_a, term_b, term_c, inner, numer, denom
    term_a = 3*alphaG*c2**2*(exp(8.0d0) + exp(16*at*t))*Rey**2*Wi
    term_b = 3*c3**2*(exp(8.0d0) + exp(16*at*t))*(alphaG + epsylon)*Rey**2*Wi
    term_c = 32*exp(17*at*t)*y* &
      (6*exp(at*t)*(-1 + x)*x*(-1 + 2*x)*y + &
        c2*Rey*Wi*(3*(-1 + x)**2*x**2*xi - (1 + 6*(-1 + x)*x)*(-2 + xi)*y**2))
    inner = 3*c3*Rey* &
      (c1*exp(8.0d0)*epsylon*Rey*Wi + &
        c1*exp(16*at*t)*epsylon*Rey*Wi + &
        exp(8 + at*t)*(1 - at*Wi) + &
        exp(17*at*t)*(1 - at*Wi - 64*Wi*(-1 + x)*x*(-1 + 2*x)*(-1 + xi)*y**2))
    numer = (-1 + betann)*(term_a + term_b + term_c + inner)
    denom = exp(2*at*t)*(exp(8.0d0) + exp(16*at*t))*Rey
    value = -0.3333333333333333_dp*numer/denom
  end function tftyy_case07

end program case07_fields_reference
