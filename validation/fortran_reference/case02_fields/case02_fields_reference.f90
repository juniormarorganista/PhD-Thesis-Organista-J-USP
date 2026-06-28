program case02_fields_reference
  !! Standalone Fortran reference exporter for MMS Case 02 exact fields and
  !! source terms.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the original expressions from:
  !!
  !!   src/fortran_original/4Order_LPOG_Euler/case_02.f
  !!     ::exact_solutions_case_02(t)
  !!     ::source_term_case_02(t)
  !!
  !! Default parameters match the Case02Parameters/Case02SourceTermParameters
  !! Python defaults (Re=1, Wi=1, beta=0.5, epsilon=0, xi=0, alphaG=0, at=0,
  !! all polynomial a/b coefficients = 1), not the transient Rey/Wi/epsylon
  !! values currently sitting in par.nn.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case02_fields_reference.f90 -o case02_fields_reference.exe
  !!   .\case02_fields_reference.exe case02_fields_5x5_default.csv 5 5
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
  outfile = 'case02_fields_5x5_default.csv'

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
        u_case02(x, y, t), v_case02(x, y, t), wz_case02(x, y, t), psi_case02(x, y, t), &
        txx_case02(x, y, t), txy_case02(x, y, t), tyy_case02(x, y, t), &
        tfwz_case02(x, y, t), tftxx_case02(x, y, t), tftxy_case02(x, y, t), tftyy_case02(x, y, t)
    end do
  end do

  close(10)

contains

  pure real(dp) function u_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (16.0_dp*(-1.0_dp + x)**2*x**2*y*(-1.0_dp + 2.0_dp*y**2))/exp(at*t)
  end function u_case02

  pure real(dp) function v_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (-16.0_dp*(-1.0_dp + x)*x*(-1.0_dp + 2.0_dp*x)*y**2*(-1.0_dp + y**2))/exp(at*t)
  end function v_case02

  pure real(dp) function wz_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (16.0_dp*(-((-1.0_dp + x)**2*x**2) + (-1.0_dp + 6.0_dp*(x - 2.0_dp*x**3 + x**4))*y**2 &
      + (1.0_dp + 6.0_dp*(-1.0_dp + x)*x)*y**4))/exp(at*t)
  end function wz_case02

  pure real(dp) function psi_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp
    value = (8.0_dp*(-1.0_dp + x)**2*x**2*y**2*(-1.0_dp + y**2))/exp(at*t)
  end function psi_case02

  pure real(dp) function txx_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    value = -(((-1.0_dp + betann)*(a1xx + x*(a2xx + a3xx*x))*(b1xx + y*(b2xx + b3xx*y)))/exp(at*t))
  end function txx_case02

  pure real(dp) function txy_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    value = -(((-1.0_dp + betann)*(a1xy + x*(a2xy + a3xy*x))*(b1xy + y*(b2xy + b3xy*y)))/exp(at*t))
  end function txy_case02

  pure real(dp) function tyy_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: at = 0.0_dp, betann = 0.5_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    value = -(((-1.0_dp + betann)*(a1yy + x*(a2yy + a3yy*x))*(b1yy + y*(b2yy + b3yy*y)))/exp(at*t))
  end function tyy_case02

  pure real(dp) function tfwz_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, betann = 0.5_dp, at = 0.0_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp

    value = (betann*exp(at*t)* &
      (-64 + a2yy*b2yy*Rey - &
        2*a1xy*b3xy*Rey + 384*x - &
        2*a3xx*b2xx*Rey*x + &
        2*a3yy*b2yy*Rey*x - &
        2*a2xy*b3xy*Rey*x - &
        192*x**2 - 384*x**3 + &
        192*x**4 + &
        2*Rey* &
         (a2yy*b3yy - &
           2*a3xx*b3xx*x + &
           2*a3yy*b3yy*x)*y + &
        192*(1 + 12*(-1 + x)*x)* &
         y**2 + 192*y**4 - &
        a2xx*Rey*(b2xx + 2*b3xx*y) + &
        2*a3xy*Rey* &
         (b1xy - b3xy*x**2 + &
           b2xy*y + b3xy*y**2)) + &
     Rey*(-512*(-1 + x)*x*(-1 + 2*x)* &
         y*((-1 + x)**2*x**2 - &
           (1 + (-1 + x)*x)* &
            (1 + 2*(-1 + x)*x)*y**2 &
            + 3*(1 + (-1 + x)*x)* &
            (1 + 2*(-1 + x)*x)*y**4 &
            - 2*(1 + 3*(-1 + x)*x)* &
            y**6) + &
        exp(at*t)* &
         (-(a2yy*b2yy) + &
           2*a1xy*b3xy + &
           a2xx*(b2xx + 2*b3xx*y) - &
           2*a3xy* &
            (b1xy - b3xy*x**2 + &
              b2xy*y + b3xy*y**2) + &
           2* &
            (x* &
               (-(a3yy*b2yy) + &
                a2xy*b3xy - &
                8*at*(-1 + x)**2*x) &
               - b3yy* &
               (a2yy + 2*a3yy*x)*y + &
              8*at* &
               (-1 + &
                6*(x - 2*x**3 + x**4) &
                )*y**2 + &
              8*at* &
               (1 + 6*(-1 + x)*x)* &
               y**4 + &
              a3xx*x* &
               (b2xx + 2*b3xx*y)))))/ &
     (exp(2*at*t)*Rey)
  end function tfwz_case02

  pure real(dp) function tftxx_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp

    value = ((64*(-1 + betann)*exp(at*t)* &
        (-1 + x)*x*(-1 + 2*x)*y* &
        (-1 + 2*y**2))/Rey + &
     alphaG*(1 - betann)*Rey*Wi* &
      ((a1xx + x*(a2xx + a3xx*x))**2* &
         (b1xx + y*(b2xx + b3xx*y))** &
          2 + &
        (a1xy + x*(a2xy + a3xy*x))** &
          2* &
         (b1xy + y*(b2xy + b3xy*y))** &
          2) + &
     Wi*(16*(-1 + betann)*(-1 + x)*x* &
         (-1 + 2*x)* &
         (a1xx + x*(a2xx + a3xx*x))* &
         y**2*(b2xx + 2*b3xx*y)* &
         (-1 + y**2) + &
        at*(-1 + betann)*exp(at*t)* &
         (a1xx + x*(a2xx + a3xx*x))* &
         (b1xx + y*(b2xx + b3xx*y)) &
         - 16*(-1 + betann)* &
         (-1 + x)**2*x**2* &
         (a2xx + 2*a3xx*x)*y* &
         (-1 + 2*y**2)* &
         (b1xx + y*(b2xx + b3xx*y)) &
         + 64*(-1 + betann)*(-1 + x)* &
         x*(-1 + 2*x)* &
         (a1xx + x*(a2xx + a3xx*x))* &
         y*(-1 + 2*y**2)* &
         (b1xx + y*(b2xx + b3xx*y)) &
         + 32*(-1 + betann)* &
         (-1 + x)**2*x**2* &
         (a1xy + x*(a2xy + a3xy*x))* &
         (-1 + 6*y**2)* &
         (b1xy + y*(b2xy + b3xy*y)) &
         + 16*(1 - betann)*xi* &
         (4*(-1 + x)*x*(-1 + 2*x)* &
            (a1xx + &
              x*(a2xx + a3xx*x))*y* &
            (-1 + 2*y**2)* &
            (b1xx + &
              y*(b2xx + b3xx*y)) - &
           (1 + 6*(-1 + x)*x)* &
            (a1xy + &
              x*(a2xy + a3xy*x))* &
            y**2*(-1 + y**2)* &
            (b1xy + &
              y*(b2xy + b3xy*y)) + &
           (-1 + x)**2*x**2* &
            (a1xy + &
              x*(a2xy + a3xy*x))* &
            (-1 + 6*y**2)* &
            (b1xy + &
              y*(b2xy + b3xy*y)))) + &
     (1 - betann)* &
      (a1xx + x*(a2xx + a3xx*x))* &
      (b1xx + y*(b2xx + b3xx*y))* &
      (exp(at*t) + &
        epsylon*Rey*Wi* &
         (a1xx* &
            (b1xx + &
              y*(b2xx + b3xx*y)) + &
           a1yy* &
            (b1yy + &
              y*(b2yy + b3yy*y)) + &
           x* &
            (a2yy*b1yy + &
              a2yy*y* &
               (b2yy + b3yy*y) + &
              a2xx* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3xx*x* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3yy*x* &
               (b1yy + &
                y*(b2yy + b3yy*y))))) &
     )/exp(2*at*t)
  end function tftxx_case02

  pure real(dp) function tftxy_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp

    value = ((-16*(-1 + betann)*exp(at*t)* &
        ((-1 + x)**2*x**2 + &
          (-1 - &
             6*(-1 + x)*x* &
              (1 + (-1 + x)*x))*y**2 &
           + (1 + 6*(-1 + x)*x)*y**4) &
        )/Rey + &
     alphaG*(1 - betann)*Rey*Wi* &
      (a1xy + x*(a2xy + a3xy*x))* &
      (b1xy + y*(b2xy + b3xy*y))* &
      ((a1xx + x*(a2xx + a3xx*x))* &
         (b1xx + y*(b2xx + b3xx*y)) &
         + (a1yy + &
           x*(a2yy + a3yy*x))* &
         (b1yy + y*(b2yy + b3yy*y))) &
      - (-1 + betann)* &
      (a1xy + x*(a2xy + a3xy*x))* &
      (b1xy + y*(b2xy + b3xy*y))* &
      (exp(at*t) + &
        epsylon*Rey*Wi* &
         (a1xx* &
            (b1xx + &
              y*(b2xx + b3xx*y)) + &
           a1yy* &
            (b1yy + &
              y*(b2yy + b3yy*y)) + &
           x* &
            (a2yy*b1yy + &
              a2yy*y* &
               (b2yy + b3yy*y) + &
              a2xx* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3xx*x* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3yy*x* &
               (b1yy + &
                y*(b2yy + b3yy*y))))) &
       + (-1 + betann)*Wi* &
      (16*(-1 + x)*x*(-1 + 2*x)* &
         (a1xy + x*(a2xy + a3xy*x))* &
         y**2*(b2xy + 2*b3xy*y)* &
         (-1 + y**2) - &
        16*(1 + 6*(-1 + x)*x)* &
         (a1xx + x*(a2xx + a3xx*x))* &
         y**2*(-1 + y**2)* &
         (b1xx + y*(b2xx + b3xx*y)) &
         + at*exp(at*t)* &
         (a1xy + x*(a2xy + a3xy*x))* &
         (b1xy + y*(b2xy + b3xy*y)) &
         - 16*(-1 + x)**2*x**2* &
         (a2xy + 2*a3xy*x)*y* &
         (-1 + 2*y**2)* &
         (b1xy + y*(b2xy + b3xy*y)) &
         + 16*(-1 + x)**2*x**2* &
         (a1yy + x*(a2yy + a3yy*x))* &
         (-1 + 6*y**2)* &
         (b1yy + y*(b2yy + b3yy*y)) &
         + 8*xi* &
         ((-1 + x)**2*x**2 + &
           (-1 - &
              6*(-1 + x)*x* &
               (1 + (-1 + x)*x))*y**2 &
             + (1 + 6*(-1 + x)*x)* &
            y**4)* &
         (a1xx* &
            (b1xx + &
              y*(b2xx + b3xx*y)) + &
           a1yy* &
            (b1yy + &
              y*(b2yy + b3yy*y)) + &
           x* &
            (a2yy*b1yy + &
              a2yy*y* &
               (b2yy + b3yy*y) + &
              a2xx* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3xx*x* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3yy*x* &
               (b1yy + &
                y*(b2yy + b3yy*y))))) &
     )/exp(2*at*t)
  end function tftxy_case02

  pure real(dp) function tftyy_case02(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp), parameter :: Rey = 1.0_dp, Wi = 1.0_dp, betann = 0.5_dp
    real(dp), parameter :: epsylon = 0.0_dp, xi = 0.0_dp, alphaG = 0.0_dp, at = 0.0_dp
    real(dp), parameter :: a1xy = 1.0_dp, a2xy = 1.0_dp, a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp, b2xy = 1.0_dp, b3xy = 1.0_dp
    real(dp), parameter :: a1yy = 1.0_dp, a2yy = 1.0_dp, a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp, b2yy = 1.0_dp, b3yy = 1.0_dp
    real(dp), parameter :: a1xx = 1.0_dp, a2xx = 1.0_dp, a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp, b2xx = 1.0_dp, b3xx = 1.0_dp

    value = ((-64*(-1 + betann)*exp(at*t)* &
        (-1 + x)*x*(-1 + 2*x)*y* &
        (-1 + 2*y**2))/Rey + &
     alphaG*(1 - betann)*Rey*Wi* &
      ((a1xy + x*(a2xy + a3xy*x))**2* &
         (b1xy + y*(b2xy + b3xy*y))** &
          2 + &
        (a1yy + x*(a2yy + a3yy*x))** &
          2* &
         (b1yy + y*(b2yy + b3yy*y))** &
          2) + &
     Wi*(16*(-1 + betann)*(-1 + x)*x* &
         (-1 + 2*x)* &
         (a1yy + x*(a2yy + a3yy*x))* &
         y**2*(b2yy + 2*b3yy*y)* &
         (-1 + y**2) - &
        32*(-1 + betann)* &
         (1 + 6*(-1 + x)*x)* &
         (a1xy + x*(a2xy + a3xy*x))* &
         y**2*(-1 + y**2)* &
         (b1xy + y*(b2xy + b3xy*y)) &
         + at*(-1 + betann)* &
         exp(at*t)* &
         (a1yy + x*(a2yy + a3yy*x))* &
         (b1yy + y*(b2yy + b3yy*y)) &
         - 16*(-1 + betann)* &
         (-1 + x)**2*x**2* &
         (a2yy + 2*a3yy*x)*y* &
         (-1 + 2*y**2)* &
         (b1yy + y*(b2yy + b3yy*y)) &
         - 64*(-1 + betann)*(-1 + x)* &
         x*(-1 + 2*x)* &
         (a1yy + x*(a2yy + a3yy*x))* &
         y*(-1 + 2*y**2)* &
         (b1yy + y*(b2yy + b3yy*y)) &
         + 16*(1 - betann)*xi* &
         (-((1 + 6*(-1 + x)*x)* &
              (a1xy + &
                x*(a2xy + a3xy*x))* &
              y**2*(-1 + y**2)* &
              (b1xy + &
                y*(b2xy + b3xy*y))) &
            + (-1 + x)**2*x**2* &
            (a1xy + &
              x*(a2xy + a3xy*x))* &
            (-1 + 6*y**2)* &
            (b1xy + &
              y*(b2xy + b3xy*y)) - &
           4*(-1 + x)*x*(-1 + 2*x)* &
            (a1yy + &
              x*(a2yy + a3yy*x))*y* &
            (-1 + 2*y**2)* &
            (b1yy + &
              y*(b2yy + b3yy*y)))) + &
     (1 - betann)* &
      (a1yy + x*(a2yy + a3yy*x))* &
      (b1yy + y*(b2yy + b3yy*y))* &
      (exp(at*t) + &
        epsylon*Rey*Wi* &
         (a1xx* &
            (b1xx + &
              y*(b2xx + b3xx*y)) + &
           a1yy* &
            (b1yy + &
              y*(b2yy + b3yy*y)) + &
           x* &
            (a2yy*b1yy + &
              a2yy*y* &
               (b2yy + b3yy*y) + &
              a2xx* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3xx*x* &
               (b1xx + &
                y*(b2xx + b3xx*y)) + &
              a3yy*x* &
               (b1yy + &
                y*(b2yy + b3yy*y))))) &
     )/exp(2*at*t)
  end function tftyy_case02

end program case02_fields_reference
