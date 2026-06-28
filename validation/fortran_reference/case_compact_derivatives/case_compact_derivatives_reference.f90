program case_compact_derivatives_reference
  !! Standalone Fortran reference for the fourth-order compact derivative
  !! operators (Etapa 13).
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the algorithm in
  !!   src/fortran_original/4Order_LPOG_Euler/derivs.f   (coef_4a, coeff_4a,
  !!                                                        rhsx_4a, rhsxx_4a,
  !!                                                        rhsy_4a, rhsyy_4a)
  !!   src/fortran_original/4Order_LPOG_Euler/solvers.f  (tridx, tridy)
  !! applied to a smooth analytic test function f(x,y) = sin(pi*x)*cos(pi*y),
  !! independent of any MMS case, since this validates the differentiation
  !! operator itself rather than a specific manufactured solution.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case_compact_derivatives_reference.f90 -o case_compact_derivatives_reference.exe
  !!   .\case_compact_derivatives_reference.exe output.csv 9 9
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer :: nx, ny, i, j
  real(dp) :: dx, dy
  character(len=512) :: outfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:), f(:,:)
  real(dp), allocatable :: exact_dfdx(:,:), exact_d2fdx2(:,:)
  real(dp), allocatable :: exact_dfdy(:,:), exact_d2fdy2(:,:)
  real(dp), allocatable :: num_dfdx(:,:), num_d2fdx2(:,:)
  real(dp), allocatable :: num_dfdy(:,:), num_d2fdy2(:,:)
  real(dp), allocatable :: a1x(:), b1x(:), c1x(:), a2x(:), b2x(:), c2x(:)
  real(dp), allocatable :: a1y(:), b1y(:), c1y(:), a2y(:), b2y(:), c2y(:)

  nx = 9
  ny = 9
  outfile = "case_compact_derivatives_9x9_default.csv"

  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) outfile = trim(arg)

  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg, *) nx

  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg, *) ny

  if (nx < 4 .or. ny < 4) error stop "nx and ny must be at least 4"

  dx = 1.0_dp / real(nx - 1, dp)
  dy = 1.0_dp / real(ny - 1, dp)

  allocate(xg(nx,ny), yg(nx,ny), f(nx,ny))
  allocate(exact_dfdx(nx,ny), exact_d2fdx2(nx,ny))
  allocate(exact_dfdy(nx,ny), exact_d2fdy2(nx,ny))
  allocate(num_dfdx(nx,ny), num_d2fdx2(nx,ny))
  allocate(num_dfdy(nx,ny), num_d2fdy2(nx,ny))

  do j = 1, ny
    do i = 1, nx
      xg(i,j) = real(i - 1, dp) * dx
      yg(i,j) = real(j - 1, dp) * dy
      f(i,j) = sin(pi * xg(i,j)) * cos(pi * yg(i,j))
      exact_dfdx(i,j) = pi * cos(pi * xg(i,j)) * cos(pi * yg(i,j))
      exact_d2fdx2(i,j) = -(pi**2) * sin(pi * xg(i,j)) * cos(pi * yg(i,j))
      exact_dfdy(i,j) = -pi * sin(pi * xg(i,j)) * sin(pi * yg(i,j))
      exact_d2fdy2(i,j) = -(pi**2) * sin(pi * xg(i,j)) * cos(pi * yg(i,j))
    end do
  end do

  allocate(a1x(nx), b1x(nx), c1x(nx), a2x(nx), b2x(nx), c2x(nx))
  allocate(a1y(ny), b1y(ny), c1y(ny), a2y(ny), b2y(ny), c2y(ny))
  call coef_4a(a1x, b1x, c1x, nx)
  call coeff_4a(a2x, b2x, c2x, nx)
  call coef_4a(a1y, b1y, c1y, ny)
  call coeff_4a(a2y, b2y, c2y, ny)

  call rhsx_4a(num_dfdx, f, nx, ny, dx)
  call tridx(a1x, b1x, c1x, num_dfdx, nx, ny)

  call rhsxx_4a(num_d2fdx2, f, nx, ny, dx)
  call tridx(a2x, b2x, c2x, num_d2fdx2, nx, ny)

  call rhsy_4a(num_dfdy, f, nx, ny, dy)
  call tridy(a1y, b1y, c1y, num_dfdy, nx, ny)

  call rhsyy_4a(num_d2fdy2, f, nx, ny, dy)
  call tridy(a2y, b2y, c2y, num_d2fdy2, nx, ny)

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,f,exact_dfdx,exact_d2fdx2,exact_dfdy,exact_d2fdy2,"// &
                   "num_dfdx,num_d2fdx2,num_dfdy,num_d2fdy2"
  do j = 1, ny
    do i = 1, nx
      write(20, '(*(ES26.17E3,:,","))') xg(i,j), yg(i,j), f(i,j), &
        exact_dfdx(i,j), exact_d2fdx2(i,j), exact_dfdy(i,j), exact_d2fdy2(i,j), &
        num_dfdx(i,j), num_d2fdx2(i,j), num_dfdy(i,j), num_d2fdy2(i,j)
    end do
  end do
  close(20)

contains

  !! Mirrors coef_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine coef_4a(a, b, c, lmax)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: a(lmax), b(lmax), c(lmax)
    integer :: l

    a(1) = 0.0_dp
    b(1) = 1.0_dp
    c(1) = 3.0_dp
    do l = 2, lmax - 1
      a(l) = 1.0_dp
      b(l) = 4.0_dp
      c(l) = 1.0_dp
    end do
    a(lmax) = 3.0_dp
    b(lmax) = 1.0_dp
    c(lmax) = 0.0_dp
  end subroutine coef_4a

  !! Mirrors coeff_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine coeff_4a(a, b, c, lmax)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: a(lmax), b(lmax), c(lmax)
    integer :: l

    a(1) = 0.0_dp
    b(1) = 1.0_dp
    c(1) = 11.0_dp
    do l = 2, lmax - 1
      a(l) = 1.0_dp
      b(l) = 10.0_dp
      c(l) = 1.0_dp
    end do
    a(lmax) = 11.0_dp
    b(lmax) = 1.0_dp
    c(lmax) = 0.0_dp
  end subroutine coeff_4a

  !! Mirrors rhsx_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine rhsx_4a(rhs, fc, n1, n2, h)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2), h
    real(dp), intent(out) :: rhs(n1,n2)
    integer :: i, j

    do j = 1, n2
      rhs(1,j) = (-17.0_dp*fc(1,j) + 9.0_dp*fc(2,j) + 9.0_dp*fc(3,j) - fc(4,j)) / (6.0_dp*h)
      do i = 2, n1 - 1
        rhs(i,j) = (3.0_dp*(fc(i+1,j) - fc(i-1,j))) / h
      end do
      rhs(n1,j) = (-17.0_dp*fc(n1,j) + 9.0_dp*fc(n1-1,j) + 9.0_dp*fc(n1-2,j) - fc(n1-3,j)) / (-6.0_dp*h)
    end do
  end subroutine rhsx_4a

  !! Mirrors rhsxx_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine rhsxx_4a(rhs, fc, n1, n2, h)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2), h
    real(dp), intent(out) :: rhs(n1,n2)
    real(dp) :: hh
    integer :: i, j

    hh = h * h
    do j = 1, n2
      rhs(1,j) = (39.0_dp*fc(1,j) - 81.0_dp*fc(2,j) + 45.0_dp*fc(3,j) - 3.0_dp*fc(4,j)) / (3.0_dp*hh)
      do i = 2, n1 - 1
        rhs(i,j) = (12.0_dp*(fc(i+1,j) + fc(i-1,j)) - 24.0_dp*fc(i,j)) / hh
      end do
      rhs(n1,j) = (39.0_dp*fc(n1,j) - 81.0_dp*fc(n1-1,j) + 45.0_dp*fc(n1-2,j) - 3.0_dp*fc(n1-3,j)) / (3.0_dp*hh)
    end do
  end subroutine rhsxx_4a

  !! Mirrors rhsy_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine rhsy_4a(rhs, fc, n1, n2, h)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2), h
    real(dp), intent(out) :: rhs(n1,n2)
    integer :: i, j

    do i = 1, n1
      rhs(i,1) = (-17.0_dp*fc(i,1) + 9.0_dp*fc(i,2) + 9.0_dp*fc(i,3) - fc(i,4)) / (6.0_dp*h)
      do j = 2, n2 - 1
        rhs(i,j) = (3.0_dp*(fc(i,j+1) - fc(i,j-1))) / h
      end do
      rhs(i,n2) = (-17.0_dp*fc(i,n2) + 9.0_dp*fc(i,n2-1) + 9.0_dp*fc(i,n2-2) - fc(i,n2-3)) / (-6.0_dp*h)
    end do
  end subroutine rhsy_4a

  !! Mirrors rhsyy_4a in src/fortran_original/4Order_LPOG_Euler/derivs.f.
  subroutine rhsyy_4a(rhs, fc, n1, n2, h)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2), h
    real(dp), intent(out) :: rhs(n1,n2)
    real(dp) :: hh
    integer :: i, j

    hh = h * h
    do i = 1, n1
      rhs(i,1) = (39.0_dp*fc(i,1) - 81.0_dp*fc(i,2) + 45.0_dp*fc(i,3) - 3.0_dp*fc(i,4)) / (3.0_dp*hh)
      do j = 2, n2 - 1
        rhs(i,j) = (12.0_dp*(fc(i,j+1) + fc(i,j-1)) - 24.0_dp*fc(i,j)) / hh
      end do
      rhs(i,n2) = (39.0_dp*fc(i,n2) - 81.0_dp*fc(i,n2-1) + 45.0_dp*fc(i,n2-2) - 3.0_dp*fc(i,n2-3)) / (3.0_dp*hh)
    end do
  end subroutine rhsyy_4a

  !! Mirrors tridx in src/fortran_original/4Order_LPOG_Euler/solvers.f.
  subroutine tridx(a, b, c, rhs, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: a(n1), b(n1), c(n1)
    real(dp), intent(inout) :: rhs(n1,n2)
    real(dp) :: gam(n1), bet, u(n1)
    integer :: i, j

    do j = 1, n2
      bet = b(1)
      u(1) = rhs(1,j) / bet
      do i = 2, n1
        gam(i) = c(i-1) / bet
        bet = b(i) - a(i)*gam(i)
        u(i) = (rhs(i,j) - a(i)*u(i-1)) / bet
      end do
      do i = n1 - 1, 1, -1
        u(i) = u(i) - gam(i+1)*u(i+1)
      end do
      do i = 1, n1
        rhs(i,j) = u(i)
      end do
    end do
  end subroutine tridx

  !! Mirrors tridy in src/fortran_original/4Order_LPOG_Euler/solvers.f.
  subroutine tridy(a, b, c, rhs, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: a(n2), b(n2), c(n2)
    real(dp), intent(inout) :: rhs(n1,n2)
    real(dp) :: gam(n2), bet, u(n2)
    integer :: i, j

    do i = 1, n1
      bet = b(1)
      u(1) = rhs(i,1) / bet
      do j = 2, n2
        gam(j) = c(j-1) / bet
        bet = b(j) - a(j)*gam(j)
        u(j) = (rhs(i,j) - a(j)*u(j-1)) / bet
      end do
      do j = n2 - 1, 1, -1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do
      do j = 1, n2
        rhs(i,j) = u(j)
      end do
    end do
  end subroutine tridy

end program case_compact_derivatives_reference
