program case_lele_filter_reference
  !! Standalone Fortran reference for the Lele-type compact spatial filter
  !! (Etapa 13).
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the algorithm in
  !!   src/fortran_original/4Order_LPOG_Euler/loop.f
  !!     (rhs_tridfx_wz, lhs_tridfx, and their y-direction counterparts)
  !!   src/fortran_original/4Order_LPOG_Euler/solvers.f  (tridx, tridy)
  !! applied to a smooth analytic test field f(x,y) = sin(3*pi*x)*cos(2*pi*y)
  !! plus a high-frequency component sin(20*pi*x), so that the filter has a
  !! non-trivial effect to check, with alphaf = 0.48 matching
  !! src/fortran_original/4Order_LPOG_Euler/par.nn.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case_lele_filter_reference.f90 -o case_lele_filter_reference.exe
  !!   .\case_lele_filter_reference.exe output.csv 13 13
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: alphaf = 0.48_dp
  real(dp), parameter :: af = (11.0_dp + 10.0_dp * alphaf) / 16.0_dp
  real(dp), parameter :: bf = (15.0_dp + 34.0_dp * alphaf) / 64.0_dp
  real(dp), parameter :: cf = (-3.0_dp + 6.0_dp * alphaf) / 32.0_dp
  real(dp), parameter :: df = (1.0_dp - 2.0_dp * alphaf) / 64.0_dp
  integer :: nx, ny, i, j
  real(dp) :: dx, dy
  character(len=512) :: outfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:), f(:,:)
  real(dp), allocatable :: filtered_x(:,:), filtered_y(:,:)
  real(dp), allocatable :: afilx(:), bfilx(:), cfilx(:)
  real(dp), allocatable :: afily(:), bfily(:), cfily(:)

  nx = 13
  ny = 13
  outfile = "case_lele_filter_13x13_default.csv"

  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) outfile = trim(arg)

  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg, *) nx

  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg, *) ny

  if (nx < 7 .or. ny < 7) error stop "nx and ny must be at least 7"

  dx = 1.0_dp / real(nx - 1, dp)
  dy = 1.0_dp / real(ny - 1, dp)

  allocate(xg(nx,ny), yg(nx,ny), f(nx,ny))
  allocate(filtered_x(nx,ny), filtered_y(nx,ny))

  do j = 1, ny
    do i = 1, nx
      xg(i,j) = real(i - 1, dp) * dx
      yg(i,j) = real(j - 1, dp) * dy
      f(i,j) = sin(3.0_dp * pi * xg(i,j)) * cos(2.0_dp * pi * yg(i,j)) &
             + 0.1_dp * sin(20.0_dp * pi * xg(i,j))
    end do
  end do

  allocate(afilx(nx), bfilx(nx), cfilx(nx))
  allocate(afily(ny), bfily(ny), cfily(ny))
  call lhs_tridfx(afilx, bfilx, cfilx, nx)
  call lhs_tridfy(afily, bfily, cfily, ny)

  call rhs_tridfx(filtered_x, f, nx, ny)
  call tridx(afilx, bfilx, cfilx, filtered_x, nx, ny)

  call rhs_tridfy(filtered_y, f, nx, ny)
  call tridy(afily, bfily, cfily, filtered_y, nx, ny)

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,f,filtered_x,filtered_y"
  do j = 1, ny
    do i = 1, nx
      write(20, '(*(ES26.17E3,:,","))') xg(i,j), yg(i,j), f(i,j), filtered_x(i,j), filtered_y(i,j)
    end do
  end do
  close(20)

contains

  !! Mirrors lhs_tridfx/lhs_tridfy in src/fortran_original/4Order_LPOG_Euler/loop.f.
  subroutine lhs_tridfx(a, b, c, lmax)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: a(lmax), b(lmax), c(lmax)

    a = 0.0_dp
    b = 1.0_dp
    c = 0.0_dp
    a(4:lmax-3) = alphaf
    c(4:lmax-3) = alphaf
  end subroutine lhs_tridfx

  subroutine lhs_tridfy(a, b, c, lmax)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: a(lmax), b(lmax), c(lmax)

    call lhs_tridfx(a, b, c, lmax)
  end subroutine lhs_tridfy

  !! Mirrors rhs_tridfx_wz in src/fortran_original/4Order_LPOG_Euler/loop.f.
  subroutine rhs_tridfx(rhs, fc, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2)
    real(dp), intent(out) :: rhs(n1,n2)
    integer :: i, j

    do j = 1, n2
      rhs(1,j) = (15.0_dp*fc(1,j) + 4.0_dp*fc(2,j) - 6.0_dp*fc(3,j) + 4.0_dp*fc(4,j) - fc(5,j)) / 16.0_dp
      rhs(2,j) = (fc(1,j) + 12.0_dp*fc(2,j) + 6.0_dp*fc(3,j) - 4.0_dp*fc(4,j) + fc(5,j)) / 16.0_dp
      rhs(3,j) = (-fc(1,j) + 4.0_dp*fc(2,j) + 10.0_dp*fc(3,j) + 4.0_dp*fc(4,j) - fc(5,j)) / 16.0_dp
      do i = 4, n1 - 3
        rhs(i,j) = af*fc(i,j) + bf*(fc(i+1,j)+fc(i-1,j)) + cf*(fc(i+2,j)+fc(i-2,j)) + df*(fc(i+3,j)+fc(i-3,j))
      end do
      rhs(n1-2,j) = (-fc(n1,j) + 4.0_dp*fc(n1-1,j) + 10.0_dp*fc(n1-2,j) + 4.0_dp*fc(n1-3,j) - fc(n1-4,j)) / 16.0_dp
      rhs(n1-1,j) = (fc(n1,j) + 12.0_dp*fc(n1-1,j) + 6.0_dp*fc(n1-2,j) - 4.0_dp*fc(n1-3,j) + fc(n1-4,j)) / 16.0_dp
      rhs(n1,j) = (15.0_dp*fc(n1,j) + 4.0_dp*fc(n1-1,j) - 6.0_dp*fc(n1-2,j) + 4.0_dp*fc(n1-3,j) - fc(n1-4,j)) / 16.0_dp
    end do
  end subroutine rhs_tridfx

  !! Mirrors rhs_tridfy_wz in src/fortran_original/4Order_LPOG_Euler/loop.f.
  subroutine rhs_tridfy(rhs, fc, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: fc(n1,n2)
    real(dp), intent(out) :: rhs(n1,n2)
    integer :: i, j

    do i = 1, n1
      rhs(i,1) = (15.0_dp*fc(i,1) + 4.0_dp*fc(i,2) - 6.0_dp*fc(i,3) + 4.0_dp*fc(i,4) - fc(i,5)) / 16.0_dp
      rhs(i,2) = (fc(i,1) + 12.0_dp*fc(i,2) + 6.0_dp*fc(i,3) - 4.0_dp*fc(i,4) + fc(i,5)) / 16.0_dp
      rhs(i,3) = (-fc(i,1) + 4.0_dp*fc(i,2) + 10.0_dp*fc(i,3) + 4.0_dp*fc(i,4) - fc(i,5)) / 16.0_dp
      do j = 4, n2 - 3
        rhs(i,j) = af*fc(i,j) + bf*(fc(i,j+1)+fc(i,j-1)) + cf*(fc(i,j+2)+fc(i,j-2)) + df*(fc(i,j+3)+fc(i,j-3))
      end do
      rhs(i,n2-2) = (-fc(i,n2) + 4.0_dp*fc(i,n2-1) + 10.0_dp*fc(i,n2-2) + 4.0_dp*fc(i,n2-3) - fc(i,n2-4)) / 16.0_dp
      rhs(i,n2-1) = (fc(i,n2) + 12.0_dp*fc(i,n2-1) + 6.0_dp*fc(i,n2-2) - 4.0_dp*fc(i,n2-3) + fc(i,n2-4)) / 16.0_dp
      rhs(i,n2) = (15.0_dp*fc(i,n2) + 4.0_dp*fc(i,n2-1) - 6.0_dp*fc(i,n2-2) + 4.0_dp*fc(i,n2-3) - fc(i,n2-4)) / 16.0_dp
    end do
  end subroutine rhs_tridfy

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

end program case_lele_filter_reference
