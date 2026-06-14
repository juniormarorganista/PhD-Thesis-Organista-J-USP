program case12_euler_step_reference
  !! Standalone Fortran reference for one explicit Euler step of the Case 12
  !! ``ty_sim == 1`` formulation (wz, Txx, Txy, Tyy).
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It is a verbatim port of the right-hand sides computed by
  !! ``subroutine drv`` (src/fortran_original/4Order_LPOG_Euler/drvs.f) and the
  !! explicit Euler update performed by ``subroutine solver_euler``
  !! (src/fortran_original/4Order_LPOG_Euler/euler_rk.f), using the same
  !! second-order finite-difference operators as
  !! mms_viscoelastic.numerics.derivatives (NOT the 4th-order LPOG compact
  !! operators nor the Lele compact filter, which are out of scope).
  !!
  !! The input CSV (generate_input_csv.py) provides the t=0 exact fields
  !! (u, v, wz, txx, txy, tyy) and source terms (tfwz, tftxx, tftxy, tftyy) on
  !! a uniform grid covering [0,1] x [0,1]. Only interior grid points are
  !! advanced by one Euler step of size dt; boundary values are copied
  !! unchanged.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case12_euler_step_reference.f90 -o case12_euler_step_reference.exe
  !!   .\case12_euler_step_reference.exe input.csv output.csv 9 9 1.0e-5
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  integer :: nx, ny, i, j
  real(dp) :: dx, dy, dt
  real(dp) :: Re, Wi, beta, epsilon_, xi, alphaG, one_minus_beta
  character(len=512) :: infile, outfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:)
  real(dp), allocatable :: ux(:,:), uy(:,:), wz(:,:), txx(:,:), txy(:,:), tyy(:,:)
  real(dp), allocatable :: tfwz(:,:), tftxx(:,:), tftxy(:,:), tftyy(:,:)
  real(dp), allocatable :: dudx(:,:), dudy(:,:), dvdx(:,:), dvdy(:,:)
  real(dp), allocatable :: duwzdx(:,:), dvwzdy(:,:), diffusion(:,:)
  real(dp), allocatable :: dtxxdx(:,:), d2txx_dxdy(:,:)
  real(dp), allocatable :: d2txy_dx2(:,:), d2txy_dy2(:,:)
  real(dp), allocatable :: dtyydx(:,:), d2tyy_dxdy(:,:)
  real(dp), allocatable :: duTxxdx(:,:), dvTxxdy(:,:)
  real(dp), allocatable :: duTxydx(:,:), dvTxydy(:,:)
  real(dp), allocatable :: duTyydx(:,:), dvTyydy(:,:)
  real(dp), allocatable :: f_tij(:,:)
  real(dp), allocatable :: dvz(:,:), dvTxx(:,:), dvTxy(:,:), dvTyy(:,:)
  real(dp), allocatable :: wz_new(:,:), txx_new(:,:), txy_new(:,:), tyy_new(:,:)

  ! Defaults match EulerStepParameters() in
  ! mms_viscoelastic.numerics.time_stepping (and the Case12*Parameters
  ! defaults validated against the Fortran field/source-term references).
  nx = 9
  ny = 9
  Re = 1.0_dp
  Wi = 1.0_dp
  beta = 0.5_dp
  epsilon_ = 0.0_dp
  xi = 0.0_dp
  alphaG = 0.0_dp
  dt = 1.0e-5_dp
  infile = "case12_euler_step_input_9x9.csv"
  outfile = "case12_euler_step_output_9x9_dt1e-05.csv"

  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) infile = trim(arg)

  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) outfile = trim(arg)

  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg, *) nx

  call get_command_argument(4, arg)
  if (len_trim(arg) > 0) read(arg, *) ny

  call get_command_argument(5, arg)
  if (len_trim(arg) > 0) read(arg, *) dt

  if (nx <= 2 .or. ny <= 2) error stop "nx and ny must be greater than 2"

  dx = 1.0_dp / real(nx - 1, dp)
  dy = 1.0_dp / real(ny - 1, dp)
  one_minus_beta = 1.0_dp - beta

  allocate(xg(nx,ny), yg(nx,ny))
  allocate(ux(nx,ny), uy(nx,ny), wz(nx,ny), txx(nx,ny), txy(nx,ny), tyy(nx,ny))
  allocate(tfwz(nx,ny), tftxx(nx,ny), tftxy(nx,ny), tftyy(nx,ny))

  open(unit=10, file=trim(infile), status="old", action="read")
  read(10, *)  ! skip header
  do j = 1, ny
    do i = 1, nx
      read(10, *) xg(i,j), yg(i,j), ux(i,j), uy(i,j), wz(i,j), txx(i,j), txy(i,j), tyy(i,j), &
                  tfwz(i,j), tftxx(i,j), tftxy(i,j), tftyy(i,j)
    end do
  end do
  close(10)

  allocate(dudx(nx,ny), dudy(nx,ny), dvdx(nx,ny), dvdy(nx,ny))
  call derivative_x(ux, dx, nx, ny, dudx)
  call derivative_y(ux, dy, nx, ny, dudy)
  call derivative_x(uy, dx, nx, ny, dvdx)
  call derivative_y(uy, dy, nx, ny, dvdy)

  ! --- dvz: right-hand side of the wz transport equation ---
  allocate(duwzdx(nx,ny), dvwzdy(nx,ny), diffusion(nx,ny))
  call derivative_x(ux*wz, dx, nx, ny, duwzdx)
  call derivative_y(uy*wz, dy, nx, ny, dvwzdy)
  call laplacian(wz, dx, dy, nx, ny, diffusion)

  allocate(dtxxdx(nx,ny), d2txx_dxdy(nx,ny))
  call derivative_x(txx, dx, nx, ny, dtxxdx)
  call derivative_y(dtxxdx, dy, nx, ny, d2txx_dxdy)

  allocate(d2txy_dx2(nx,ny), d2txy_dy2(nx,ny))
  call second_derivative_x(txy, dx, nx, ny, d2txy_dx2)
  call second_derivative_y(txy, dy, nx, ny, d2txy_dy2)

  allocate(dtyydx(nx,ny), d2tyy_dxdy(nx,ny))
  call derivative_x(tyy, dx, nx, ny, dtyydx)
  call derivative_y(dtyydx, dy, nx, ny, d2tyy_dxdy)

  allocate(dvz(nx,ny))
  dvz = -duwzdx - dvwzdy + (beta/Re)*diffusion &
      + d2txx_dxdy + d2txy_dy2 - d2txy_dx2 - d2tyy_dxdy - tfwz

  ! --- dvTxx, dvTxy, dvTyy: right-hand sides of the stress equations ---
  allocate(duTxxdx(nx,ny), dvTxxdy(nx,ny))
  call derivative_x(ux*txx, dx, nx, ny, duTxxdx)
  call derivative_y(uy*txx, dy, nx, ny, dvTxxdy)

  allocate(duTxydx(nx,ny), dvTxydy(nx,ny))
  call derivative_x(ux*txy, dx, nx, ny, duTxydx)
  call derivative_y(uy*txy, dy, nx, ny, dvTxydy)

  allocate(duTyydx(nx,ny), dvTyydy(nx,ny))
  call derivative_x(ux*tyy, dx, nx, ny, duTyydx)
  call derivative_y(uy*tyy, dy, nx, ny, dvTyydy)

  allocate(f_tij(nx,ny))
  f_tij = 1.0_dp + (epsilon_*Re*Wi/one_minus_beta)*(txx + tyy)

  allocate(dvTxx(nx,ny), dvTxy(nx,ny), dvTyy(nx,ny))

  dvTxx = -f_tij*txx/Wi - duTxxdx - dvTxxdy &
        + 2.0_dp*txx*dudx + 2.0_dp*txy*dudy &
        - 2.0_dp*xi*txx*dudx - xi*txy*(dudy + dvdx) &
        - (alphaG*Re/one_minus_beta)*(txx*txx + txy*txy) &
        + (2.0_dp*one_minus_beta/(Re*Wi))*dudx + tftxx/Wi

  dvTxy = -f_tij*txy/Wi - duTxydx - dvTxydy &
        + txx*dvdx + tyy*dudy &
        - xi*txy*dudx - 0.5_dp*xi*tyy*(dudy + dvdx) &
        - 0.5_dp*xi*txx*(dudy + dvdx) - xi*txy*dvdy &
        - (alphaG*Re/one_minus_beta)*(txx*txy + txy*tyy) &
        + (one_minus_beta/(Re*Wi))*(dudy + dvdx) + tftxy/Wi

  dvTyy = -f_tij*tyy/Wi - duTyydx - dvTyydy &
        + 2.0_dp*txy*dvdx + 2.0_dp*tyy*dvdy &
        - xi*txy*(dvdx + dudy) - 2.0_dp*xi*tyy*dvdy &
        - (alphaG*Re/one_minus_beta)*(txy*txy + tyy*tyy) &
        + (2.0_dp*one_minus_beta/(Re*Wi))*dvdy + tftyy/Wi

  ! --- explicit Euler update, interior points only ---
  allocate(wz_new(nx,ny), txx_new(nx,ny), txy_new(nx,ny), tyy_new(nx,ny))
  wz_new = wz
  txx_new = txx
  txy_new = txy
  tyy_new = tyy

  do j = 2, ny - 1
    do i = 2, nx - 1
      wz_new(i,j) = wz(i,j) + dvz(i,j) * dt
      txx_new(i,j) = txx(i,j) + dvTxx(i,j) * dt
      txy_new(i,j) = txy(i,j) + dvTxy(i,j) * dt
      tyy_new(i,j) = tyy(i,j) + dvTyy(i,j) * dt
    end do
  end do

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,wz,txx,txy,tyy"
  do j = 1, ny
    do i = 1, nx
      write(20, '(ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3)') &
        xg(i,j), ',', yg(i,j), ',', wz_new(i,j), ',', txx_new(i,j), ',', txy_new(i,j), ',', tyy_new(i,j)
    end do
  end do
  close(20)

contains

  subroutine derivative_x(f, dxval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval
    real(dp), intent(out) :: out(n1,n2)
    integer :: ii, jj
    do jj = 1, n2
      do ii = 2, n1 - 1
        out(ii,jj) = (f(ii+1,jj) - f(ii-1,jj)) / (2.0_dp*dxval)
      end do
      out(1,jj) = (-3.0_dp*f(1,jj) + 4.0_dp*f(2,jj) - f(3,jj)) / (2.0_dp*dxval)
      out(n1,jj) = (3.0_dp*f(n1,jj) - 4.0_dp*f(n1-1,jj) + f(n1-2,jj)) / (2.0_dp*dxval)
    end do
  end subroutine derivative_x

  subroutine derivative_y(f, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dyval
    real(dp), intent(out) :: out(n1,n2)
    integer :: ii, jj
    do ii = 1, n1
      do jj = 2, n2 - 1
        out(ii,jj) = (f(ii,jj+1) - f(ii,jj-1)) / (2.0_dp*dyval)
      end do
      out(ii,1) = (-3.0_dp*f(ii,1) + 4.0_dp*f(ii,2) - f(ii,3)) / (2.0_dp*dyval)
      out(ii,n2) = (3.0_dp*f(ii,n2) - 4.0_dp*f(ii,n2-1) + f(ii,n2-2)) / (2.0_dp*dyval)
    end do
  end subroutine derivative_y

  subroutine second_derivative_x(f, dxval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval
    real(dp), intent(out) :: out(n1,n2)
    integer :: ii, jj
    out = 0.0_dp
    do jj = 1, n2
      do ii = 2, n1 - 1
        out(ii,jj) = (f(ii+1,jj) - 2.0_dp*f(ii,jj) + f(ii-1,jj)) / dxval**2
      end do
    end do
  end subroutine second_derivative_x

  subroutine second_derivative_y(f, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dyval
    real(dp), intent(out) :: out(n1,n2)
    integer :: ii, jj
    out = 0.0_dp
    do ii = 1, n1
      do jj = 2, n2 - 1
        out(ii,jj) = (f(ii,jj+1) - 2.0_dp*f(ii,jj) + f(ii,jj-1)) / dyval**2
      end do
    end do
  end subroutine second_derivative_y

  subroutine laplacian(f, dxval, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval, dyval
    real(dp), intent(out) :: out(n1,n2)
    integer :: ii, jj
    out = 0.0_dp
    do jj = 2, n2 - 1
      do ii = 2, n1 - 1
        out(ii,jj) = (f(ii+1,jj) - 2.0_dp*f(ii,jj) + f(ii-1,jj)) / dxval**2 &
                    + (f(ii,jj+1) - 2.0_dp*f(ii,jj) + f(ii,jj-1)) / dyval**2
      end do
    end do
  end subroutine laplacian

end program case12_euler_step_reference
