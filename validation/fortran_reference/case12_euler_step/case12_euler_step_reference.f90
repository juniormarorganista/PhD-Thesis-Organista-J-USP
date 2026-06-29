program case12_euler_step_reference
  !! Standalone Fortran reference for one explicit Euler step of the Case 12
  !! ``ty_sim == 1`` formulation (wz, Txx, Txy, Tyy).
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It is a verbatim port of the right-hand sides computed by
  !! ``subroutine drv`` (src/fortran_original/4Order_LPOG_Euler/drvs.f) and the
  !! explicit Euler update performed by ``subroutine solver_euler``
  !! (src/fortran_original/4Order_LPOG_Euler/euler_rk.f), using the same
  !! 4th-order LPOG compact (Pade-type) derivative operators that ``drv``
  !! itself uses (``derx``/``dery``/``derxx``/``deryy`` in
  !! src/fortran_original/4Order_LPOG_Euler/derivs.f are the compact
  !! operators; ``derivative_x``/``derivative_y``/``second_derivative_x``/
  !! ``second_derivative_y`` below reproduce them, replacing the 2nd-order
  !! placeholder used before Etapa 17). The Lele compact filter is a
  !! separate per-time-step operation (``subroutine loop``, ``main.f``), not
  !! part of these right-hand sides, and is out of scope for a single Euler
  !! step in isolation (see ``case12_time_step`` for where it applies).
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

  ! 4th-order LPOG compact (Pade-type) operators (Etapa 17), replacing the
  ! 2nd-order placeholder used before. Same name and signature as before, so
  ! every call site above is unchanged.
  subroutine derivative_x(f, dxval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval
    real(dp), intent(out) :: out(n1,n2)
    real(dp) :: a(n1), b(n1), c(n1)
    call coef_4a(a, b, c, n1)
    call rhsx_4a(out, f, n1, n2, dxval)
    call tridx(a, b, c, out, n1, n2)
  end subroutine derivative_x

  subroutine derivative_y(f, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dyval
    real(dp), intent(out) :: out(n1,n2)
    real(dp) :: a(n2), b(n2), c(n2)
    call coef_4a(a, b, c, n2)
    call rhsy_4a(out, f, n1, n2, dyval)
    call tridy(a, b, c, out, n1, n2)
  end subroutine derivative_y

  subroutine second_derivative_x(f, dxval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval
    real(dp), intent(out) :: out(n1,n2)
    real(dp) :: a(n1), b(n1), c(n1)
    call coeff_4a(a, b, c, n1)
    call rhsxx_4a(out, f, n1, n2, dxval)
    call tridx(a, b, c, out, n1, n2)
  end subroutine second_derivative_x

  subroutine second_derivative_y(f, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dyval
    real(dp), intent(out) :: out(n1,n2)
    real(dp) :: a(n2), b(n2), c(n2)
    call coeff_4a(a, b, c, n2)
    call rhsyy_4a(out, f, n1, n2, dyval)
    call tridy(a, b, c, out, n1, n2)
  end subroutine second_derivative_y

  subroutine laplacian(f, dxval, dyval, n1, n2, out)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: f(n1,n2), dxval, dyval
    real(dp), intent(out) :: out(n1,n2)
    real(dp) :: dxx_term(n1,n2), dyy_term(n1,n2)
    call second_derivative_x(f, dxval, n1, n2, dxx_term)
    call second_derivative_y(f, dyval, n1, n2, dyy_term)
    out = dxx_term + dyy_term
  end subroutine laplacian

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

end program case12_euler_step_reference
