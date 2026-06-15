program case12_full_march_reference
  !! Standalone Fortran reference for the Case 12 full time march (Etapa 12):
  !! repeats the full time-marching step from Etapa 10/11 (explicit Euler
  !! step + boundary-condition re-imposition at t+dt + Poisson re-solve for
  !! psi + recompute of u, v from the updated psi) until the
  !! convergence-to-steady-state criterion of subroutine solver_euler
  !! (src/fortran_original/4Order_LPOG_Euler/euler_rk.f, lines 352-408) is
  !! met, instead of a fixed n_steps.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. Each iteration mirrors the sequence
  !!   call source_term(time + dt)
  !!   call boundaries_condictions(time + dt)
  !!   call loop(2)
  !! from subroutine solver_euler, using the same second-order
  !! finite-difference operators as mms_viscoelastic.numerics (NOT the
  !! 4th-order LPOG compact operators nor the Lele compact filter, which are
  !! out of scope, as in Etapas 8-11).
  !!
  !! wz boundary condition: as in case12_time_step_reference.f90, the wz
  !! boundary is re-imposed from the exact manufactured solution at t+dt
  !! instead of the 4th-order one-sided stencil used by
  !! boundarie_vorticity_cavity_flow.f.
  !!
  !! Stopping criterion (mirrors solver_euler)
  !! -------------------------------------------------
  !! At each step t (1-indexed), with time = (t-1)*dt before the update:
  !!   errowz  = max|wz_new  - wz_old | / max(max|wz_new |, 1)
  !!   erroTxx = max|Txx_new - Txx_old| / max(max|Txx_new|, 1)
  !!   erroTxy = max|Txy_new - Txy_old| / max(max|Txy_new|, 1)
  !!   erroTyy = max|Tyy_new - Tyy_old| / max(max|Tyy_new|, 1)
  !!   emv = max(errowz, erroTxx, erroTxy, erroTyy)
  !! The march stops as soon as (emv < 1e-7 .and. time > 0.01), or after
  !! max_steps iterations (safety cap; solver_euler's own upper bound,
  !! time < timef - dt with timef = 50, corresponds to ~5e6 steps at
  !! dt = 1e-5 and is not used directly as a practical cap here).
  !!
  !! Time-independence at the default "at" parameter
  !! -------------------------------------------------
  !! With the default manufactured-solution parameter at = 0 (used
  !! throughout this project), both the exact fields
  !! exact_solutions_case_12(t) and the manufactured source terms
  !! tfwz/tfTxx/tfTxy/tfTyy are independent of t. The boundary values and
  !! source terms read from the t=0 input CSV are therefore reused unchanged
  !! at every step, as in case12_multistep_reference.f90.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case12_full_march_reference.f90 -o case12_full_march_reference.exe
  !!   .\case12_full_march_reference.exe input.csv output.csv 9 9 1.0e-5 history.csv 200000 50.0 [hist_start_step] [hist_stride]
  !!
  !! hist_start_step (default 1) and hist_stride (default 1) restrict the
  !! history file to rows with t_step >= hist_start_step and
  !! (t_step - hist_start_step) a multiple of hist_stride, keeping the file
  !! small for long runs. The history file also reports the per-step maximum
  !! absolute values of wz, Txx, Txy, Tyy and psi, used to detect blow-up.
  !!
  !! Divergence detection
  !! -------------------------------------------------
  !! If any of errowz/erroTxx/erroTxy/erroTyy/emv becomes NaN (e.g. because a
  !! field overflowed to +/-Infinity), the march stops immediately with
  !! diverged = .true.. Without this guard, gfortran's max(NaN, x) = x
  !! semantics make emv collapse to 0.0 once a NaN appears, which would be
  !! misreported as converged = .true..
  use iso_fortran_env, only: real64, int64
  use ieee_arithmetic, only: ieee_is_nan
  implicit none

  integer, parameter :: dp = real64
  integer :: nx, ny, i, j, ix, jy, nix, niy, n, k
  integer(int64) :: t_step, max_steps, hist_start_step, hist_stride
  integer(int64) :: clock_start, clock_end, clock_rate
  real(dp) :: dx, dy, dt, dx2, dy2, t_next, bterm, timef
  real(dp) :: time, errowz, erroTxx, erroTxy, erroTyy, emv
  real(dp) :: Re, Wi, beta, epsilon_, xi, alphaG, one_minus_beta
  logical :: converged, diverged
  character(len=512) :: infile, outfile, histfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:)
  real(dp), allocatable :: ux(:,:), uy(:,:), wz(:,:), psi(:,:), txx(:,:), txy(:,:), tyy(:,:)
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
  real(dp), allocatable :: u_ex(:,:), v_ex(:,:), wz_ex(:,:), psi_ex(:,:)
  real(dp), allocatable :: txx_ex(:,:), txy_ex(:,:), tyy_ex(:,:)
  real(dp), allocatable :: psi_new(:,:), u_full(:,:), v_full(:,:), u_new(:,:), v_new(:,:)
  real(dp), allocatable :: amat(:,:), bvec(:), xvec(:)

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
  max_steps = 200000_int64
  timef = 50.0_dp
  hist_start_step = 1_int64
  hist_stride = 1_int64
  infile = "case12_time_step_input_9x9.csv"
  outfile = "case12_full_march_output_9x9.csv"
  histfile = ""

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

  call get_command_argument(6, arg)
  if (len_trim(arg) > 0) histfile = trim(arg)

  call get_command_argument(7, arg)
  if (len_trim(arg) > 0) read(arg, *) max_steps

  call get_command_argument(8, arg)
  if (len_trim(arg) > 0) read(arg, *) timef

  call get_command_argument(9, arg)
  if (len_trim(arg) > 0) read(arg, *) hist_start_step

  call get_command_argument(10, arg)
  if (len_trim(arg) > 0) read(arg, *) hist_stride

  if (nx <= 2 .or. ny <= 2) error stop "nx and ny must be greater than 2"
  if (max_steps < 1) error stop "max_steps must be positive"

  dx = 1.0_dp / real(nx - 1, dp)
  dy = 1.0_dp / real(ny - 1, dp)
  dx2 = dx * dx
  dy2 = dy * dy
  one_minus_beta = 1.0_dp - beta
  ! Exact boundary fields and source terms are independent of t for the
  ! default at = 0 (see module-level note above); t_next is only used to
  ! evaluate the (constant) exact-field functions.
  t_next = dt

  allocate(xg(nx,ny), yg(nx,ny))
  allocate(ux(nx,ny), uy(nx,ny), wz(nx,ny), psi(nx,ny), txx(nx,ny), txy(nx,ny), tyy(nx,ny))
  allocate(tfwz(nx,ny), tftxx(nx,ny), tftxy(nx,ny), tftyy(nx,ny))

  open(unit=10, file=trim(infile), status="old", action="read")
  read(10, *)  ! skip header
  do j = 1, ny
    do i = 1, nx
      read(10, *) xg(i,j), yg(i,j), ux(i,j), uy(i,j), wz(i,j), psi(i,j), txx(i,j), txy(i,j), tyy(i,j), &
                  tfwz(i,j), tftxx(i,j), tftxy(i,j), tftyy(i,j)
    end do
  end do
  close(10)

  ! --- exact manufactured fields, used for the boundary conditions at every step ---
  allocate(u_ex(nx,ny), v_ex(nx,ny), wz_ex(nx,ny), psi_ex(nx,ny))
  allocate(txx_ex(nx,ny), txy_ex(nx,ny), tyy_ex(nx,ny))
  do j = 1, ny
    do i = 1, nx
      u_ex(i,j) = u_case12_default(xg(i,j), yg(i,j), t_next)
      v_ex(i,j) = v_case12_default(xg(i,j), yg(i,j), t_next)
      wz_ex(i,j) = wz_case12_default(xg(i,j), yg(i,j), t_next)
      psi_ex(i,j) = psi_case12_default(xg(i,j), yg(i,j), t_next)
      txx_ex(i,j) = txx_case12_default(xg(i,j), yg(i,j), t_next)
      txy_ex(i,j) = txy_case12_default(xg(i,j), yg(i,j), t_next)
      tyy_ex(i,j) = tyy_case12_default(xg(i,j), yg(i,j), t_next)
    end do
  end do

  allocate(dudx(nx,ny), dudy(nx,ny), dvdx(nx,ny), dvdy(nx,ny))
  allocate(duwzdx(nx,ny), dvwzdy(nx,ny), diffusion(nx,ny))
  allocate(dtxxdx(nx,ny), d2txx_dxdy(nx,ny))
  allocate(d2txy_dx2(nx,ny), d2txy_dy2(nx,ny))
  allocate(dtyydx(nx,ny), d2tyy_dxdy(nx,ny))
  allocate(dvz(nx,ny))
  allocate(duTxxdx(nx,ny), dvTxxdy(nx,ny))
  allocate(duTxydx(nx,ny), dvTxydy(nx,ny))
  allocate(duTyydx(nx,ny), dvTyydy(nx,ny))
  allocate(f_tij(nx,ny))
  allocate(dvTxx(nx,ny), dvTxy(nx,ny), dvTyy(nx,ny))
  allocate(wz_new(nx,ny), txx_new(nx,ny), txy_new(nx,ny), tyy_new(nx,ny))
  allocate(psi_new(nx,ny), u_full(nx,ny), v_full(nx,ny), u_new(nx,ny), v_new(nx,ny))

  nix = nx - 2
  niy = ny - 2
  n = nix * niy
  allocate(amat(n,n), bvec(n), xvec(n))

  if (len_trim(histfile) > 0) then
    open(unit=30, file=trim(histfile), status="replace", action="write")
    write(30, '(A)') "step,time,errowz,erroTxx,erroTxy,erroTyy,emv,maxwz,maxtxx,maxtxy,maxtyy,maxpsi"
  end if

  converged = .false.
  diverged = .false.
  call system_clock(count_rate=clock_rate)
  call system_clock(count=clock_start)

  do t_step = 1, max_steps
    time = real(t_step - 1, dp) * dt

    call derivative_x(ux, dx, nx, ny, dudx)
    call derivative_y(ux, dy, nx, ny, dudy)
    call derivative_x(uy, dx, nx, ny, dvdx)
    call derivative_y(uy, dy, nx, ny, dvdy)

    ! --- dvz: right-hand side of the wz transport equation ---
    call derivative_x(ux*wz, dx, nx, ny, duwzdx)
    call derivative_y(uy*wz, dy, nx, ny, dvwzdy)
    call laplacian(wz, dx, dy, nx, ny, diffusion)

    call derivative_x(txx, dx, nx, ny, dtxxdx)
    call derivative_y(dtxxdx, dy, nx, ny, d2txx_dxdy)

    call second_derivative_x(txy, dx, nx, ny, d2txy_dx2)
    call second_derivative_y(txy, dy, nx, ny, d2txy_dy2)

    call derivative_x(tyy, dx, nx, ny, dtyydx)
    call derivative_y(dtyydx, dy, nx, ny, d2tyy_dxdy)

    dvz = -duwzdx - dvwzdy + (beta/Re)*diffusion &
        + d2txx_dxdy + d2txy_dy2 - d2txy_dx2 - d2tyy_dxdy - tfwz

    ! --- dvTxx, dvTxy, dvTyy: right-hand sides of the stress equations ---
    call derivative_x(ux*txx, dx, nx, ny, duTxxdx)
    call derivative_y(uy*txx, dy, nx, ny, dvTxxdy)

    call derivative_x(ux*txy, dx, nx, ny, duTxydx)
    call derivative_y(uy*txy, dy, nx, ny, dvTxydy)

    call derivative_x(ux*tyy, dx, nx, ny, duTyydx)
    call derivative_y(uy*tyy, dy, nx, ny, dvTyydy)

    f_tij = 1.0_dp + (epsilon_*Re*Wi/one_minus_beta)*(txx + tyy)

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

    ! --- re-impose Dirichlet boundary conditions at t+dt ---
    call apply_boundary(wz_new, wz_ex, nx, ny)
    call apply_boundary(txx_new, txx_ex, nx, ny)
    call apply_boundary(txy_new, txy_ex, nx, ny)
    call apply_boundary(tyy_new, tyy_ex, nx, ny)

    ! --- re-solve laplacian(psi) = wz_new for the interior of psi ---
    amat = 0.0_dp
    do jy = 1, niy
      do ix = 1, nix
        k = (jy - 1) * nix + ix
        i = ix + 1
        j = jy + 1
        amat(k,k) = -2.0_dp / dx2 - 2.0_dp / dy2
        bterm = 0.0_dp

        if (ix > 1) then
          amat(k, k-1) = amat(k, k-1) + 1.0_dp / dx2
        else
          bterm = bterm + psi_ex(1, j) / dx2
        end if

        if (ix < nix) then
          amat(k, k+1) = amat(k, k+1) + 1.0_dp / dx2
        else
          bterm = bterm + psi_ex(nx, j) / dx2
        end if

        if (jy > 1) then
          amat(k, k-nix) = amat(k, k-nix) + 1.0_dp / dy2
        else
          bterm = bterm + psi_ex(i, 1) / dy2
        end if

        if (jy < niy) then
          amat(k, k+nix) = amat(k, k+nix) + 1.0_dp / dy2
        else
          bterm = bterm + psi_ex(i, ny) / dy2
        end if

        bvec(k) = wz_new(i,j) - bterm
      end do
    end do

    call solve_linear_system(amat, bvec, xvec, n)

    psi_new = psi_ex
    do jy = 1, niy
      do ix = 1, nix
        k = (jy - 1) * nix + ix
        psi_new(ix+1, jy+1) = xvec(k)
      end do
    end do

    ! --- recompute u = d(psi)/dy, v = -d(psi)/dx from the updated psi ---
    call derivative_y(psi_new, dy, nx, ny, u_full)
    call derivative_x(psi_new, dx, nx, ny, v_full)
    v_full = -v_full

    u_new = u_full
    v_new = v_full
    call apply_boundary(u_new, u_ex, nx, ny)
    call apply_boundary(v_new, v_ex, nx, ny)

    ! --- convergence-to-steady-state criterion (solver_euler, lines 373-401) ---
    errowz  = maxval(abs(wz_new  - wz )) / max(maxval(abs(wz_new )), 1.0_dp)
    erroTxx = maxval(abs(txx_new - txx)) / max(maxval(abs(txx_new)), 1.0_dp)
    erroTxy = maxval(abs(txy_new - txy)) / max(maxval(abs(txy_new)), 1.0_dp)
    erroTyy = maxval(abs(tyy_new - tyy)) / max(maxval(abs(tyy_new)), 1.0_dp)
    emv = max(errowz, max(erroTxx, max(erroTxy, erroTyy)))

    if (len_trim(histfile) > 0 .and. t_step >= hist_start_step &
        .and. mod(t_step - hist_start_step, hist_stride) == 0_int64) then
      write(30, '(I0,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,' // &
                 'A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3)') &
        t_step, ',', time + dt, ',', errowz, ',', erroTxx, ',', erroTxy, ',', erroTyy, ',', emv, &
        ',', maxval(abs(wz_new)), ',', maxval(abs(txx_new)), ',', maxval(abs(txy_new)), &
        ',', maxval(abs(tyy_new)), ',', maxval(abs(psi_new))
    end if

    ! --- state at t + dt becomes the state at the start of the next step ---
    ux = u_new
    uy = v_new
    wz = wz_new
    txx = txx_new
    txy = txy_new
    tyy = tyy_new
    psi = psi_new

    if (ieee_is_nan(emv) .or. ieee_is_nan(errowz) .or. ieee_is_nan(erroTxx) &
        .or. ieee_is_nan(erroTxy) .or. ieee_is_nan(erroTyy)) then
      diverged = .true.
      exit
    end if

    if ((emv < 1.0e-7_dp) .and. (time > 0.01_dp)) then
      converged = .true.
      exit
    end if

    if (time >= timef - dt) exit
  end do

  call system_clock(count=clock_end)

  if (len_trim(histfile) > 0) close(30)

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,wz,txx,txy,tyy,psi,u,v"
  do j = 1, ny
    do i = 1, nx
      write(20, '(ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3)') &
        xg(i,j), ',', yg(i,j), ',', wz(i,j), ',', txx(i,j), ',', txy(i,j), ',', tyy(i,j), &
        ',', psi(i,j), ',', ux(i,j), ',', uy(i,j)
    end do
  end do
  close(20)

  write(*, '(A,I0,A,I0)') "grid = ", nx, " x ", ny
  write(*, '(A,L1)') "converged = ", converged
  write(*, '(A,L1)') "diverged = ", diverged
  write(*, '(A,I0)') "steps = ", t_step
  write(*, '(A,ES12.5)') "final time = ", time + dt
  write(*, '(A,ES12.5)') "final emv = ", emv
  write(*, '(A,F10.3,A)') "wall time = ", real(clock_end - clock_start, dp) / real(clock_rate, dp), " s"

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

  subroutine apply_boundary(field, bvals, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(inout) :: field(n1,n2)
    real(dp), intent(in) :: bvals(n1,n2)
    integer :: ii, jj
    do jj = 1, n2
      field(1,jj) = bvals(1,jj)
      field(n1,jj) = bvals(n1,jj)
    end do
    do ii = 1, n1
      field(ii,1) = bvals(ii,1)
      field(ii,n2) = bvals(ii,n2)
    end do
  end subroutine apply_boundary

  !! Dense Gaussian elimination with partial pivoting, A x = b. Suitable for
  !! the small interior Poisson system assembled above (n = (nx-2)*(ny-2)).
  subroutine solve_linear_system(a, b, x, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n,n), b(n)
    real(dp), intent(out) :: x(n)
    real(dp) :: aug(n,n+1), factor, pivot_val, tmp
    integer :: ii, jj, kk, piv

    aug(:,1:n) = a
    aug(:,n+1) = b

    do kk = 1, n
      piv = kk
      pivot_val = abs(aug(kk,kk))
      do ii = kk+1, n
        if (abs(aug(ii,kk)) > pivot_val) then
          pivot_val = abs(aug(ii,kk))
          piv = ii
        end if
      end do
      if (piv /= kk) then
        do jj = 1, n+1
          tmp = aug(kk,jj)
          aug(kk,jj) = aug(piv,jj)
          aug(piv,jj) = tmp
        end do
      end if
      do ii = kk+1, n
        factor = aug(ii,kk) / aug(kk,kk)
        aug(ii,kk:n+1) = aug(ii,kk:n+1) - factor * aug(kk,kk:n+1)
      end do
    end do

    do ii = n, 1, -1
      x(ii) = aug(ii,n+1)
      do jj = ii+1, n
        x(ii) = x(ii) - aug(ii,jj) * x(jj)
      end do
      x(ii) = x(ii) / aug(ii,ii)
    end do
  end subroutine solve_linear_system

  pure function u_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp

    value = ((-1.0_dp + x) * x * sin(pii * x) &
        * (pii * y * cos(pii * y) + sin(pii * y))) / exp(at * t)
  end function u_case12_default

  pure function v_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp

    value = -((y * (pii * (-1.0_dp + x) * x * cos(pii * x) &
        + (-1.0_dp + 2.0_dp * x) * sin(pii * x)) * sin(pii * y)) / exp(at * t))
  end function v_case12_default

  pure function wz_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp

    value = (2.0_dp * (pii * (-1.0_dp + x) * x * cos(pii * y) * sin(pii * x) &
        + y * (pii * (-1.0_dp + 2.0_dp * x) * cos(pii * x) &
        + (1.0_dp - pii**2 * (-1.0_dp + x) * x) * sin(pii * x)) * sin(pii * y))) &
        / exp(at * t)
  end function wz_case12_default

  pure function psi_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp

    value = ((-1.0_dp + x) * x * y * sin(pii * x) * sin(pii * y)) / exp(at * t)
  end function psi_case12_default

  pure function txx_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: betann = 0.5_dp

    real(dp), parameter :: a1xx = 1.0_dp
    real(dp), parameter :: a2xx = 1.0_dp
    real(dp), parameter :: a3xx = 1.0_dp
    real(dp), parameter :: b1xx = 1.0_dp
    real(dp), parameter :: b2xx = 1.0_dp
    real(dp), parameter :: b3xx = 1.0_dp

    value = ((1.0_dp - betann) * (a1xx + a3xx * cos(pii * x) + a2xx * sin(pii * x)) &
        * (b1xx + b3xx * cos(pii * y) + b2xx * sin(pii * y))) / exp(at * t)
  end function txx_case12_default

  pure function txy_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: betann = 0.5_dp

    real(dp), parameter :: a1xy = 1.0_dp
    real(dp), parameter :: a2xy = 1.0_dp
    real(dp), parameter :: a3xy = 1.0_dp
    real(dp), parameter :: b1xy = 1.0_dp
    real(dp), parameter :: b2xy = 1.0_dp
    real(dp), parameter :: b3xy = 1.0_dp

    value = ((1.0_dp - betann) * (a1xy + a3xy * cos(pii * x) + a2xy * sin(pii * x)) &
        * (b1xy + b3xy * cos(pii * y) + b2xy * sin(pii * y))) / exp(at * t)
  end function txy_case12_default

  pure function tyy_case12_default(x, y, t) result(value)
    real(dp), intent(in) :: x, y, t
    real(dp) :: value

    real(dp), parameter :: pii = acos(-1.0_dp)
    real(dp), parameter :: at = 0.0_dp
    real(dp), parameter :: betann = 0.5_dp

    real(dp), parameter :: a1yy = 1.0_dp
    real(dp), parameter :: a2yy = 1.0_dp
    real(dp), parameter :: a3yy = 1.0_dp
    real(dp), parameter :: b1yy = 1.0_dp
    real(dp), parameter :: b2yy = 1.0_dp
    real(dp), parameter :: b3yy = 1.0_dp

    value = ((1.0_dp - betann) * (a1yy + a3yy * cos(pii * x) + a2yy * sin(pii * x)) &
        * (b1yy + b3yy * cos(pii * y) + b2yy * sin(pii * y))) / exp(at * t)
  end function tyy_case12_default

end program case12_full_march_reference
