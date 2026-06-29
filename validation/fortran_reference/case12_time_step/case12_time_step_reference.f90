program case12_time_step_reference
  !! Standalone Fortran reference for one full Case 12 time-marching step
  !! (Etapa 10): explicit Euler step (Etapa 9) + boundary-condition
  !! re-imposition at t+dt + Poisson re-solve for psi (Etapa 8) + recompute
  !! of u, v from the updated psi.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It mirrors the sequence
  !!   call source_term(time + dt)
  !!   call boundaries_condictions(time + dt)
  !!   call loop(2)
  !! from subroutine solver_euler (src/fortran_original/4Order_LPOG_Euler/euler_rk.f),
  !! where subroutine loop(num) (main.f) is, in order: call filter_trid_x,
  !! call filter_trid_y, call poi_psi4a, call calc_ux, call calc_uy. As of
  !! Etapa 17, this reference uses the real pipeline: the 4th-order LPOG
  !! compact derivative operators (derivative_x/derivative_y/
  !! second_derivative_x/second_derivative_y below), the Lele compact filter
  !! (apply_lele_filter), and the FAS multigrid V-cycle Poisson solver
  !! (poi_psi4a and dependents, mirroring
  !! src/fortran_original/4Order_LPOG_Euler/poisson.f) -- not the 2nd-order
  !! placeholder and the dense direct Poisson solve used before Etapa 17.
  !!
  !! wz boundary condition: as documented in
  !! src/mms_viscoelastic/numerics/time_marching.py, the wz boundary is
  !! re-imposed from the exact manufactured solution at t+dt, instead of the
  !! 4th-order one-sided stencil used by
  !! boundarie_vorticity_cavity_flow.f. This is a deliberate deviation from
  !! the Fortran reference.
  !!
  !! The input CSV (generate_input_csv.py) provides the t=0 exact fields
  !! (u, v, wz, psi, txx, txy, tyy) and source terms (tfwz, tftxx, tftxy,
  !! tftyy) on a uniform grid covering [0,1] x [0,1].
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case12_time_step_reference.f90 -o case12_time_step_reference.exe
  !!   .\case12_time_step_reference.exe input.csv output.csv 9 9 1.0e-5
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: alphaf = 0.48_dp
  real(dp), parameter :: af = (11.0_dp + 10.0_dp * alphaf) / 16.0_dp
  real(dp), parameter :: bf = (15.0_dp + 34.0_dp * alphaf) / 64.0_dp
  real(dp), parameter :: cf = (-3.0_dp + 6.0_dp * alphaf) / 32.0_dp
  real(dp), parameter :: df = (1.0_dp - 2.0_dp * alphaf) / 64.0_dp
  integer, parameter :: msh = 4
  integer :: nx, ny, i, j
  real(dp) :: dx, dy, dt, dx2, dy2, t_next
  real(dp) :: Re, Wi, beta, epsilon_, xi, alphaG, one_minus_beta
  character(len=512) :: infile, outfile, arg

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
  real(dp), allocatable :: fc(:,:,:), sc(:,:,:), fct(:,:,:)
  real(dp), allocatable :: u_ex(:,:), v_ex(:,:), wz_ex(:,:), psi_ex(:,:)
  real(dp), allocatable :: txx_ex(:,:), txy_ex(:,:), tyy_ex(:,:)
  real(dp), allocatable :: psi_new(:,:), u_full(:,:), v_full(:,:), u_new(:,:), v_new(:,:)

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
  infile = "case12_time_step_input_9x9.csv"
  outfile = "case12_time_step_output_9x9_dt1e-05.csv"

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
  dx2 = dx * dx
  dy2 = dy * dy
  one_minus_beta = 1.0_dp - beta
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

  ! --- exact manufactured fields at t+dt, used for the boundary conditions ---
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

  ! --- re-impose Dirichlet boundary conditions at t+dt ---
  call apply_boundary(wz_new, wz_ex, nx, ny)
  call apply_boundary(txx_new, txx_ex, nx, ny)
  call apply_boundary(txy_new, txy_ex, nx, ny)
  call apply_boundary(tyy_new, tyy_ex, nx, ny)

  ! --- apply the Lele compact filter (x-pass then y-pass), interior only ---
  call apply_lele_filter(wz_new, nx, ny)
  call apply_lele_filter(txx_new, nx, ny)
  call apply_lele_filter(txy_new, nx, ny)
  call apply_lele_filter(tyy_new, nx, ny)

  ! --- re-solve laplacian(psi) = wz_new for the interior of psi (4th-order
  !     compact multigrid V-cycle, poi_psi4a) ---
  allocate(psi_new(nx,ny))
  psi_new = psi_ex
  allocate(fc(msh,nx,ny), sc(msh,nx,ny), fct(msh,nx,ny))
  call poi_psi4a(nx, ny, msh, dx2, dy2, wz_new, psi_new, fc, sc, fct)

  ! --- recompute u = d(psi)/dy, v = -d(psi)/dx from the updated psi ---
  allocate(u_full(nx,ny), v_full(nx,ny))
  call derivative_y(psi_new, dy, nx, ny, u_full)
  call derivative_x(psi_new, dx, nx, ny, v_full)
  v_full = -v_full

  allocate(u_new(nx,ny), v_new(nx,ny))
  u_new = u_full
  v_new = v_full
  call apply_boundary(u_new, u_ex, nx, ny)
  call apply_boundary(v_new, v_ex, nx, ny)

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,wz,txx,txy,tyy,psi,u,v"
  do j = 1, ny
    do i = 1, nx
      write(20, '(ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3)') &
        xg(i,j), ',', yg(i,j), ',', wz_new(i,j), ',', txx_new(i,j), ',', txy_new(i,j), ',', tyy_new(i,j), &
        ',', psi_new(i,j), ',', u_new(i,j), ',', v_new(i,j)
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

  !! Mirrors rhs_tridfx_wz/rhs_tridfx_T in src/fortran_original/4Order_LPOG_Euler/loop.f.
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

  !! Mirrors rhs_tridfy_wz/rhs_tridfy_T in src/fortran_original/4Order_LPOG_Euler/loop.f.
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

  !! Mirrors lhs_tridfx/lhs_tridfy in src/fortran_original/4Order_LPOG_Euler/loop.f.
  subroutine lhs_tridf(a, b, c, lmax)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: a(lmax), b(lmax), c(lmax)
    a = 0.0_dp
    b = 1.0_dp
    c = 0.0_dp
    a(4:lmax-3) = alphaf
    c(4:lmax-3) = alphaf
  end subroutine lhs_tridf

  !! Port of the relevant part of subroutine loop (main.f): call
  !! filter_trid_x followed by call filter_trid_y, applied to a single
  !! field. Only the interior is overwritten by the filtered value; the
  !! boundary is left untouched (matches
  !! wz(2:imax-1,2:jmax-1) = rhswz(2:imax-1,2:jmax-1) in the Fortran).
  subroutine apply_lele_filter(field, n1, n2)
    integer, intent(in) :: n1, n2
    real(dp), intent(inout) :: field(n1,n2)
    real(dp) :: rhs(n1,n2), a1(n1), b1(n1), c1(n1), a2(n2), b2(n2), c2(n2)

    call rhs_tridfx(rhs, field, n1, n2)
    call lhs_tridf(a1, b1, c1, n1)
    call tridx(a1, b1, c1, rhs, n1, n2)
    field(2:n1-1, 2:n2-1) = rhs(2:n1-1, 2:n2-1)

    call rhs_tridfy(rhs, field, n1, n2)
    call lhs_tridf(a2, b2, c2, n2)
    call tridy(a2, b2, c2, rhs, n1, n2)
    field(2:n1-1, 2:n2-1) = rhs(2:n1-1, 2:n2-1)
  end subroutine apply_lele_filter

  !! Verbatim port of poi_psi4a and dependents
  !! (src/fortran_original/4Order_LPOG_Euler/poisson.f), see
  !! validation/fortran_reference/case12_poisson_multigrid/ for the
  !! standalone, separately validated version of this solver.
  subroutine poi_psi4a(imax, jmax, msh_, dxx, dyy, wzf, psif, fcf, scf, fctf)
    integer, intent(in) :: imax, jmax, msh_
    real(dp), intent(in) :: dxx, dyy
    real(dp), intent(in) :: wzf(imax,jmax)
    real(dp), intent(inout) :: psif(imax,jmax)
    real(dp), intent(inout) :: fcf(msh_,imax,jmax), scf(msh_,imax,jmax), fctf(msh_,imax,jmax)

    logical :: continua
    integer :: itera, nptsx, nptsy, mesh, vcycles
    real(dp) :: dx2, dy2, err, erromax, dyypdxx

    erromax = 1.0e-06_dp
    nptsx = imax
    nptsy = jmax
    dx2 = dxx
    dy2 = dyy
    itera = 2

    call function_source(imax, jmax, msh_, wzf, psif, fcf, scf, fctf)
    continua = .true.

    do vcycles = 1, 81
      do mesh = 1, msh_ - 1
        dyypdxx = dy2 / dx2
        call sor(imax, jmax, msh_, fcf, scf, nptsx, nptsy, dy2, dyypdxx, itera, mesh)
        call defect(imax, jmax, msh_, fcf, scf, fctf, nptsx, nptsy, dx2, dy2, -1.0_dp, mesh)
        if (mesh == 1) then
          call verifica_poi(imax, jmax, msh_, fctf, err, erromax, continua)
          if (.not. continua) goto 22
        end if
        call restrict(imax, jmax, msh_, fcf, scf, fctf, nptsx, nptsy, dx2, dy2, mesh)
        dyypdxx = dy2 / dx2
        call defect(imax, jmax, msh_, fcf, scf, fctf, nptsx, nptsy, dx2, dy2, 1.0_dp, mesh+1)
      end do

      dyypdxx = dy2 / dx2
      call sor(imax, jmax, msh_, fcf, scf, nptsx, nptsy, dy2, dyypdxx, 100*itera, msh_)

      do mesh = msh_, 2, -1
        call correction(imax, jmax, msh_, fcf, fctf, nptsx, nptsy, mesh)
        call intpol(imax, jmax, msh_, fcf, nptsx, nptsy, dx2, dy2, mesh)
        dyypdxx = dy2 / dx2
        call sor(imax, jmax, msh_, fcf, scf, nptsx, nptsy, dy2, dyypdxx, 1, mesh-1)
      end do
    end do

22  continue
    psif(:,:) = fcf(1,:,:)

    if (vcycles > 80) then
      write(*,*) 'Vcycles = ', vcycles
      error stop
    end if
  end subroutine poi_psi4a

  subroutine verifica_poi(imax, jmax, msh_, fctf, err, erromax, continua)
    integer, intent(in) :: imax, jmax, msh_
    real(dp), intent(in) :: fctf(msh_,imax,jmax), erromax
    real(dp), intent(out) :: err
    logical, intent(inout) :: continua
    integer :: i, j
    err = 0.0_dp
    do j = 2, jmax - 1
      do i = 2, imax - 1
        err = max(err, abs(fctf(1,i,j)))
      end do
    end do
    if (err < erromax) continua = .false.
  end subroutine verifica_poi

  subroutine function_source(imax, jmax, msh_, wzf, psif, fcf, scf, fctf)
    integer, intent(in) :: imax, jmax, msh_
    real(dp), intent(in) :: wzf(imax,jmax), psif(imax,jmax)
    real(dp), intent(out) :: fcf(msh_,imax,jmax), scf(msh_,imax,jmax), fctf(msh_,imax,jmax)
    integer :: i, j
    fcf = 0.0_dp
    scf = 0.0_dp
    fctf = 0.0_dp
    do j = 2, jmax - 1
      do i = 2, imax - 1
        scf(1,i,j) = wzf(i,j-1) + wzf(i,j+1) + wzf(i+1,j) + wzf(i-1,j) + 8.0_dp*wzf(i,j)
      end do
    end do
    fcf(1,:,:) = psif(:,:)
  end subroutine function_source

  subroutine sor(imax, jmax, msh_, fcf, scf, nptsx, nptsy, dy2, dyypdxx, itera, mesh)
    integer, intent(in) :: imax, jmax, msh_, nptsx, nptsy, itera, mesh
    real(dp), intent(in) :: dy2, dyypdxx
    real(dp), intent(inout) :: fcf(msh_,imax,jmax)
    real(dp), intent(in) :: scf(msh_,imax,jmax)
    integer :: i, j, it
    real(dp) :: a(jmax), b(jmax), c(jmax), rhs(jmax), rf

    rf = 1.2_dp
    if (itera == 1) rf = 1.0_dp

    a(1) = 0.0_dp
    b(1) = 1.0_dp
    c(1) = 0.0_dp
    do j = 2, nptsy - 1
      a(j) = 10.0_dp - 2.0_dp*dyypdxx
      b(j) = -20.0_dp * (1.0_dp + dyypdxx)
      c(j) = 10.0_dp - 2.0_dp*dyypdxx
    end do
    a(nptsy) = 0.0_dp
    b(nptsy) = 1.0_dp
    c(nptsy) = 0.0_dp

    do it = 1, itera
      do i = 2, nptsx - 1
        rhs(1) = fcf(mesh,i,1)
        do j = 2, nptsy - 1
          rhs(j) = scf(mesh,i,j)*dy2 - ( (1.0_dp+dyypdxx)*(fcf(mesh,i+1,j+1) + fcf(mesh,i-1,j+1) &
                    + fcf(mesh,i+1,j-1) + fcf(mesh,i-1,j-1)) &
                    + (-2.0_dp+10.0_dp*dyypdxx)*(fcf(mesh,i+1,j) + fcf(mesh,i-1,j)) )
        end do
        rhs(nptsy) = fcf(mesh,i,nptsy)
        call tridv(rhs, nptsy, a, b, c, jmax)
        do j = 2, nptsy - 1
          fcf(mesh,i,j) = fcf(mesh,i,j) + rf * (rhs(j) - fcf(mesh,i,j))
        end do
      end do
    end do
  end subroutine sor

  subroutine defect(imax, jmax, msh_, fcf, scf, fctf, nptsx, nptsy, dx2, dy2, si, mesh)
    integer, intent(in) :: imax, jmax, msh_, nptsx, nptsy, mesh
    real(dp), intent(in) :: dx2, dy2, si
    real(dp), intent(in) :: fcf(msh_,imax,jmax)
    real(dp), intent(inout) :: scf(msh_,imax,jmax)
    real(dp), intent(inout) :: fctf(msh_,imax,jmax)
    integer :: i, j
    real(dp) :: b, c

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        b = ( fcf(mesh,i+1,j+1) + fcf(mesh,i+1,j-1) + fcf(mesh,i-1,j+1) + fcf(mesh,i-1,j-1) &
              + 10.0_dp*(fcf(mesh,i+1,j) + fcf(mesh,i-1,j)) &
              - 2.0_dp*(fcf(mesh,i,j+1) + fcf(mesh,i,j-1)) &
              - 20.0_dp*fcf(mesh,i,j) )
        c = ( fcf(mesh,i+1,j+1) + fcf(mesh,i+1,j-1) + fcf(mesh,i-1,j+1) + fcf(mesh,i-1,j-1) &
              + 10.0_dp*(fcf(mesh,i,j+1) + fcf(mesh,i,j-1)) &
              - 2.0_dp*(fcf(mesh,i+1,j) + fcf(mesh,i-1,j)) &
              - 20.0_dp*fcf(mesh,i,j) )
        fctf(1,i,j) = scf(mesh,i,j) + si * (b/dx2 + c/dy2)
      end do
    end do

    if (si < 0) return

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        scf(mesh,i,j) = fctf(1,i,j)
      end do
    end do
  end subroutine defect

  subroutine restrict(imax, jmax, msh_, fcf, scf, fctf, nptsx, nptsy, dx2, dy2, mesh)
    integer, intent(in) :: imax, jmax, msh_, mesh
    integer, intent(inout) :: nptsx, nptsy
    real(dp), intent(inout) :: dx2, dy2
    real(dp), intent(inout) :: fcf(msh_,imax,jmax)
    real(dp), intent(inout) :: scf(msh_,imax,jmax)
    real(dp), intent(inout) :: fctf(msh_,imax,jmax)
    integer :: i, j, fi, fj

    nptsx = (nptsx + 1) / 2
    nptsy = (nptsy + 1) / 2
    dx2 = 4.0_dp * dx2
    dy2 = 4.0_dp * dy2

    fj = 3
    do j = 2, nptsy - 1
      fi = 3
      do i = 2, nptsx - 1
        fcf(mesh+1,i,j) = fcf(mesh,fi,fj)
        fctf(mesh+1,i,j) = fcf(mesh+1,i,j)
        scf(mesh+1,i,j) = 6.25e-2_dp * ( 2.0_dp * ( fctf(1,fi-1,fj) + fctf(1,fi+1,fj) &
                          + fctf(1,fi,fj-1) + fctf(1,fi,fj+1) ) + 4.0_dp*fctf(1,fi,fj) &
                          + fctf(1,fi+1,fj+1) + fctf(1,fi+1,fj-1) &
                          + fctf(1,fi-1,fj+1) + fctf(1,fi-1,fj-1) )
        fi = fi + 2
      end do
      fj = fj + 2
    end do
  end subroutine restrict

  subroutine intpol(imax, jmax, msh_, fcf, nptsx, nptsy, dx2, dy2, mesh)
    integer, intent(in) :: imax, jmax, msh_, mesh
    integer, intent(inout) :: nptsx, nptsy
    real(dp), intent(inout) :: dx2, dy2
    real(dp), intent(inout) :: fcf(msh_,imax,jmax)
    integer :: i, j, fi, fj
    real(dp) :: temp(imax,jmax)

    dx2 = 0.25_dp * dx2
    dy2 = 0.25_dp * dy2

    temp = 0.0_dp
    fj = 1
    do j = 1, nptsy
      fi = 1
      do i = 1, nptsx
        temp(fi,fj) = fcf(mesh,i,j)
        fi = fi + 2
      end do
      fj = fj + 2
    end do

    nptsx = 2*nptsx - 1
    nptsy = 2*nptsy - 1

    do j = 1, nptsy, 2
      do i = 2, nptsx - 1, 2
        temp(i,j) = 0.5_dp * (temp(i-1,j) + temp(i+1,j))
      end do
    end do

    do j = 2, nptsy - 1, 2
      do i = 1, nptsx
        temp(i,j) = 0.5_dp * (temp(i,j+1) + temp(i,j-1))
      end do
    end do

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        fcf(mesh-1,i,j) = fcf(mesh-1,i,j) + temp(i,j)
      end do
    end do
  end subroutine intpol

  subroutine correction(imax, jmax, msh_, fcf, fctf, nptsx, nptsy, mesh)
    integer, intent(in) :: imax, jmax, msh_, nptsx, nptsy, mesh
    real(dp), intent(inout) :: fcf(msh_,imax,jmax)
    real(dp), intent(in) :: fctf(msh_,imax,jmax)
    integer :: i, j
    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        fcf(mesh,i,j) = fcf(mesh,i,j) - fctf(mesh,i,j)
      end do
    end do
  end subroutine correction

  subroutine tridv(r, lmax, a, b, c, jmax)
    integer, intent(in) :: lmax, jmax
    real(dp), intent(in) :: a(jmax), b(jmax), c(jmax)
    real(dp), intent(inout) :: r(jmax)
    integer :: j
    real(dp) :: bet, gam(jmax), u(jmax)
    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, lmax
      gam(j) = c(j-1) / bet
      bet = b(j) - a(j)*gam(j)
      u(j) = (r(j) - a(j)*u(j-1)) / bet
    end do
    do j = lmax - 1, 1, -1
      u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    r(1:lmax) = u(1:lmax)
  end subroutine tridv

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

end program case12_time_step_reference
