program case12_log_euler_step_reference
  !! Standalone Fortran reference for one explicit Euler step of the Case 12
  !! log-conformation formulation (ty_sim == 22): wz, Psixx, Psixy, Psiyy.
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It is a verbatim port of:
  !!   - decomp_matrix, multi_matrix, T_to_Psi, Psi_to_T, kernel_conformation
  !!     (src/fortran_original/4Order_LPOG_Euler/conform.f), specialized to
  !!     dimen = 2 (the only case ever used by the original code), but kept
  !!     as the same Jacobi-rotation algorithm rather than a closed form, so
  !!     that this reference is independent from the Python closed-form
  !!     implementation it validates.
  !!   - drv_Psi, nlterms_logPsi (src/fortran_original/4Order_LPOG_Euler/drvs.f)
  !!   - the explicit Euler update in subroutine solver_log_euler
  !!     (src/fortran_original/4Order_LPOG_Euler/euler_rk.f)
  !!
  !! using the same 4th-order LPOG compact derivative operators that drv_Psi
  !! itself uses (derx/dery in drvs.f are the compact operators from
  !! derivs.f; derivative_x/derivative_y/second_derivative_x/
  !! second_derivative_y below reproduce them, replacing the 2nd-order
  !! placeholder used before Etapa 17, mirroring case12_euler_step). The
  !! Lele compact filter is a separate per-time-step operation, out of
  !! scope for a single Euler step in isolation.
  !!
  !! The input CSV provides u, v, wz, txx, txy, tyy and the manufactured
  !! source terms tfwz, tftxx, tftxy, tftyy on a uniform grid covering
  !! [0,1] x [0,1]. The polymeric stress is converted to the log-conformation
  !! representation Psi before marching (T_to_Psi), advanced one Euler step,
  !! and the output CSV reports wz and the recovered (Txx, Txy, Tyy)
  !! (Psi_to_T of the new Psi) for direct comparison with the non-log
  !! reference's output schema. Only interior grid points are advanced;
  !! boundary values are copied unchanged.
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case12_log_euler_step_reference.f90 -o case12_log_euler_step_reference.exe
  !!   .\case12_log_euler_step_reference.exe input.csv output.csv 9 9 1.0e-5
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: dimen = 2
  integer :: nx, ny, i, j
  real(dp) :: dx, dy, dt
  real(dp) :: Re, Wi, beta, alphaG, etap, etapinv
  character(len=512) :: infile, outfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:)
  real(dp), allocatable :: ux(:,:), uy(:,:), wz(:,:), txx(:,:), txy(:,:), tyy(:,:)
  real(dp), allocatable :: psixx(:,:), psixy(:,:), psiyy(:,:)
  real(dp), allocatable :: tfwz(:,:), tftxx(:,:), tftxy(:,:), tftyy(:,:)
  real(dp), allocatable :: tfpsixx(:,:), tfpsixy(:,:), tfpsiyy(:,:)
  real(dp), allocatable :: dudx(:,:), dudy(:,:), dvdx(:,:), dvdy(:,:)
  real(dp), allocatable :: duwzdx(:,:), dvwzdy(:,:), diffusion(:,:)
  real(dp), allocatable :: dtxxdx(:,:), d2txx_dxdy(:,:)
  real(dp), allocatable :: d2txy_dx2(:,:), d2txy_dy2(:,:)
  real(dp), allocatable :: dtyydx(:,:), d2tyy_dxdy(:,:)
  real(dp), allocatable :: omegalog(:,:), bxx(:,:), bxy(:,:), byy(:,:)
  real(dp), allocatable :: mxx(:,:), mxy(:,:), myy(:,:)
  real(dp), allocatable :: uPsixx(:,:), vPsixx(:,:), uPsixy(:,:), vPsixy(:,:)
  real(dp), allocatable :: uPsiyy(:,:), vPsiyy(:,:)
  real(dp), allocatable :: duPsixxdx(:,:), dvPsixxdy(:,:)
  real(dp), allocatable :: duPsixydx(:,:), dvPsixydy(:,:)
  real(dp), allocatable :: duPsiyydx(:,:), dvPsiyydy(:,:)
  real(dp), allocatable :: dvz(:,:), dvPsixx(:,:), dvPsixy(:,:), dvPsiyy(:,:)
  real(dp), allocatable :: wz_new(:,:), psixx_new(:,:), psixy_new(:,:), psiyy_new(:,:)
  real(dp), allocatable :: txx_out(:,:), txy_out(:,:), tyy_out(:,:)

  ! Defaults match EulerStepParameters() in mms_viscoelastic.numerics.time_stepping.
  nx = 9
  ny = 9
  Re = 1.0_dp
  Wi = 1.0_dp
  beta = 0.5_dp
  alphaG = 0.0_dp
  dt = 1.0e-5_dp
  infile = "case12_log_euler_step_input_9x9.csv"
  outfile = "case12_log_euler_step_output_9x9_dt1e-05.csv"

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
  etap = (1.0_dp - beta) / (Re * Wi)
  etapinv = 1.0_dp / etap

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

  allocate(psixx(nx,ny), psixy(nx,ny), psiyy(nx,ny))
  allocate(tfpsixx(nx,ny), tfpsixy(nx,ny), tfpsiyy(nx,ny))
  tfpsixx = tftxx * etapinv
  tfpsixy = tftxy * etapinv
  tfpsiyy = tftyy * etapinv

  do j = 1, ny
    do i = 1, nx
      call T_to_Psi_point(txx(i,j), txy(i,j), tyy(i,j), etap, &
                           psixx(i,j), psixy(i,j), psiyy(i,j))
    end do
  end do

  allocate(dudx(nx,ny), dudy(nx,ny), dvdx(nx,ny), dvdy(nx,ny))
  call derivative_x(ux, dx, nx, ny, dudx)
  call derivative_y(ux, dy, nx, ny, dudy)
  call derivative_x(uy, dx, nx, ny, dvdx)
  call derivative_y(uy, dy, nx, ny, dvdy)

  ! --- recover Txx, Txy, Tyy from Psi (call Psi_to_T) ---
  do j = 1, ny
    do i = 1, nx
      call Psi_to_T_point(psixx(i,j), psixy(i,j), psiyy(i,j), etap, &
                           txx(i,j), txy(i,j), tyy(i,j))
    end do
  end do

  ! --- kernel_conformation(Omegalog, Bxx, Bxy, Byy, Mxx, Mxy, Myy) ---
  allocate(omegalog(nx,ny), bxx(nx,ny), bxy(nx,ny), byy(nx,ny))
  allocate(mxx(nx,ny), mxy(nx,ny), myy(nx,ny))
  do j = 1, ny
    do i = 1, nx
      call kernel_conformation_point(txx(i,j), txy(i,j), tyy(i,j), &
                                      dudx(i,j), dudy(i,j), dvdx(i,j), dvdy(i,j), &
                                      etapinv, alphaG, &
                                      omegalog(i,j), bxx(i,j), bxy(i,j), byy(i,j), &
                                      mxx(i,j), mxy(i,j), myy(i,j))
    end do
  end do

  ! --- dvz: right-hand side of the wz transport equation (uses recovered Txx, Txy, Tyy) ---
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

  ! --- nlterms_logPsi + dvPsixx, dvPsixy, dvPsiyy ---
  allocate(uPsixx(nx,ny), vPsixx(nx,ny), uPsixy(nx,ny), vPsixy(nx,ny))
  allocate(uPsiyy(nx,ny), vPsiyy(nx,ny))
  uPsixx = ux * psixx
  vPsixx = uy * psixx
  uPsixy = ux * psixy
  vPsixy = uy * psixy
  uPsiyy = ux * psiyy
  vPsiyy = uy * psiyy

  allocate(duPsixxdx(nx,ny), dvPsixxdy(nx,ny))
  call derivative_x(uPsixx, dx, nx, ny, duPsixxdx)
  call derivative_y(vPsixx, dy, nx, ny, dvPsixxdy)

  allocate(duPsixydx(nx,ny), dvPsixydy(nx,ny))
  call derivative_x(uPsixy, dx, nx, ny, duPsixydx)
  call derivative_y(vPsixy, dy, nx, ny, dvPsixydy)

  allocate(duPsiyydx(nx,ny), dvPsiyydy(nx,ny))
  call derivative_x(uPsiyy, dx, nx, ny, duPsiyydx)
  call derivative_y(vPsiyy, dy, nx, ny, dvPsiyydy)

  allocate(dvPsixx(nx,ny), dvPsixy(nx,ny), dvPsiyy(nx,ny))
  dvPsixx = -duPsixxdx - dvPsixxdy + 2.0_dp*omegalog*psixy + 2.0_dp*bxx + mxx/Wi - tfpsixx
  dvPsixy = -duPsixydx - dvPsixydy + omegalog*(psiyy - psixx) + 2.0_dp*bxy + mxy/Wi - tfpsixy
  dvPsiyy = -duPsiyydx - dvPsiyydy - 2.0_dp*omegalog*psixy + 2.0_dp*byy + myy/Wi - tfpsiyy

  ! --- explicit Euler update, interior points only ---
  allocate(wz_new(nx,ny), psixx_new(nx,ny), psixy_new(nx,ny), psiyy_new(nx,ny))
  wz_new = wz
  psixx_new = psixx
  psixy_new = psixy
  psiyy_new = psiyy

  do j = 2, ny - 1
    do i = 2, nx - 1
      wz_new(i,j) = wz(i,j) + dvz(i,j) * dt
      psixx_new(i,j) = psixx(i,j) + dvPsixx(i,j) * dt
      psixy_new(i,j) = psixy(i,j) + dvPsixy(i,j) * dt
      psiyy_new(i,j) = psiyy(i,j) + dvPsiyy(i,j) * dt
    end do
  end do

  allocate(txx_out(nx,ny), txy_out(nx,ny), tyy_out(nx,ny))
  do j = 1, ny
    do i = 1, nx
      call Psi_to_T_point(psixx_new(i,j), psixy_new(i,j), psiyy_new(i,j), etap, &
                           txx_out(i,j), txy_out(i,j), tyy_out(i,j))
    end do
  end do

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,wz,txx,txy,tyy"
  do j = 1, ny
    do i = 1, nx
      write(20, '(ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3,A,ES26.17E3)') &
        xg(i,j), ',', yg(i,j), ',', wz_new(i,j), ',', txx_out(i,j), ',', txy_out(i,j), ',', tyy_out(i,j)
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

  ! Verbatim port of subroutine decomp_matrix (conform.f), specialized to a
  ! single dimen=2 matrix (no grid loop, no par.nn include).
  subroutine decomp_matrix(a, eigval, eigvec)
    real(dp), intent(inout) :: a(dimen,dimen)
    real(dp), intent(out) :: eigval(dimen), eigvec(dimen,dimen)
    integer :: ii, ip, iq, jj
    integer :: nrot
    real(dp) :: c, g, h, s, sm, t, tau, theta, tresh
    real(dp) :: b(dimen), z(dimen)

    eigval = 0.0_dp
    eigvec = 0.0_dp
    do ip = 1, dimen
      eigvec(ip,ip) = 1.0_dp
    end do

    do ip = 1, dimen
      b(ip) = a(ip,ip)
      eigval(ip) = b(ip)
      z(ip) = 0.0_dp
    end do

    nrot = 0
    do ii = 1, 500
      sm = 0.0_dp
      do ip = 1, dimen - 1
        do iq = ip + 1, dimen
          sm = sm + abs(a(ip,iq))
        end do
      end do

      if (sm == 0.0_dp) return

      if (ii < 4) then
        tresh = 0.2_dp * sm / real(dimen, dp)**2
      else
        tresh = 0.0_dp
      end if

      do ip = 1, dimen - 1
        do iq = ip + 1, dimen
          g = 100.0_dp * abs(a(ip,iq))
          if ((ii > 4) .and. (abs(eigval(ip)) + g == abs(eigval(ip))) &
              .and. (abs(eigval(iq)) + g == abs(eigval(iq)))) then
            a(ip,iq) = 0.0_dp
          else if (abs(a(ip,iq)) > tresh) then
            h = eigval(iq) - eigval(ip)
            if (abs(h) + g == abs(h)) then
              t = a(ip,iq) / h
            else
              theta = 0.5_dp * h / a(ip,iq)
              t = 1.0_dp / (abs(theta) + sqrt(1.0_dp + theta**2))
              if (theta < 0.0_dp) t = -t
            end if

            c = 1.0_dp / sqrt(1.0_dp + t**2)
            s = t * c
            tau = s / (1.0_dp + c)
            h = t * a(ip,iq)
            z(ip) = z(ip) - h
            z(iq) = z(iq) + h
            eigval(ip) = eigval(ip) - h
            eigval(iq) = eigval(iq) + h
            a(ip,iq) = 0.0_dp

            do jj = 1, ip - 1
              g = a(jj,ip)
              h = a(jj,iq)
              a(jj,ip) = g - s*(h + g*tau)
              a(jj,iq) = h + s*(g - h*tau)
            end do
            do jj = ip + 1, iq - 1
              g = a(ip,jj)
              h = a(jj,iq)
              a(ip,jj) = g - s*(h + g*tau)
              a(jj,iq) = h + s*(g - h*tau)
            end do
            do jj = iq + 1, dimen
              g = a(ip,jj)
              h = a(iq,jj)
              a(ip,jj) = g - s*(h + g*tau)
              a(iq,jj) = h + s*(g - h*tau)
            end do
            do jj = 1, dimen
              g = eigvec(jj,ip)
              h = eigvec(jj,iq)
              eigvec(jj,ip) = g - s*(h + g*tau)
              eigvec(jj,iq) = h + s*(g - h*tau)
            end do
            nrot = nrot + 1
          end if
        end do
      end do

      do ip = 1, dimen
        b(ip) = b(ip) + z(ip)
        eigval(ip) = b(ip)
        z(ip) = 0.0_dp
      end do
    end do
  end subroutine decomp_matrix

  ! Verbatim port of subroutine multi_matrix (conform.f), dimen=2.
  subroutine multi_matrix(a, b, m)
    real(dp), intent(in) :: a(dimen,dimen), b(dimen,dimen)
    real(dp), intent(out) :: m(dimen,dimen)
    integer :: ii, jj, kk
    do ii = 1, dimen
      do jj = 1, dimen
        m(ii,jj) = 0.0_dp
        do kk = 1, dimen
          m(ii,jj) = m(ii,jj) + a(ii,kk) * b(kk,jj)
        end do
      end do
    end do
  end subroutine multi_matrix

  ! Pointwise port of subroutine T_to_Psi (conform.f).
  subroutine T_to_Psi_point(txx_p, txy_p, tyy_p, etap_p, psixx_p, psixy_p, psiyy_p)
    real(dp), intent(in) :: txx_p, txy_p, tyy_p, etap_p
    real(dp), intent(out) :: psixx_p, psixy_p, psiyy_p
    real(dp) :: a(dimen,dimen), eigvec(dimen,dimen), eigvec_t(dimen,dimen)
    real(dp) :: psic(dimen,dimen), psii(dimen,dimen), llambda(dimen,dimen)
    real(dp) :: eigval(dimen)
    integer :: k

    a(1,1) = txx_p / etap_p + 1.0_dp
    a(1,2) = txy_p / etap_p
    a(2,1) = a(1,2)
    a(2,2) = tyy_p / etap_p + 1.0_dp

    call decomp_matrix(a, eigval, eigvec)
    eigvec_t = transpose(eigvec)

    llambda = 0.0_dp
    do k = 1, dimen
      llambda(k,k) = log(eigval(k))
    end do

    call multi_matrix(eigvec, llambda, psii)
    call multi_matrix(psii, eigvec_t, psic)

    psixx_p = psic(1,1)
    psixy_p = psic(1,2)
    psiyy_p = psic(2,2)
  end subroutine T_to_Psi_point

  ! Pointwise port of subroutine Psi_to_T (conform.f).
  subroutine Psi_to_T_point(psixx_p, psixy_p, psiyy_p, etap_p, txx_p, txy_p, tyy_p)
    real(dp), intent(in) :: psixx_p, psixy_p, psiyy_p, etap_p
    real(dp), intent(out) :: txx_p, txy_p, tyy_p
    real(dp) :: psic(dimen,dimen), ai(dimen,dimen), a(dimen,dimen)
    real(dp) :: elambda(dimen,dimen), eigvec(dimen,dimen), eigvec_t(dimen,dimen)
    real(dp) :: eigval(dimen)
    integer :: k

    psic(1,1) = psixx_p
    psic(1,2) = psixy_p
    psic(2,1) = psixy_p
    psic(2,2) = psiyy_p

    call decomp_matrix(psic, eigval, eigvec)
    eigvec_t = transpose(eigvec)

    elambda = 0.0_dp
    do k = 1, dimen
      elambda(k,k) = exp(eigval(k))
    end do

    call multi_matrix(eigvec, elambda, ai)
    call multi_matrix(ai, eigvec_t, a)

    txx_p = etap_p * (a(1,1) - 1.0_dp)
    txy_p = etap_p * a(1,2)
    tyy_p = etap_p * (a(2,2) - 1.0_dp)
  end subroutine Psi_to_T_point

  ! Pointwise port of subroutine kernel_conformation (conform.f).
  subroutine kernel_conformation_point(txx_p, txy_p, tyy_p, dudx_p, dudy_p, dvdx_p, dvdy_p, &
                                        etapinv_p, alphaG_p, &
                                        omegalog_p, bxx_p, bxy_p, byy_p, mxx_p, mxy_p, myy_p)
    real(dp), intent(in) :: txx_p, txy_p, tyy_p, dudx_p, dudy_p, dvdx_p, dvdy_p
    real(dp), intent(in) :: etapinv_p, alphaG_p
    real(dp), intent(out) :: omegalog_p, bxx_p, bxy_p, byy_p, mxx_p, mxy_p, myy_p

    real(dp) :: a(dimen,dimen), eigvec(dimen,dimen), eigvec_t(dimen,dimen)
    real(dp) :: id(dimen,dimen), jac(dimen,dimen), gradu(dimen,dimen)
    real(dp) :: mi(dimen,dimen), m(dimen,dimen), eigval(dimen,dimen)
    real(dp) :: m_lambda(dimen,dimen), bmat(dimen,dimen)
    real(dp) :: eigval_vec(dimen), div, omegal
    integer :: k

    a(1,1) = txx_p * etapinv_p + 1.0_dp
    a(1,2) = txy_p * etapinv_p
    a(2,1) = a(1,2)
    a(2,2) = tyy_p * etapinv_p + 1.0_dp

    call decomp_matrix(a, eigval_vec, eigvec)
    eigvec_t = transpose(eigvec)

    eigval = 0.0_dp
    eigval(1,1) = eigval_vec(1)
    eigval(2,2) = eigval_vec(2)

    jac = 0.0_dp
    jac(1,1) = 1.0_dp / eigval_vec(1)
    jac(2,2) = 1.0_dp / eigval_vec(2)

    gradu(1,1) = dudx_p
    gradu(1,2) = dudy_p
    gradu(2,1) = dvdx_p
    gradu(2,2) = dvdy_p
    call multi_matrix(eigvec_t, gradu, mi)
    call multi_matrix(mi, eigvec, m)

    div = eigval_vec(2) - eigval_vec(1)
    if (abs(div) < 1.0e-14_dp) div = 1.0e-14_dp
    omegal = (m(1,2) * eigval_vec(2) + m(2,1) * eigval_vec(1)) / div

    m(1,2) = 0.0_dp
    m(2,1) = 0.0_dp
    call multi_matrix(eigvec, m, mi)
    call multi_matrix(mi, eigval, m)
    call multi_matrix(m, jac, mi)
    call multi_matrix(mi, eigvec_t, bmat)
    bxx_p = bmat(1,1)
    bxy_p = bmat(1,2)
    byy_p = bmat(2,2)

    id = 0.0_dp
    do k = 1, dimen
      id(k,k) = 1.0_dp
    end do
    call multi_matrix(id - eigval, id + alphaG_p * (eigval - id), m_lambda)

    call multi_matrix(eigvec, m_lambda, mi)
    call multi_matrix(mi, jac, m)
    call multi_matrix(m, eigvec_t, mi)
    mxx_p = mi(1,1)
    mxy_p = mi(1,2)
    myy_p = mi(2,2)

    m(1,1) = 0.0_dp
    m(1,2) = omegal
    m(2,1) = -omegal
    m(2,2) = 0.0_dp
    call multi_matrix(eigvec, m, mi)
    call multi_matrix(mi, eigvec_t, m)
    omegalog_p = m(1,2)
  end subroutine kernel_conformation_point

end program case12_log_euler_step_reference
