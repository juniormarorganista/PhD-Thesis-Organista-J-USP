program case12_poisson_multigrid_reference
  !! Standalone Fortran reference for the multigrid V-cycle Poisson solver
  !! (Etapa 17).
  !!
  !! This utility does not modify or depend on the original thesis Fortran
  !! files. It is a verbatim port of:
  !!   src/fortran_original/4Order_LPOG_Euler/poisson.f
  !!     poi_psi4a, verifica_poi, function_source, sor, defect, restrict,
  !!     intpol, correction, tridv
  !!
  !! laplacian(psi) = wz is solved (this package's sign convention, see
  !! mms_viscoelastic.numerics.convergence) using the 4th-order compact
  !! (Mehrstellen-type) 9-point discretization and the FAS (Full
  !! Approximation Scheme) multigrid V-cycle implemented by poi_psi4a, NOT
  !! the pentadiagonal band5_poi/banbk5_poi routines in solvers.f, which are
  !! never called anywhere in the original Fortran source (dead code).
  !!
  !! The input CSV provides wz and the boundary/initial-guess field psi on a
  !! uniform grid covering [0,1] x [0,1] (nx = ny = 33 by default, matching
  !! par.nn's imax = jmax = 33). Only the interior of psi is solved for; the
  !! boundary is copied unchanged (psi's boundary values are never
  !! overwritten by sor, matching the Fortran: the tridiagonal system pins
  !! row 1 and row nptsy to their current value at every sweep).
  !!
  !! Usage:
  !!   gfortran -O2 -std=f2008 case12_poisson_multigrid_reference.f90 -o case12_poisson_multigrid_reference.exe
  !!   .\case12_poisson_multigrid_reference.exe input.csv output.csv 33 33 4
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  integer :: nx, ny, msh, i, j
  real(dp) :: dx, dy
  character(len=512) :: infile, outfile, arg

  real(dp), allocatable :: xg(:,:), yg(:,:), wz(:,:), psi(:,:)
  real(dp), allocatable :: fc(:,:,:), sc(:,:,:), fct(:,:,:)

  nx = 33
  ny = 33
  msh = 4
  infile = "case12_poisson_multigrid_input_33x33.csv"
  outfile = "case12_poisson_multigrid_output_33x33.csv"

  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) infile = trim(arg)
  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) outfile = trim(arg)
  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg, *) nx
  call get_command_argument(4, arg)
  if (len_trim(arg) > 0) read(arg, *) ny
  call get_command_argument(5, arg)
  if (len_trim(arg) > 0) read(arg, *) msh

  if (nx <= 2 .or. ny <= 2) error stop "nx and ny must be greater than 2"

  dx = 1.0_dp / real(nx - 1, dp)
  dy = 1.0_dp / real(ny - 1, dp)

  allocate(xg(nx,ny), yg(nx,ny), wz(nx,ny), psi(nx,ny))

  open(unit=10, file=trim(infile), status="old", action="read")
  read(10, *)  ! skip header
  do j = 1, ny
    do i = 1, nx
      read(10, *) xg(i,j), yg(i,j), wz(i,j), psi(i,j)
    end do
  end do
  close(10)

  allocate(fc(msh,nx,ny), sc(msh,nx,ny), fct(msh,nx,ny))

  call poi_psi4a(nx, ny, msh, dx*dx, dy*dy, wz, psi, fc, sc, fct)

  open(unit=20, file=trim(outfile), status="replace", action="write")
  write(20, '(A)') "x,y,psi"
  do j = 1, ny
    do i = 1, nx
      write(20, '(ES26.17E3,A,ES26.17E3,A,ES26.17E3)') xg(i,j), ',', yg(i,j), ',', psi(i,j)
    end do
  end do
  close(20)

contains

  subroutine poi_psi4a(imax, jmax, msh, dxx, dyy, wz, psi, fc, sc, fct)
    integer, intent(in) :: imax, jmax, msh
    real(dp), intent(in) :: dxx, dyy
    real(dp), intent(in) :: wz(imax,jmax)
    real(dp), intent(inout) :: psi(imax,jmax)
    real(dp), intent(inout) :: fc(msh,imax,jmax), sc(msh,imax,jmax), fct(msh,imax,jmax)

    logical :: continua
    integer :: itera, nptsx, nptsy, mesh, vcycles
    real(dp) :: dx2, dy2, err, erromax, dyypdxx

    erromax = 1.0e-06_dp
    nptsx = imax
    nptsy = jmax
    dx2 = dxx
    dy2 = dyy
    itera = 2

    call function_source(imax, jmax, msh, wz, psi, fc, sc, fct)
    continua = .true.

    do vcycles = 1, 81
      do mesh = 1, msh - 1
        dyypdxx = dy2 / dx2
        call sor(imax, jmax, msh, fc, sc, nptsx, nptsy, dy2, dyypdxx, itera, mesh)
        call defect(imax, jmax, msh, fc, sc, fct, nptsx, nptsy, dx2, dy2, -1.0_dp, mesh)
        if (mesh == 1) then
          call verifica_poi(imax, jmax, msh, fct, err, erromax, continua)
          if (.not. continua) goto 22
        end if
        call restrict(imax, jmax, msh, fc, sc, fct, nptsx, nptsy, dx2, dy2, mesh)
        dyypdxx = dy2 / dx2
        call defect(imax, jmax, msh, fc, sc, fct, nptsx, nptsy, dx2, dy2, 1.0_dp, mesh+1)
      end do

      dyypdxx = dy2 / dx2
      call sor(imax, jmax, msh, fc, sc, nptsx, nptsy, dy2, dyypdxx, 100*itera, msh)

      do mesh = msh, 2, -1
        call correction(imax, jmax, msh, fc, fct, nptsx, nptsy, mesh)
        call intpol(imax, jmax, msh, fc, nptsx, nptsy, dx2, dy2, mesh)
        dyypdxx = dy2 / dx2
        call sor(imax, jmax, msh, fc, sc, nptsx, nptsy, dy2, dyypdxx, 1, mesh-1)
      end do
    end do

22  continue
    psi(:,:) = fc(1,:,:)

    if (vcycles > 80) then
      write(*,*) 'Vcycles = ', vcycles
      error stop
    end if
  end subroutine poi_psi4a

  subroutine verifica_poi(imax, jmax, msh, fct, err, erromax, continua)
    integer, intent(in) :: imax, jmax, msh
    real(dp), intent(in) :: fct(msh,imax,jmax), erromax
    real(dp), intent(out) :: err
    logical, intent(inout) :: continua
    integer :: i, j

    err = 0.0_dp
    do j = 2, jmax - 1
      do i = 2, imax - 1
        err = max(err, abs(fct(1,i,j)))
      end do
    end do
    if (err < erromax) continua = .false.
  end subroutine verifica_poi

  subroutine function_source(imax, jmax, msh, wz, psi, fc, sc, fct)
    integer, intent(in) :: imax, jmax, msh
    real(dp), intent(in) :: wz(imax,jmax), psi(imax,jmax)
    real(dp), intent(out) :: fc(msh,imax,jmax), sc(msh,imax,jmax), fct(msh,imax,jmax)
    integer :: i, j

    fc = 0.0_dp
    sc = 0.0_dp
    fct = 0.0_dp

    do j = 2, jmax - 1
      do i = 2, imax - 1
        sc(1,i,j) = wz(i,j-1) + wz(i,j+1) + wz(i+1,j) + wz(i-1,j) + 8.0_dp*wz(i,j)
      end do
    end do

    fc(1,:,:) = psi(:,:)
  end subroutine function_source

  subroutine sor(imax, jmax, msh, fc, sc, nptsx, nptsy, dy2, dyypdxx, itera, mesh)
    integer, intent(in) :: imax, jmax, msh, nptsx, nptsy, itera, mesh
    real(dp), intent(in) :: dy2, dyypdxx
    real(dp), intent(inout) :: fc(msh,imax,jmax)
    real(dp), intent(in) :: sc(msh,imax,jmax)
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
        rhs(1) = fc(mesh,i,1)
        do j = 2, nptsy - 1
          rhs(j) = sc(mesh,i,j)*dy2 - ( (1.0_dp+dyypdxx)*(fc(mesh,i+1,j+1) + fc(mesh,i-1,j+1) &
                    + fc(mesh,i+1,j-1) + fc(mesh,i-1,j-1)) &
                    + (-2.0_dp+10.0_dp*dyypdxx)*(fc(mesh,i+1,j) + fc(mesh,i-1,j)) )
        end do
        rhs(nptsy) = fc(mesh,i,nptsy)
        call tridv(rhs, nptsy, a, b, c, jmax)
        do j = 2, nptsy - 1
          fc(mesh,i,j) = fc(mesh,i,j) + rf * (rhs(j) - fc(mesh,i,j))
        end do
      end do
    end do
  end subroutine sor

  subroutine defect(imax, jmax, msh, fc, sc, fct, nptsx, nptsy, dx2, dy2, si, mesh)
    integer, intent(in) :: imax, jmax, msh, nptsx, nptsy, mesh
    real(dp), intent(in) :: dx2, dy2, si
    real(dp), intent(in) :: fc(msh,imax,jmax)
    real(dp), intent(inout) :: sc(msh,imax,jmax)
    real(dp), intent(inout) :: fct(msh,imax,jmax)
    integer :: i, j
    real(dp) :: b, c

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        b = ( fc(mesh,i+1,j+1) + fc(mesh,i+1,j-1) + fc(mesh,i-1,j+1) + fc(mesh,i-1,j-1) &
              + 10.0_dp*(fc(mesh,i+1,j) + fc(mesh,i-1,j)) &
              - 2.0_dp*(fc(mesh,i,j+1) + fc(mesh,i,j-1)) &
              - 20.0_dp*fc(mesh,i,j) )
        c = ( fc(mesh,i+1,j+1) + fc(mesh,i+1,j-1) + fc(mesh,i-1,j+1) + fc(mesh,i-1,j-1) &
              + 10.0_dp*(fc(mesh,i,j+1) + fc(mesh,i,j-1)) &
              - 2.0_dp*(fc(mesh,i+1,j) + fc(mesh,i-1,j)) &
              - 20.0_dp*fc(mesh,i,j) )
        fct(1,i,j) = sc(mesh,i,j) + si * (b/dx2 + c/dy2)
      end do
    end do

    if (si < 0) return

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        sc(mesh,i,j) = fct(1,i,j)
      end do
    end do
  end subroutine defect

  subroutine restrict(imax, jmax, msh, fc, sc, fct, nptsx, nptsy, dx2, dy2, mesh)
    integer, intent(in) :: imax, jmax, msh, mesh
    integer, intent(inout) :: nptsx, nptsy
    real(dp), intent(inout) :: dx2, dy2
    real(dp), intent(inout) :: fc(msh,imax,jmax)
    real(dp), intent(inout) :: sc(msh,imax,jmax)
    real(dp), intent(inout) :: fct(msh,imax,jmax)
    integer :: i, j, fi, fj

    nptsx = (nptsx + 1) / 2
    nptsy = (nptsy + 1) / 2
    dx2 = 4.0_dp * dx2
    dy2 = 4.0_dp * dy2

    fj = 3
    do j = 2, nptsy - 1
      fi = 3
      do i = 2, nptsx - 1
        fc(mesh+1,i,j) = fc(mesh,fi,fj)
        fct(mesh+1,i,j) = fc(mesh+1,i,j)
        sc(mesh+1,i,j) = 6.25e-2_dp * ( 2.0_dp * ( fct(1,fi-1,fj) + fct(1,fi+1,fj) &
                          + fct(1,fi,fj-1) + fct(1,fi,fj+1) ) + 4.0_dp*fct(1,fi,fj) &
                          + fct(1,fi+1,fj+1) + fct(1,fi+1,fj-1) &
                          + fct(1,fi-1,fj+1) + fct(1,fi-1,fj-1) )
        fi = fi + 2
      end do
      fj = fj + 2
    end do
  end subroutine restrict

  subroutine intpol(imax, jmax, msh, fc, nptsx, nptsy, dx2, dy2, mesh)
    integer, intent(in) :: imax, jmax, msh, mesh
    integer, intent(inout) :: nptsx, nptsy
    real(dp), intent(inout) :: dx2, dy2
    real(dp), intent(inout) :: fc(msh,imax,jmax)
    integer :: i, j, fi, fj
    real(dp) :: temp(imax,jmax)

    dx2 = 0.25_dp * dx2
    dy2 = 0.25_dp * dy2

    temp = 0.0_dp
    fj = 1
    do j = 1, nptsy
      fi = 1
      do i = 1, nptsx
        temp(fi,fj) = fc(mesh,i,j)
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
        fc(mesh-1,i,j) = fc(mesh-1,i,j) + temp(i,j)
      end do
    end do
  end subroutine intpol

  subroutine correction(imax, jmax, msh, fc, fct, nptsx, nptsy, mesh)
    integer, intent(in) :: imax, jmax, msh, nptsx, nptsy, mesh
    real(dp), intent(inout) :: fc(msh,imax,jmax)
    real(dp), intent(in) :: fct(msh,imax,jmax)
    integer :: i, j

    do j = 2, nptsy - 1
      do i = 2, nptsx - 1
        fc(mesh,i,j) = fc(mesh,i,j) - fct(mesh,i,j)
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

end program case12_poisson_multigrid_reference
