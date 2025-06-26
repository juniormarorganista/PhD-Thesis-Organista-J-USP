ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             Poisson solver subroutine                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Parallel multigrid Poisson solver for the v-poisson. It makes v cicles as shown in the figure below until the defect is lower then the prescribed error (erromax).
!> \image html FAS.jpg 'Multigrid method FAS - V cicle'
!> where \f$ N_{1},~N_{2},~N_{3} \f$ correspond to the iteractions number for each mesh, \f$ R \f$ represents the restriction process to pass the results from the finest to the coarsest mesh and \f$ P \f$ represents the prolongation process that is responsible to pass the results from the coarsest to the finest mesh.
      subroutine poi_psi4a

      implicit none
      include 'par.nn'
      include 'comm.var'
      logical continua
      integer itera, nptsx, nptsy, i_ini, i_fim, mesh, vcycles
      real*8 dx2, dy2, err, erromax, fc(msh,imax,jmax), 
     & sc(msh,imax,jmax), fct(msh,imax,jmax)

      erromax = 1.d-06

      ! Initialization of variables
      nptsx = imax
      nptsy = jmax
      dx2   = dxx
      dy2   = dyy
      itera = 2
      i_ini = 2
      i_fim = imax - 1

      call function_source(fc, sc, fct)
      continua = .true.
      do vcycles = 1, 81
        do mesh = 1, msh - 1  ! from finest to coarsest mesh <<<<<
          call sor(fc, sc, nptsx, nptsy, dy2, itera, mesh)
          call defect(fc,sc,fct,nptsx,nptsy,dx2,dy2,-1.d0,mesh)
          if (mesh.eq.1) then ! stop condition
            call verifica_poi(fct, err, erromax, continua)
            if (.not.continua) goto 22
          end if              ! stop condition
          call restrict(fc,sc,fct,nptsx,nptsy,dx2,dy2,mesh)
          call defect(fc,sc,fct,nptsx,nptsy,dx2,dy2,1.d0,mesh+1)
        end do                ! from finest to coarsest mesh >>>>>
        ! coarsest mesh
        call sor(fc, sc, nptsx, nptsy, dy2, 100*itera, msh)
        do mesh = msh, 2, - 1 ! from coarsest to finest mesh <<<<<
          call correction(fc,fct,nptsx,nptsy,mesh)
          call intpol(fc,nptsx,nptsy,dx2,dy2, mesh)
          call sor(fc, sc, nptsx, nptsy, dy2, 1, mesh-1)
        end do                ! from coarsest to finest mesh >>>>>
      end do
  22  psi = fc(1,:,:)

      if (vcycles.gt.80) then
         write(*,*) 'Vcycles = ', vcycles
         stop
      end if 

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Verifies the max error (err) of the process and compare with a tolerance error (errmax)
      subroutine verifica_poi(fct, err, erromax, continua)

      implicit none
      include 'par.nn'
      logical continua
      integer i, j
      real*8 err, erromax, fct(msh,imax,jmax)

      err = 0.d0

      do j = 2, jmax - 1
        do i = 2, imax - 1
          err = max( err, abs(fct(1,i,j)) )
        end do
      end do

      if (err .lt. erromax) continua = .false.

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine function_source(fc, sc, fct)

      ! calculates the rhs weighted of the finest grid Poisson equation
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 fc(msh,imax,jmax), sc(msh,imax,jmax), fct(msh,imax,jmax)

      ! Source term of the Poisson Equation
      fc  = 0.d0
      sc  = 0.d0
      fct = 0.d0

      do j = 2, jmax - 1
        do i = 2, imax - 1
          ! weight in the middle
          sc(1,i,j)=( wz(i,j-1)+wz(i,j+1)+wz(i+1,j)+wz(i-1,j)+
     &           8.d0*wz(i,j) )
        end do
      end do

      ! function of the Poisson equation
      fc(1,:,:) = psi(:,:)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Calculatez the SOR method to solve the poisson equation
      subroutine sor(fc, sc, nptsx, nptsy, dy2, itera, mesh)

      implicit none
      include 'par.nn'
      integer i, j, it, itera, nptsx, nptsy, mesh
      real*8 a(jmax), b(jmax), c(jmax),
     &       dy2, rf, dyc, fc(msh,imax,jmax), sc(msh,imax,jmax), 
     &       rhs(jmax)

      dyc  = dsqrt(dy2)
      rf   = 1.2d0
      if (itera.eq.1) rf = 1.d0

      ! LHS calculation
      a(1) = 0.d0
      b(1) = 1.d0
      c(1) = 0.d0
      do j = 2, nptsy - 1
        a(j) =   10.d0 - 2.d0*dyypdxx
        b(j) = - 20.d0 * (1.d0+dyypdxx)
        c(j) =   10.d0 - 2.d0*dyypdxx
      end do
      a(nptsy) = 0.d0
      b(nptsy) = 1.d0
      c(nptsy) = 0.d0

      do it = 1, itera

        do i = 2, nptsx - 1
          rhs(1) = fc(mesh,i,1)
          do j = 2, nptsy - 1
            rhs(j)=sc(mesh,i,j)*dy2-
     &    ( (1.d0+dyypdxx)*(fc(mesh,i+1,j+1) + fc(mesh,i-1,j+1)+
     &                      fc(mesh,i+1,j-1) + fc(mesh,i-1,j-1))+
     &    (-2.d0+10.d0*dyypdxx)*(fc(mesh,i+1,j)   + fc(mesh,i-1,j)) )
          end do
          rhs(nptsy) = fc(mesh,i,nptsy)
          call tridv(rhs,nptsy,a,b,c) ! solve the matrix
          do j = 2, nptsy - 1
            fc(mesh,i,j) = fc(mesh,i,j) + rf * (rhs(j) - fc(mesh,i,j))
          end do
        end do
      end do
    
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defect(fc,sc,fct,nptsx,nptsy,dx2,dy2,si,mesh)

      implicit none
      include 'par.nn'
      integer i, j, nptsx, nptsy, mesh
      real*8 dx2, dy2, si, dyc, b, c, fc(msh,imax,jmax), 
     &       sc(msh,imax,jmax), fct(msh,imax,jmax)

      dyc    =   dsqrt(dy2)

      do j = 2, nptsy - 1
        do i = 2, nptsx - 1
          b=(  fc(mesh,i+1,j+1) + fc(mesh,i+1,j-1)
     &       + fc(mesh,i-1,j+1) + fc(mesh,i-1,j-1)
     & +10.d0*(fc(mesh,i+1,j)   + fc(mesh,i-1,j) )
     & - 2.d0*(fc(mesh,i,j+1)   + fc(mesh,i,j-1) )
     & -20.d0* fc(mesh,i,j) )
          c=(  fc(mesh,i+1,j+1) + fc(mesh,i+1,j-1)
     &       + fc(mesh,i-1,j+1) + fc(mesh,i-1,j-1)
     & +10.d0*(fc(mesh,i,j+1)  +  fc(mesh,i,j-1) )
     & - 2.d0*(fc(mesh,i+1,j)  +  fc(mesh,i-1,j) )
     & -20.d0* fc(mesh,i,j) )
          fct(1,i,j) = sc(mesh,i,j) + si * ( b / dx2 + c / dy2 )
        end do 
      end do

      if (si.lt.0) return

      do j = 2, nptsy - 1
        do i = 2, nptsx - 1
          sc(mesh,i,j) = fct(1,i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Makes the restriction process of the multigrid method, where the defect is transmited from the finest mesh \f$ h \f$ to another less refined mesh \f$2h\f$ by means of a weighting (FW - Full Weight). The velocity is passed from one mesh to another by direct injection (SI - Straight Injection).
!!\f{eqnarray*}{v_{h} \Rightarrow v_{2h}\ (SI)\f}
!!\f{eqnarray*}{d_{h} \Rightarrow d_{2h}\ (FW)\f}
      subroutine restrict(fc,sc,fct,nptsx,nptsy,dx2,dy2,mesh)

      implicit none
      include 'par.nn'
      integer i, j, nptsx, nptsy, fi, fj, mesh
      real*8 dx2, dy2, fc(msh,imax,jmax), sc(msh,imax,jmax),
     &       fct(msh,imax,jmax)

      ! new values for dx2, dy2, nptsx e nptsy
      nptsx = (nptsx + 1)/ 2
      nptsy = (nptsy + 1)/ 2
      dx2   = 4.d0 * dx2
      dy2   = 4.d0 * dy2

      ! restriction in the middle
      fj = 3
      do j = 2, nptsy - 1
        fi = 3
        do i = 2, nptsx - 1
          fc(mesh+1,i,j)  = fc(mesh,fi,fj)
          fct(mesh+1,i,j) = fc(mesh+1,i,j)
          sc(mesh+1,i,j)  = 6.25d-2 * ( 2.d0 *   ( fct(1,fi-1,fj)  +
     &                      fct(1,fi+1,fj)  +      fct(1,fi,fj-1)  +
     &                      fct(1,fi,fj+1) )+ 4.d0*fct(1,fi,fj)    +
     &                      fct(1,fi+1,fj+1)+      fct(1,fi+1,fj-1)+
     &                      fct(1,fi-1,fj+1)+      fct(1,fi-1,fj-1) )
          fi = fi + 2
        end do
        fj = fj + 2
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Makes the bilinear interpolation responsible to return to a finest level mesh, where
!!\f{eqnarray*}{corr_{8h} \Rightarrow corr_{4h}\f}
!!\f{eqnarray*}{v_{4h} \Leftarrow v_{4h} + corr_{4h}\f}
      subroutine intpol(fc,nptsx,nptsy,dx2,dy2,mesh)

      implicit none
      include 'par.nn'
      integer i, j, nptsx, nptsy, fi, fj, mesh
      real*8 dx2, dy2, fc(msh,imax,jmax), temp(imax,jmax)

      ! new values for dx2 e dy2
      dx2 = 0.25d0 * dx2
      dy2 = 0.25d0 * dy2

      ! copying fcc to temp
      fj = 1
      do j = 1, nptsy
        fi = 1
        do i = 1, nptsx
          temp(fi,fj) = fc(mesh,i,j)
          fi = fi + 2
        end do
        fj = fj + 2
      end do

      ! new values for nptsx and nptsy
      nptsx = 2 * nptsx - 1
      nptsy = 2 * nptsy - 1

      ! Interpolation in x
      do j = 1, nptsy, 2
        do i = 2, nptsx - 1, 2
          temp(i,j) = 0.5d0 * ( temp(i-1,j) + temp(i+1,j) )
        end do
      end do

      ! Interpolation in y
      do j = 2, nptsy - 1, 2
        do i = 1, nptsx
          temp(i,j) = 0.5d0 * ( temp(i,j+1) + temp(i,j-1) )
        end do
      end do

      ! adition to new function
      do j = 2, nptsy - 1
        do i = 2, nptsx - 1
          fc(mesh-1,i,j) = fc(mesh-1,i,j) + temp(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Calculate the correction for the coarsest mesh by 
!!\f{eqnarray*}{corr_{8h} = v_{8h}^{n} - v_{8h}^{'}\f}
!!
!!where \f$v_{8h}^{'}\f$ represents the approximation generated by the restriction operation and \f$v_{8h}^{n}\f$ represents the approximation that was calculated recently.
      subroutine correction(fc, fct, nptsx, nptsy, mesh)

      implicit none
      include 'par.nn'
      integer i, j, nptsx, nptsy, mesh
      real*8 fc(msh,imax,jmax), fct(msh,imax,jmax)

      do j = 2, nptsy - 1
        do i = 2, nptsx - 1
          fc(mesh,i,j) = fc(mesh,i,j) - fct(mesh,i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridv(r,lmax,a,b,c) ! retirado do Numerical Recipes
      
      ! resolve um sistema tridiagonal
      implicit none
      include 'par.nn'
      integer j,lmax
      real*8 a(jmax),b(jmax),c(jmax),r(jmax),u(jmax)
      real*8 bet,gam(jmax)
      
      bet  = b(1)
      u(1) = r(1)/bet
      do j = 2, lmax
        gam(j) = c(j-1)/bet
        bet    = b(j)-a(j)*gam(j)
        u(j)   = (r(j)-a(j)*u(j-1))/bet
      end do
      do j = lmax - 1, 1, -1
        u(j) = u(j)-gam(j+1)*u(j+1)
      end do
      do j = 1, lmax
        r(j) = u(j)
      end do

      return
      end
