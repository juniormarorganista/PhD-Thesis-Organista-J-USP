cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine T_to_Psi

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j, k
      real*8   A(dimen,dimen),  eigvec(dimen,dimen),
     &  eigvec_t(dimen,dimen),    Psic(dimen,dimen),
     &      PsiI(dimen,dimen), Llambda(dimen,dimen), 
     &    eigval(dimen)

      Llambda = 0.d0

      do j = 1, jmax
        do i = 1, imax

          ! Conformation matrix
          A(1,1) = Txx(i,j) * etapinv + 1.d0
          A(1,2) = Txy(i,j) * etapinv
          A(2,1) =   A(1,2)
          A(2,2) = Tyy(i,j) * etapinv + 1.d0

          call decomp_matrix(A,eigval,eigvec)
!         eigvec(:,1) = - eigvec(:,1)
          eigvec_t = transpose(eigvec)

          do k = 1, dimen
            Llambda(k,k) = dlog(eigval(k))
!           Llambda(k,k) = dsqrt(eigval(k))
!           Llambda(k,k) = eigval(k)
          end do

          call multi_matrix(eigvec,Llambda,PsiI)
          call multi_matrix(PsiI,eigvec_t,Psic)

          Psixx(i,j) = Psic(1,1)
          Psixy(i,j) = Psic(1,2)
          Psiyy(i,j) = Psic(2,2)

        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cc_Psi

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j, k
      real*8   A(dimen,dimen),  eigvec(dimen,dimen),
     &  eigvec_t(dimen,dimen),    Psic(dimen,dimen),
     &      PsiI(dimen,dimen), Llambda(dimen,dimen), 
     &    eigval(dimen)

      Llambda = 0.d0

      do j = 1, jmax
        do i = 1, imax, imax - 1

          ! Conformation matrix
          A(1,1) = Txx(i,j) * etapinv + 1.d0
          A(1,2) = Txy(i,j) * etapinv
          A(2,1) = Txy(i,j) * etapinv
          A(2,2) = Tyy(i,j) * etapinv + 1.d0

          call decomp_matrix(A,eigval,eigvec)
!          eigvec(:,1) = - eigvec(:,1)
          eigvec_t = transpose(eigvec)

          do k = 1, dimen
            Llambda(k,k) = dlog(eigval(k))
!           Llambda(k,k) = dsqrt(eigval(k))
!           Llambda(k,k) = eigval(k)
          end do

          call multi_matrix(eigvec,Llambda,PsiI)
          call multi_matrix(PsiI,eigvec_t,Psic)

          Psixx(i,j) = Psic(1,1)
          Psixy(i,j) = Psic(1,2)
          Psiyy(i,j) = Psic(2,2)

        end do
      end do

      do j = 1, jmax, jmax - 1
        do i = 1, imax

          ! Conformation matrix
          A(1,1) = Txx(i,j) * etapinv + 1.d0
          A(1,2) = Txy(i,j) * etapinv
          A(2,1) = Txy(i,j) * etapinv
          A(2,2) = Tyy(i,j) * etapinv + 1.d0

          call decomp_matrix(A,eigval,eigvec)
!          eigvec(:,1) = - eigvec(:,1)
          eigvec_t = transpose(eigvec)

          do k = 1, dimen
            Llambda(k,k) = dlog(eigval(k))
!           Llambda(k,k) = dsqrt(eigval(k))
!           Llambda(k,k) = eigval(k)
          end do

          call multi_matrix(eigvec,Llambda,PsiI)
          call multi_matrix(PsiI,eigvec_t,Psic)

          Psixx(i,j) = Psic(1,1)
          Psixy(i,j) = Psic(1,2)
          Psiyy(i,j) = Psic(2,2)

        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine Psi_to_T

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j, k
      real*8 Psic(dimen,dimen),       AI(dimen,dimen),  
     &          A(dimen,dimen),  elambda(dimen,dimen),
     &     eigvec(dimen,dimen), eigvec_t(dimen,dimen),
     &     eigval(dimen)

      elambda = 0.d0

      do j = 1, jmax
        do i = 1, imax

          Psic(1,1) = Psixx(i,j)
          Psic(1,2) = Psixy(i,j)
          Psic(2,1) = Psic(1,2)
          Psic(2,2) = Psiyy(i,j)
         
          call decomp_matrix(Psic,eigval,eigvec)
!         eigvec(:,1) = - eigvec(:,1)
          eigvec_t = transpose(eigvec)

          do k = 1, dimen
            elambda(k,k) = dexp(eigval(k))
!           elambda(k,k) = (eigval(k))*(eigval(k))
!           elambda(k,k) = eigval(k)
          end do

          call multi_matrix(eigvec,elambda,AI)
          call multi_matrix(AI,eigvec_t,A)

          ! Componentes do Tensor T
          Txx(i,j) = etap * (A(1,1) - 1.d0)
          Txy(i,j) = etap *  A(1,2)
          Tyy(i,j) = etap * (A(2,2) - 1.d0)

        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tf_conformation

      implicit none
      include 'par.nn'
      include 'comm.var'

      tfpsixx = tfTxx * etapinv
      tfpsixy = tfTxy * etapinv
      tfpsiyy = tfTyy * etapinv

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine kernel_conformation(Omegalog,Bxx,Bxy,Byy,Mxx,Mxy,Myy)

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j, k
      real*8 Omegalog(imax,jmax),       Bxx(imax,jmax),
     &            Bxy(imax,jmax),       Byy(imax,jmax),
     &           dudx(imax,jmax),      dudy(imax,jmax),
     &           dvdx(imax,jmax),      dvdy(imax,jmax),
     &            Mxx(imax,jmax),       Mxy(imax,jmax),
     &            Myy(imax,jmax),
     &            A(dimen,dimen),    eigvec(dimen,dimen),
     &           Id(dimen,dimen),       Jac(dimen,dimen),
     &     eigvec_t(dimen,dimen),     gradu(dimen,dimen),
     &           MI(dimen,dimen),         M(dimen,dimen),
     &                               eigval(dimen,dimen),
     &     M_lambda(dimen,dimen),         B(dimen,dimen),
     &   eigval_vec(dimen), div, omegal

      ! derivative calculations
      call derx(dudx, ux)
      call dery(dudy, ux)
      call derx(dvdx, uy)
      call dery(dvdy, uy)

      eigval = 0.d0
      Jac    = 0.d0
      Id     = 0.d0
      do k = 1, dimen
        Id(k,k) = 1.d0
      end do

      do j = 1, jmax
        do i = 1, imax

          ! Conformation matrix
          A(1,1) = Txx(i,j) * etapinv + 1.d0
          A(1,2) = Txy(i,j) * etapinv
          A(2,1) =   A(1,2)
          A(2,2) = Tyy(i,j) * etapinv + 1.d0

          call decomp_matrix(A,eigval_vec,eigvec)
!          eigvec(:,1) = - eigvec(:,1)
          eigvec_t    = transpose(eigvec)

          eigval(1,1) = eigval_vec(1)
          eigval(2,2) = eigval_vec(2)

          ! Jacobiana para transformada log
          Jac(1,1)   = 1.d0/eigval_vec(1)
          Jac(2,2)   = 1.d0/eigval_vec(2)

          ! Jacobiana para transformada sqrt
!         Jac(1,1)   = 0.5d0*(1.d0/sqrt(eigval_vec(1)))
!         Jac(2,2)   = 0.5d0*(1.d0/sqrt(eigval_vec(2)))

          ! Jacobiana para transformada 1
!         Jac(1,1)   = 1.d0
!         Jac(2,2)   = 1.d0

          gradu(1,1) = dudx(i,j)
          gradu(1,2) = dudy(i,j)
          gradu(2,1) = dvdx(i,j)
          gradu(2,2) = dvdy(i,j)
          call multi_matrix(eigvec_t,gradu,MI)
          call multi_matrix(MI,eigvec,M)

          div = eigval_vec(2) - eigval_vec(1)
          if (abs(div).lt.1d-14) div = 1.d-14
          omegal = ( M(1,2) * eigval_vec(2) +
     &               M(2,1) * eigval_vec(1) ) / div

          M(1,2) = 0.d0
          M(2,1) = 0.d0
          call multi_matrix(eigvec,M,MI)
          call multi_matrix(MI,eigval,M)
          call multi_matrix(M,Jac,MI)
          call multi_matrix(Mi,eigvec_t,B)
          Bxx(i,j) = B(1,1)
          Bxy(i,j) = B(1,2)
          Byy(i,j) = B(2,2)

          ! Cálculo de M(lambda)
!         call multi_matrix(Id-eigval,Id+alphaG*(eigval-Id),M_lambda)
          call multi_matrix(Id-eigval,Id+alphaG*(eigval-Id),M_lambda)

          call multi_matrix(eigvec,M_lambda,MI)
          call multi_matrix(MI,Jac,M)
          call multi_matrix(M,eigvec_t,MI)
          Mxx(i,j) = MI(1,1)
          Mxy(i,j) = MI(1,2)
          Myy(i,j) = MI(2,2)

          M(1,1)   =   0.d0
          M(1,2)   =   omegal
          M(2,1)   = - omegal
          M(2,2)   =   0.d0
          call multi_matrix(eigvec,M,MI)
          call multi_matrix(MI,eigvec_t,M)
          Omegalog(i,j) = M(1,2)

        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc !
! Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which
! is of size n by n, stored in a physical np by np array. On output, elements of
! a above the diagonal are destroyed. d returns the eigenvalues of a in its first
! n elements. v is a matrix with the same logical and physical dimensions as a,
! whose columns contain, on output, the normalized eigenvectors of a. nrot
! returns the number of Jacobi rotations that were required.

      subroutine decomp_matrix(a,eigval,eigvec)

      implicit none
      include 'par.nn'
      integer i, ip, iq, j
      integer nrot, NMAX
      real*8 a(dimen,dimen),eigval(dimen),eigvec(dimen,dimen)
      parameter (NMAX=500)
      real*8 c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)

      eigval = 0.d0
      ! Initialize to the identity matrix.
      eigvec = 0.d0
      do ip = 1, dimen 
        eigvec(ip,ip) = 1.d0
      end do

      do ip = 1, dimen
        b(ip)      = a(ip,ip) ! Initialize b and d to the diagonal of a.
        eigval(ip) = b(ip)
        z(ip)      = 0.d0 ! This vector will accumulate terms of the form tapq
      end do ! as in equation (11.1.14).

      nrot = 0
      do i = 1, 500
        sm = 0.d0
        do ip = 1, dimen - 1 ! Sum off-diagonal elements. 
          do iq = ip + 1, dimen
            sm = sm + dabs(a(ip,iq))
          end do
        end do

        if (sm.eq.0.d0) return

        if (i.lt.4) then
          tresh = 0.2d0 * sm / dimen**2 
        else
          tresh = 0.d0
        end if

        do ip = 1, dimen - 1
          do iq = ip + 1, dimen
            g = 100.d0 * dabs(a(ip,iq))
            if ((i.gt.4).and.(dabs(eigval(ip))+
     &        g.eq.dabs(eigval(ip))).and.
     &        (dabs(eigval(iq))+g.eq.dabs(eigval(iq)))) then
              a(ip,iq)=0.d0
             else if (dabs(a(ip,iq)).gt.tresh) then 
              h = eigval(iq) - eigval(ip)
              if (dabs(h)+g.eq.dabs(h)) then 
                t = a(ip,iq) / h
               else
                theta = 0.5d0 * h / a(ip,iq)
                t = 1.d0 / (dabs(theta)+dsqrt(1.d0+theta**2))
                if (theta.lt.0.d0) t = - t
              end if

              c          = 1.d0 / sqrt(1.d0 + t**2)
              s          = t * c
              tau        = s / (1.d0 + c)
              h          = t * a(ip,iq)
              z(ip)      = z(ip) - h
              z(iq)      = z(iq) + h
              eigval(ip) = eigval(ip) - h
              eigval(iq) = eigval(iq) + h
              a(ip,iq)   = 0.d0

              do j = 1, ip - 1
                g       = a(j,ip)
                h       = a(j,iq)
                a(j,ip) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
              end do
              do j = ip + 1, iq - 1
                g       = a(ip,j)
                h       = a(j,iq)
                a(ip,j) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
              end do
              do j = iq + 1, dimen 
                g       = a(ip,j)
                h       = a(iq,j)
                a(ip,j) = g-s*(h+g*tau)
                a(iq,j) = h+s*(g-h*tau)
              end do
              do j = 1, dimen
                g            = eigvec(j,ip)
                h            = eigvec(j,iq)
                eigvec(j,ip) = g-s*(h+g*tau)
                eigvec(j,iq) = h+s*(g-h*tau)
              enddo
              nrot = nrot + 1
            end if
          end do
        end do
        do ip = 1, dimen
          b(ip)      = b(ip)+z(ip)
          eigval(ip) = b(ip)
          z(ip)      = 0.d0
        end do
      end do 

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine multi_matrix(A,B,M)

      implicit none
      include 'par.nn'
      integer i, j, k
      real*8 M(dimen,dimen), A(dimen,dimen), B(dimen,dimen)

      do i = 1, dimen
        do j = 1, dimen
          M(i,j) = 0
          do k = 1, dimen
            M(i,j) = M(i,j) + A(i,k) * B(k,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
