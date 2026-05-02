cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_10(t)
      
      implicit none
      include 'par.nn'
      include 'comm.vare'
      integer i, j
      real*8 x, y, t

      do j = 1, jmax
         y = y0 + dble(j-1)*dy
         do i = 1, imax
             x    = x0 + dble(i-1)*dx

         !!! Velocitys 
         uxe(i,j) = ((-1.0d0 + x)*x*dsin(pii*x)*
     & (pii*y*dcos(pii*y) + dsin(pii*y)))/dexp(at*t)

         uye(i,j) = -((y*(pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     & (-1.0d0 + 2.0d0*x)*dsin(pii*x))*dsin(pii*y))/dexp(at*t))

         wze(i,j) = (2.0d0*(pii*(-1.0d0 + x)*x*dcos(pii*y)*
     & dsin(pii*x) + 
     & y*(pii*(-1.0d0 + 2.0d0*x)*dcos(pii*x) + 
     & (1.0d0 - pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
     & dsin(pii*y)))/dexp(at*t)

        psie(i,j) = ((-1.0d0 + x)*x*y*dsin(pii*x)*
     & dsin(pii*y))/dexp(at*t)

         !!! Tensors
        Txxe(i,j) = -(((-1.0d0 + betann)*c1)/dexp(at*t))

        Txye(i,j) = -(((-1.0d0 + betann)*c2)/dexp(at*t))

        Tyye(i,j) = -(((-1.0d0 + betann)*c3)/dexp(at*t))

         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_10(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 x, y, t

      do j = 1, jmax
         y    = y0 + dble(j-1)*dy
         do i = 1, imax
            x    = x0 + dble(i-1)*dx

        !!! Velocitys 
        tfwz(i,j) = (2*pii**2*Rey*(-1 + x)*x*(-1 + 2*x)*
     & y*dcos(pii*x)**2*
     & dsin(pii*y)*(pii*y*dcos(pii*y) + dsin(pii*y)) + 
     & 2*pii*dcos(pii*x)*
     & (pii**2*Rey*x**2*y*
     & (-1 + x - x**2 + x*dcos(2*pii*y))*dsin(pii*x) - 
     & dexp(at*t)*(-1 + 2*x)*
     & (-4*betann*pii*dcos(pii*y) + 
     & (4*betann*pii**2 - at*Rey)*y*dsin(pii*y))) + 
     & dsin(pii*x)*(2*dexp(at*t)*
     & (pii*(at*Rey*(-1 + x)*x + 
     & betann*(4 - 4*pii**2*(-1 + x)*x))*
     & dcos(pii*y) + (at*Rey + pii**2*
     & (-(at*Rey*(-1 + x)*x) + 
     & 2*betann*(-4 + pii**2*(-1 + x)*x)))*y*
     & dsin(pii*y)) + 
     & pii*Rey*dcos(pii*x)*
     &  (2*(1 + 2*x*(-1 + x + pii**2*x**2))*y*
     &     dsin(pii*y)**2 - 
     &    pii*((-1 + x)**2*x**2 + 
     &       (-1 - 2*(-1 + x)*x)*y**2)*dsin(2*pii*y))
     &  + Rey*(-1 + 2*x)*dsin(pii*x)*
     &  (y - pii**2*(-1 + x)*x*y + 
     &    (-1 - pii**2*(-1 + x)*x)*y*dcos(2*pii*y) + 
     &    pii*(y**2 + (-1 + x)*x*(-1 + pii**2*y**2))*
     &     dsin(2*pii*y))))/(dexp(2*at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) = ((-1 + betann)*(Rey*
     & (-((alphaG*(c1**2 + c2**2) + 
     &         c1*(c1 + c3)*epsylon)*Rey*Wi) + 
     &    c1*dexp(at*t)*(-1 + at*Wi)) + 
     & 2*pii*dcos(pii*x)*
     &  (pii*(-1 + x)*x*
     &     (dexp(at*t) - c1*Rey*Wi*(-1 + xi))*y*
     &     dcos(pii*y) + 
     &    (dexp(at*t)*(-1 + x)*x + 
     &       Rey*Wi*
     &        (-(c1*(-1 + x)*x*(-1 + xi)) + 
     &          c2*(-1 + 2*x)*xi*y))*dsin(pii*y)) + 
     & 2*dsin(pii*x)*
     &  (pii*(-(c2*Rey*Wi*(-1 + x)*x*(-2 + xi)) + 
     &       (-1 + 2*x)*
     &        (dexp(at*t) - c1*Rey*Wi*(-1 + xi))*y)*
     &     dcos(pii*y) + 
     &    (dexp(at*t)*(-1 + 2*x) + 
     &       Rey*Wi*
     &        (-(c1*(-1 + 2*x)*(-1 + xi)) + 
     &          c2*(-(pii**2*(-1 + x)*x) + xi)*y))*
     &     dsin(pii*y))))/(dexp(2*at*t)*Rey)

       tfTxy(i,j) = ((-1 + betann)*(c2*Rey*
     &  (-((c1 + c3)*(alphaG + epsylon)*Rey*Wi) + 
     &    dexp(at*t)*(-1 + at*Wi)) + 
     & pii*(-1 + x)*x*
     &  (2*dexp(at*t) + Rey*Wi*(2*c3 - (c1 + c3)*xi))*
     &  dcos(pii*y)*dsin(pii*x) + 
     & y*(-(pii*(-1 + 2*x)*
     &       (2*dexp(at*t) + 
     &         Rey*Wi*(2*c1 - (c1 + c3)*xi))*dcos(pii*x)
     &       ) + (-2*dexp(at*t) + 
     &       Rey*Wi*
     &        (c3*(-(pii**2*(-1 + x)*x) + xi) + 
     &          c1*(-2 + pii**2*(-1 + x)*x + xi)))*
     &     dsin(pii*x))*dsin(pii*y)))/(dexp(2*at*t)*Rey)

       tfTyy(i,j) = ((-1 + betann)*(Rey*
     &  (-((alphaG*(c2**2 + c3**2) + 
     &         c3*(c1 + c3)*epsylon)*Rey*Wi) + 
     &    c3*dexp(at*t)*(-1 + at*Wi)) + 
     & 2*pii*dcos(pii*x)*
     &  (-(pii*(-1 + x)*x*
     &       (dexp(at*t) - c3*Rey*Wi*(-1 + xi))*y*
     &       dcos(pii*y)) + 
     &    (-(dexp(at*t)*(-1 + x)*x) + 
     &       Rey*Wi*
     &        (c3*(-1 + x)*x*(-1 + xi) + 
     &          c2*(-1 + 2*x)*(-2 + xi)*y))*dsin(pii*y))
     &   + 2*dsin(pii*x)*
     &  (pii*(-(c2*Rey*Wi*(-1 + x)*x*xi) - 
     &       (-1 + 2*x)*
     &        (dexp(at*t) - c3*Rey*Wi*(-1 + xi))*y)*
     &     dcos(pii*y) + 
     &    (dexp(at*t)*(1 - 2*x) + 
     &       Rey*Wi*
     &        (c3*(-1 + 2*x)*(-1 + xi) + 
     &          c2*(-2 + pii**2*(-1 + x)*x + xi)*y))*
     &     dsin(pii*y))))/(dexp(2*at*t)*Rey)

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_10(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j
      real*8 x, y, t
      
      !!! Boundary: left and Right
      do j = 1, jmax
         y = y0 + dble(j-1)*dy
         do i = 1, imax, imax -  1
            x     = x0 + dble(i-1)*dx

         !!! Velocitys 
         ux(i,j) = 0.0d0 

         uy(i,j) = 0.0d0

c       wz(1,j)    = (-85.d0*psi(1,j)+108.d0*psi(2,j)-
c    &                 27.d0*psi(3,j)+  4.d0*psi(4,j))/dxx
c       wz(imax,j) = (-85.d0*psi(imax,j)+108.d0*psi(imax-1,j)-
c    &                 27.d0*psi(imax-2,j)+4.d0*psi(imax-3,j))/dxx

c        wz(i,j) = (2.0d0*(pii*(-1.0d0 + x)*x*dcos(pii*y)*
c    & dsin(pii*x) + 
c    & y*(pii*(-1.0d0 + 2.0d0*x)*dcos(pii*x) + 
c    & (1.0d0 - pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
c    & dsin(pii*y)))/dexp(at*t)

        psi(i,j) = ((-1.0d0 + x)*x*y*dsin(pii*x)*
     & dsin(pii*y))/dexp(at*t)

         !!! Tensors
        Txx(i,j) = -(((-1.0d0 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1.0d0 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1.0d0 + betann)*c3)/dexp(at*t))


         end do
      end do
      !!! Boundary: top and botton
      do j = 1, jmax, jmax - 1
         y = y0 + dble(j-1)*dy
         do i = 1, imax
            x     = x0 + dble(i-1)*dx

         !!! Velocitys 
         ux(i,1) = 0.0d0 

         ux(i,jmax) = ((-1.0d0 + x)*x*dsin(pii*x)*
     &    (pii*y*dcos(pii*y) + dsin(pii*y)))/dexp(at*t)

         uy(i,j) = 0.0d0 

c        wz(i,1)    = (-85.d0*psi(i,1)+108.d0*psi(i,2)-
c    &                 27.d0*psi(i,3)+  4.d0*psi(i,4))/dyy
c        wz(i,jmax) = (11.d0*ux(i,jmax) -18.d0*ux(i,jmax-1)+
c    &                 9.d0*ux(i,jmax-2)-2.d0*ux(i,jmax-3))/(6.d0*dy)

c        wz(i,j) = (2.0d0*(pii*(-1.0d0 + x)*x*dcos(pii*y)*
c    & dsin(pii*x) + 
c    & y*(pii*(-1.0d0 + 2.0d0*x)*dcos(pii*x) + 
c    & (1.0d0 - pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
c    & dsin(pii*y)))/dexp(at*t)

        psi(i,j) = ((-1.0d0 + x)*x*y*dsin(pii*x)*
     & dsin(pii*y))/dexp(at*t)

         !!! Tensors
        Txx(i,j) = -(((-1.0d0 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1.0d0 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1.0d0 + betann)*c3)/dexp(at*t))

         end do
      end do
      
      call boundarie_vorticity_cavity_flow 
      
      return
      end
