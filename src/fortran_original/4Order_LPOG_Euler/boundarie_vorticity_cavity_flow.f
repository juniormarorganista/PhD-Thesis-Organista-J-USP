      subroutine boundarie_vorticity_cavity_flow()

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j    
      real*8 dx2, dy2

      dy2 = 18.d0*dy*dy
      dx2 = 18.d0*dx*dx

      do i = 2, imax - 1
        wz(i,1)    = (-85.d0*psi(i,1)+108.d0*psi(i,2)-
     &                 27.d0*psi(i,3)+  4.d0*psi(i,4))/dy2
        wz(i,jmax) = (11.d0*ux(i,jmax) -18.d0*ux(i,jmax-1)+
     &                 9.d0*ux(i,jmax-2)-2.d0*ux(i,jmax-3))/(6.d0*dy)
      end do

      do j = 2, jmax - 1
        wz(1,j)    = (-85.d0*psi(1,j)+108.d0*psi(2,j)-
     &                 27.d0*psi(3,j)+  4.d0*psi(4,j))/dx2
        wz(imax,j) = (-85.d0*psi(imax,j)+108.d0*psi(imax-1,j)-
     &                 27.d0*psi(imax-2,j)+4.d0*psi(imax-3,j))/dx2
      end do
      
      return

      end
