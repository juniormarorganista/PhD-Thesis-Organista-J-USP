cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_13(t)
      
      implicit none
      include 'par.nn'
      include 'comm.vare'
      integer i, j
      real*8 x, y, t, r
      r = t
      !write(*,*) '********************'
      do j = 1, jmax
         y = y0 + dble(j-1)*dy
         do i = 1, imax
             x    = x0 + dble(i-1)*dx
         !!! Velocitys 
         uxe(i,j) = 0.0
         
         uye(i,j) = 0.0
         
         wze(i,j) = 0.0
         
        psie(i,j) = 0.0
        
         !!! Tensors
        Txxe(i,j) = 0.0
        
        Txye(i,j) = 0.0
        
        Tyye(i,j) = 0.0
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_13(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 x, y, t, r
      r = t

      do j = 1, jmax
         y    = y0 + dble(j-1)*dy
         do i = 1, imax
            x    = x0 + dble(i-1)*dx

        !!! Velocitys 
        tfwz(i,j) = 0.0

       !!! Tensors
       tfTxx(i,j) = 0.0

       tfTxy(i,j) = 0.0 

       tfTyy(i,j) = 0.0

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_13(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j
      real*8 x, y, t, r
      r = t

      do i = 2, imax - 1
      end do

      do j = 2, jmax - 1
      end do
      
      !!! Boundary: left and Right
      do j = 1, jmax
         y = y0 + dble(j-1)*dy
         do i = 1, imax, imax -  1
            x     = x0 + dble(i-1)*dx

         !!! Velocitys 
         ux(i,j) = 0.0 

         uy(i,j) = 0.0 

      wz(1,:)    = (-85.d0*psi(1,:)+108.d0*psi(2,:)-
     &               27.d0*psi(3,:)+  4.d0*psi(4,:))/(18.d0*dxx)
      wz(imax,:) =(-85.d0*psi(imax,:)+108.d0*psi(imax-1,:)-
     &             27.d0*psi(imax-2,:)+4.d0*psi(imax-3,:))/(18.d0*dxx)

        psi(i,j) = 0.0 

         !!! Tensors
        Txx(i,j) = 0.0

        Txy(i,j) = 0.0 

        Tyy(i,j) = 0.0
         end do
      end do

      !!! Boundary: top and botton
      do j = 1, jmax, jmax - 1
         y = y0 + dble(j-1)*dy
         do i = 1, imax
            x     = x0 + dble(i-1)*dx

         !!! Velocitys 
         ux(i,1) = 0.0 

         ux(i,jmax) = 1.0 

         uy(i,j) = 0.0 

      wz(:,1)    = (-85.d0*psi(:,1)+108.d0*psi(:,2)-
     &               27.d0*psi(:,3)+  4.d0*psi(:,4))/(18.d0*dyy)
      wz(:,jmax) = (11.d0*ux(:,jmax) -18.d0*ux(:,jmax-1)+
     &               9.d0*ux(:,jmax-2)-2.d0*ux(:,jmax-3))/(6.d0*dy)

        psi(i,j) = 0.0 

         !!! Tensors
        Txx(i,j) = 0.0 

        Txy(i,j) = 0.0

        Tyy(i,j) = 0.0 

         end do
      end do
      return
      end
