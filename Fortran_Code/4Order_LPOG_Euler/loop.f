ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 loop calculations                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_ux

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 dpsidy(imax,jmax)
       
      call dery(dpsidy, psi)

      ux(2:imax-1,2:jmax-1) = dpsidy(2:imax-1,2:jmax-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_uy
      
      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 dpsidx(imax,jmax)
      
      call derx(dpsidx, psi)

      uy(2:imax-1,2:jmax-1) = - dpsidx(2:imax-1,2:jmax-1)
       
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Applies the filter to remove short length scales at the x-direction
!! for vorticity and non-Newtonian extra-tensor components
      subroutine filter_trid_x

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 afilx(imax), bfilx(imax), cfilx(imax),
     &        rhswz(imax,jmax),   rhsTxx(imax,jmax), 
     &       rhsTxy(imax,jmax),   rhsTyy(imax,jmax),
     &     rhsPsixx(imax,jmax), rhsPsixy(imax,jmax), 
     &     rhsPsiyy(imax,jmax)
      common/filtx/ afilx, bfilx, cfilx

      call rhs_tridfx_wz(rhswz)
      call tridx(afilx, bfilx, cfilx, rhswz)
!     wz(2:imax-1,2:jmax-1) = rhswz(2:imax-1,2:jmax-1)
!     wz = rhswz

      if (ty_sim.eq.1) then
        call rhs_tridfx_T(rhsTxx, rhsTxy, rhsTyy)
        call tridx(afilx, bfilx, cfilx, rhsTxx)
        call tridx(afilx, bfilx, cfilx, rhsTxy)
        call tridx(afilx, bfilx, cfilx, rhsTyy)
         wz(2:imax-1,2:jmax-1) =  rhswz(2:imax-1,2:jmax-1)
        Txx(2:imax-1,2:jmax-1) = rhsTxx(2:imax-1,2:jmax-1)
        Txy(2:imax-1,2:jmax-1) = rhsTxy(2:imax-1,2:jmax-1)
        Tyy(2:imax-1,2:jmax-1) = rhsTyy(2:imax-1,2:jmax-1)
!       Txx = rhsTxx
!       Txy = rhsTxy
!       Tyy = rhsTyy
        return
       elseif (ty_sim.eq.2) then
        call rhs_tridfx_Psi(rhsPsixx, rhsPsixy, rhsPsiyy)
        call tridx(afilx, bfilx, cfilx, rhsPsixx)
        call tridx(afilx, bfilx, cfilx, rhsPsixy)
        call tridx(afilx, bfilx, cfilx, rhsPsiyy)
           wz(2:imax-1,2:jmax-1) =    rhswz(2:imax-1,2:jmax-1)
        Psixx(2:imax-1,2:jmax-1) = rhsPsixx(2:imax-1,2:jmax-1)
        Psixy(2:imax-1,2:jmax-1) = rhsPsixy(2:imax-1,2:jmax-1)
        Psiyy(2:imax-1,2:jmax-1) = rhsPsiyy(2:imax-1,2:jmax-1)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Applies the filter to remove short length scales at the y-direction 
!! for vorticity and non-Newtonian extra-tensor components
      subroutine filter_trid_y

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 afily(jmax), bfily(jmax), cfily(jmax),
     &        rhswz(imax,jmax),   rhsTxx(imax,jmax), 
     &       rhsTxy(imax,jmax),   rhsTyy(imax,jmax),
     &     rhsPsixx(imax,jmax), rhsPsixy(imax,jmax), 
     &     rhsPsiyy(imax,jmax)
      common/filty/ afily, bfily, cfily

      call rhs_tridfy_wz(rhswz)
      call tridy(afily, bfily, cfily, rhswz)
!     wz(2:imax-1,2:jmax-1) = rhswz(2:imax-1,2:jmax-1)
!     wz = rhswz

      if (ty_sim.eq.1) then
        call rhs_tridfy_T(rhsTxx, rhsTxy, rhsTyy)
        call tridy(afily, bfily, cfily, rhsTxx)
        call tridy(afily, bfily, cfily, rhsTxy)
        call tridy(afily, bfily, cfily, rhsTyy)
         wz(2:imax-1,2:jmax-1) =  rhswz(2:imax-1,2:jmax-1)
        Txx(2:imax-1,2:jmax-1) = rhsTxx(2:imax-1,2:jmax-1)
        Txy(2:imax-1,2:jmax-1) = rhsTxy(2:imax-1,2:jmax-1)
        Tyy(2:imax-1,2:jmax-1) = rhsTyy(2:imax-1,2:jmax-1)
!       Txx = rhsTxx
!       Txy = rhsTxy
!       Tyy = rhsTyy
        return
       elseif (ty_sim.eq.2) then
        call rhs_tridfy_Psi(rhsPsixx, rhsPsixy, rhsPsiyy)
        call tridy(afily, bfily, cfily, rhsPsixx)
        call tridy(afily, bfily, cfily, rhsPsixy)
        call tridy(afily, bfily, cfily, rhsPsiyy)
           wz(2:imax-1,2:jmax-1) =    rhswz(2:imax-1,2:jmax-1)
        Psixx(2:imax-1,2:jmax-1) = rhsPsixx(2:imax-1,2:jmax-1)
        Psixy(2:imax-1,2:jmax-1) = rhsPsixy(2:imax-1,2:jmax-1)
        Psiyy(2:imax-1,2:jmax-1) = rhsPsiyy(2:imax-1,2:jmax-1)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>Calculates the right hand side matrix of the application of the filter used to remove the short length scales for the vorticity \f$ \omega_{z} \f$ and the tensors \f$ T^{xx} \f$, \f$ T^{xy} \f$ and \f$ T^{yy} \f$, at the x-direction
!! - Near the boundary nodes
!! \f{eqnarray*}{\displaystyle \widehat{f}_1 = \frac{15}{16}f_{1} + \frac{1}{16}\big(4 f_{2} - 6f_{3} - 4f_{4} - f_{5}\big)\f} 
!! \f{eqnarray*}{\displaystyle \widehat{f}_2 = \frac{3}{4}f_{2} + \frac{1}{16}\big(f_{1} + 6f_{3} - 4f_{4} + f_{5}\big)\f} 
!! \f{eqnarray*}{\displaystyle \widehat{f}_3 = \frac{5}{8}f_{3} + \frac{1}{16}\big(- f_{1} + 4f_{2} - 4f_{4} - f_{5}\big)\f}
!! - At the center nodes 
!! \f{eqnarray*}{\displaystyle \widehat{f}_i = \left(\frac{11}{16} + \frac{10}{16}\alpha\right) f_{i} + \left(\frac{15}{64} + \frac{34}{64}\alpha\right) \left(f_{i+1}+f_{i-1}\right) + \left(\frac{-3}{32} + \frac{6}{32}\alpha\right) \left(f_{i+2}+f_{i-2}\right) + \left(\frac{1}{64} - \frac{2}{64}\alpha\right) \left(f_{i+3}+f_{i-3}\right) \f}
!! where \f$ \widehat{f}_i \f$ represents the filtered values at node \f$ x_{i} \f$. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfx_wz(rhswz)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8  rhswz(imax,jmax)
      do j = 1, jmax
        rhswz(1,j) = (     15.d0 * wz(1,j) +  4.d0 * wz(2,j)
     &                   -  6.d0 * wz(3,j) +  4.d0 * wz(4,j)
     &                   -         wz(5,j) ) / 16.d0
        rhswz(2,j) = (             wz(1,j) + 12.d0 * wz(2,j)
     &                   +  6.d0 * wz(3,j) -  4.d0 * wz(4,j)
     &                   +         wz(5,j) ) / 16.d0
        rhswz(3,j) = (   -         wz(1,j) +  4.d0 * wz(2,j)
     &                   + 10.d0 * wz(3,j) +  4.d0 * wz(4,j)
     &                   -         wz(5,j) ) / 16.d0
        do i = 4, imax - 3
          rhswz(i,j)  = af *   wz(i,j)
     &                + bf * ( wz(i+1,j) + wz(i-1,j) )
     &                + cf * ( wz(i+2,j) + wz(i-2,j) )
     &                + df * ( wz(i+3,j) + wz(i-3,j) )
        end do
        rhswz(imax-2,j) = (  -         wz(imax,j)
     &                       +  4.d0 * wz(imax-1,j)
     &                       + 10.d0 * wz(imax-2,j)
     &                       +  4.d0 * wz(imax-3,j)
     &                       -         wz(imax-4,j) ) / 16.d0 
        rhswz(imax-1,j) = (            wz(imax,j)
     &                       + 12.d0 * wz(imax-1,j)
     &                       +  6.d0 * wz(imax-2,j)
     &                       -  4.d0 * wz(imax-3,j)
     &                       +         wz(imax-4,j) ) / 16.d0 
        rhswz(imax,j)   = (    15.d0 * wz(imax,j)
     &                       +  4.d0 * wz(imax-1,j)
     &                       -  6.d0 * wz(imax-2,j)
     &                       +  4.d0 * wz(imax-3,j)
     &                       -         wz(imax-4,j) ) / 16.d0 
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfx_T(rhsTxx, rhsTxy, rhsTyy)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 rhsTxx(imax,jmax), rhsTxy(imax,jmax), 
     &       rhsTyy(imax,jmax)

      do j = 1, jmax

        rhsTxx(1,j) = (     15.d0 * Txx(1,j) +  4.d0 * Txx(2,j)
     &                    -  6.d0 * Txx(3,j) +  4.d0 * Txx(4,j)
     &                    -         Txx(5,j) ) / 16.d0
        rhsTxx(2,j) = (             Txx(1,j) + 12.d0 * Txx(2,j)
     &                    +  6.d0 * Txx(3,j) -  4.d0 * Txx(4,j)
     &                    +         Txx(5,j) ) / 16.d0 
        rhsTxx(3,j) = ( -           Txx(1,j) +  4.d0 * Txx(2,j)
     &                    + 10.d0 * Txx(3,j) +
     &                       4.d0 * Txx(4,j) -         Txx(5,j) )
     &                    / 16.d0 

        rhsTxy(1,j) = (     15.d0 * Txy(1,j) +  4.d0 * Txy(2,j)
     &                    -  6.d0 * Txy(3,j) +  4.d0 * Txy(4,j)
     &                    -         Txy(5,j) ) / 16.d0
        rhsTxy(2,j) = (             Txy(1,j) + 12.d0 * Txy(2,j)
     &                    +  6.d0 * Txy(3,j) -  4.d0 * Txy(4,j)
     &                    +         Txy(5,j) ) / 16.d0 
        rhsTxy(3,j) = ( -           Txy(1,j) +  4.d0 * Txy(2,j)
     &                    + 10.d0 * Txy(3,j) +
     &                       4.d0 * Txy(4,j) -         Txy(5,j) )
     &                    / 16.d0 

        rhsTyy(1,j) = (     15.d0 * Tyy(1,j) +  4.d0 * Tyy(2,j)
     &                    -  6.d0 * Tyy(3,j) +  4.d0 * Tyy(4,j)
     &                    -         Tyy(5,j) ) / 16.d0
        rhsTyy(2,j) = (             Tyy(1,j) + 12.d0 * Tyy(2,j)
     &                    +  6.d0 * Tyy(3,j) -  4.d0 * Tyy(4,j)
     &                    +         Tyy(5,j) ) / 16.d0 
        rhsTyy(3,j) = ( -           Tyy(1,j) +  4.d0 * Tyy(2,j)
     &                    + 10.d0 * Tyy(3,j) +
     &                       4.d0 * Tyy(4,j) -         Tyy(5,j) )
     &                    / 16.d0 


        do i = 4, imax - 3
          rhsTxx(i,j) = af *   Txx(i,j)
     &                + bf * ( Txx(i+1,j) + Txx(i-1,j) )
     &                + cf * ( Txx(i+2,j) + Txx(i-2,j) )
     &                + df * ( Txx(i+3,j) + Txx(i-3,j) )
          rhsTxy(i,j) = af *   Txy(i,j)
     &                + bf * ( Txy(i+1,j) + Txy(i-1,j) )
     &                + cf * ( Txy(i+2,j) + Txy(i-2,j) )
     &                + df * ( Txy(i+3,j) + Txy(i-3,j) )
          rhsTyy(i,j) = af *   Tyy(i,j)
     &                + bf * ( Tyy(i+1,j) + Tyy(i-1,j) )
     &                + cf * ( Tyy(i+2,j) + Tyy(i-2,j) )
     &                + df * ( Tyy(i+3,j) + Tyy(i-3,j) )
        end do

        rhsTxx(imax-2,j) = ( -         Txx(imax,j)
     &                       +  4.d0 * Txx(imax-1,j)
     &                       + 10.d0 * Txx(imax-2,j)
     &                       +  4.d0 * Txx(imax-3,j)
     &                       -         Txx(imax-4,j) ) / 16.d0 
        rhsTxx(imax-1,j) = (           Txx(imax,j)
     &                       + 12.d0 * Txx(imax-1,j)
     &                       +  6.d0 * Txx(imax-2,j)
     &                       -  4.d0 * Txx(imax-3,j)
     &                       +         Txx(imax-4,j) ) / 16.d0 
        rhsTxx(imax,j)   = (   15.d0 * Txx(imax,j)
     &                       +  4.d0 * Txx(imax-1,j)
     &                       -  6.d0 * Txx(imax-2,j)
     &                       +  4.d0 * Txx(imax-3,j)
     &                       -         Txx(imax-4,j) ) / 16.d0 

        rhsTxy(imax-2,j) = ( -         Txy(imax,j)
     &                       +  4.d0 * Txy(imax-1,j)
     &                       + 10.d0 * Txy(imax-2,j)
     &                       +  4.d0 * Txy(imax-3,j)
     &                       -         Txy(imax-4,j) ) / 16.d0 
        rhsTxy(imax-1,j) = (           Txy(imax,j)
     &                       + 12.d0 * Txy(imax-1,j)
     &                       +  6.d0 * Txy(imax-2,j)
     &                       -  4.d0 * Txy(imax-3,j)
     &                       +         Txy(imax-4,j) ) / 16.d0 
        rhsTxy(imax,j)   = (   15.d0 * Txy(imax,j)
     &                       +  4.d0 * Txy(imax-1,j)
     &                       -  6.d0 * Txy(imax-2,j)
     &                       +  4.d0 * Txy(imax-3,j)
     &                       -         Txy(imax-4,j) ) / 16.d0 

        rhsTyy(imax-2,j) = ( -         Tyy(imax,j)
     &                       +  4.d0 * Tyy(imax-1,j)
     &                       + 10.d0 * Tyy(imax-2,j)
     &                       +  4.d0 * Tyy(imax-3,j)
     &                       -         Tyy(imax-4,j) ) / 16.d0 
        rhsTyy(imax-1,j) = (           Tyy(imax,j)
     &                       + 12.d0 * Tyy(imax-1,j)
     &                       +  6.d0 * Tyy(imax-2,j)
     &                       -  4.d0 * Tyy(imax-3,j)
     &                       +         Tyy(imax-4,j) ) / 16.d0 
        rhsTyy(imax,j)   = (   15.d0 * Tyy(imax,j)
     &                       +  4.d0 * Tyy(imax-1,j)
     &                       -  6.d0 * Tyy(imax-2,j)
     &                       +  4.d0 * Tyy(imax-3,j)
     &                       -         Tyy(imax-4,j) ) / 16.d0 

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfx_Psi(rhsPsixx, rhsPsixy, rhsPsiyy)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 rhsPsixx(imax,jmax), rhsPsixy(imax,jmax), 
     &       rhsPsiyy(imax,jmax)

      do j = 1, jmax

        rhsPsixx(1,j) = (   15.d0 * Psixx(1,j) +  4.d0 * Psixx(2,j)
     &                    -  6.d0 * Psixx(3,j) +  4.d0 * Psixx(4,j)
     &                    -         Psixx(5,j) ) / 16.d0
        rhsPsixx(2,j) = (           Psixx(1,j) + 12.d0 * Psixx(2,j)
     &                    +  6.d0 * Psixx(3,j) -  4.d0 * Psixx(4,j)
     &                    +         Psixx(5,j) ) / 16.d0 
        rhsPsixx(3,j) = ( -         Psixx(1,j) +  4.d0 * Psixx(2,j)
     &                    + 10.d0 * Psixx(3,j) +
     &                       4.d0 * Psixx(4,j) -         Psixx(5,j) )
     &                    / 16.d0 

        rhsPsixy(1,j) = (   15.d0 * Psixy(1,j) +  4.d0 * Psixy(2,j)
     &                    -  6.d0 * Psixy(3,j) +  4.d0 * Psixy(4,j)
     &                    -         Psixy(5,j) ) / 16.d0
        rhsPsixy(2,j) = (           Psixy(1,j) + 12.d0 * Psixy(2,j)
     &                    +  6.d0 * Psixy(3,j) -  4.d0 * Psixy(4,j)
     &                    +         Psixy(5,j) ) / 16.d0 
        rhsPsixy(3,j) = (           Psixy(1,j) +  4.d0 * Psixy(2,j)
     &                    + 10.d0 * Psixy(3,j) +
     &                       4.d0 * Psixy(4,j) -         Psixy(5,j) )
     &                    / 16.d0 

        rhsPsiyy(1,j) = (   15.d0 * Psiyy(1,j) +  4.d0 * Psiyy(2,j)
     &                    -  6.d0 * Psiyy(3,j) +  4.d0 * Psiyy(4,j)
     &                    -         Psiyy(5,j) ) / 16.d0
        rhsPsiyy(2,j) = (           Psiyy(1,j) + 12.d0 * Psiyy(2,j)
     &                    +  6.d0 * Psiyy(3,j) -  4.d0 * Psiyy(4,j)
     &                    +         Psiyy(5,j) ) / 16.d0 
        rhsPsiyy(3,j) = ( -         Psiyy(1,j) +  4.d0 * Psiyy(2,j)
     &                    + 10.d0 * Psiyy(3,j) +
     &                       4.d0 * Psiyy(4,j) -         Psiyy(5,j) )
     &                    / 16.d0 


        do i = 4, imax - 3
          rhsPsixx(i,j) = af * Psixx(i  ,j)
     &                + bf * ( Psixx(i+1,j) + Psixx(i-1,j) )
     &                + cf * ( Psixx(i+2,j) + Psixx(i-2,j) )
     &                + df * ( Psixx(i+3,j) + Psixx(i-3,j) )
          rhsPsixy(i,j) = af * Psixy(i  ,j)
     &                + bf * ( Psixy(i+1,j) + Psixy(i-1,j) )
     &                + cf * ( Psixy(i+2,j) + Psixy(i-2,j) )
     &                + df * ( Psixy(i+3,j) + Psixy(i-3,j) )
          rhsPsiyy(i,j) = af * Psiyy(i  ,j)
     &                + bf * ( Psiyy(i+1,j) + Psiyy(i-1,j) )
     &                + cf * ( Psiyy(i+2,j) + Psiyy(i-2,j) )
     &                + df * ( Psiyy(i+3,j) + Psiyy(i-3,j) )
        end do

        rhsPsixx(imax-2,j) = ( -       Psixx(imax  ,j)
     &                       +  4.d0 * Psixx(imax-1,j)
     &                       + 10.d0 * Psixx(imax-2,j)
     &                       +  4.d0 * Psixx(imax-3,j)
     &                       -         Psixx(imax-4,j) ) / 16.d0 
        rhsPsixx(imax-1,j) = (         Psixx(imax  ,j)
     &                       + 12.d0 * Psixx(imax-1,j)
     &                       +  6.d0 * Psixx(imax-2,j)
     &                       -  4.d0 * Psixx(imax-3,j)
     &                       +         Psixx(imax-4,j) ) / 16.d0 
        rhsPsixx(imax,j)   = ( 15.d0 * Psixx(imax  ,j)
     &                       +  4.d0 * Psixx(imax-1,j)
     &                       -  6.d0 * Psixx(imax-2,j)
     &                       +  4.d0 * Psixx(imax-3,j)
     &                       -         Psixx(imax-4,j) ) / 16.d0 

        rhsPsixy(imax-2,j) = ( -       Psixy(imax  ,j)
     &                       +  4.d0 * Psixy(imax-1,j)
     &                       + 10.d0 * Psixy(imax-2,j)
     &                       +  4.d0 * Psixy(imax-3,j)
     &                       -         Psixy(imax-4,j) ) / 16.d0 
        rhsPsixy(imax-1,j) = (         Psixy(imax  ,j)
     &                       + 12.d0 * Psixy(imax-1,j)
     &                       +  6.d0 * Psixy(imax-2,j)
     &                       -  4.d0 * Psixy(imax-3,j)
     &                       +         Psixy(imax-4,j) ) / 16.d0 
        rhsPsixy(imax,j)   = ( 15.d0 * Psixy(imax  ,j)
     &                       +  4.d0 * Psixy(imax-1,j)
     &                       -  6.d0 * Psixy(imax-2,j)
     &                       +  4.d0 * Psixy(imax-3,j)
     &                       -         Psixy(imax-4,j) ) / 16.d0 

        rhsPsiyy(imax-2,j) = ( -       Psiyy(imax  ,j)
     &                       +  4.d0 * Psiyy(imax-1,j)
     &                       + 10.d0 * Psiyy(imax-2,j)
     &                       +  4.d0 * Psiyy(imax-3,j)
     &                       -         Psiyy(imax-4,j) ) / 16.d0 
        rhsPsiyy(imax-1,j) = (         Psiyy(imax,j)
     &                       + 12.d0 * Psiyy(imax-1,j)
     &                       +  6.d0 * Psiyy(imax-2,j)
     &                       -  4.d0 * Psiyy(imax-3,j)
     &                       +         Psiyy(imax-4,j) ) / 16.d0 
        rhsPsiyy(imax,j)   = ( 15.d0 * Psiyy(imax,j)
     &                       +  4.d0 * Psiyy(imax-1,j)
     &                       -  6.d0 * Psiyy(imax-2,j)
     &                       +  4.d0 * Psiyy(imax-3,j)
     &                       -         Psiyy(imax-4,j) ) / 16.d0 

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>Calculates the right hand side matrix of the application of the filter used to remove the short length scales for the vorticity \f$ \omega_{z} \f$ and the tensors \f$ T^{xx} \f$, \f$ T^{xy} \f$ and \f$ T^{yy} \f$, at the y-direction
!! - Near the boundary nodes
!! \f{eqnarray*}{\displaystyle \widehat{f}_1 = \frac{15}{16}f_{1} + \frac{1}{16}\big(4 f_{2} - 6f_{3} - 4f_{4} - f_{5}\big)\f} 
!! \f{eqnarray*}{\displaystyle \widehat{f}_2 = \frac{3}{4}f_{2} + \frac{1}{16}\big(f_{1} + 6f_{3} - 4f_{4} + f_{5}\big)\f} 
!! \f{eqnarray*}{\displaystyle \widehat{f}_3 = \frac{5}{8}f_{3} + \frac{1}{16}\big(- f_{1} + 4f_{2} - 4f_{4} - f_{5}\big)\f}
!! - At the center nodes 
!! \f{eqnarray*}{\displaystyle \widehat{f}_j = \left(\frac{11}{16} + \frac{10}{16}\alpha\right) f_{j} + \left(\frac{15}{64} + \frac{34}{64}\alpha\right) \left(f_{j+1}+f_{j-1}\right) + \left(\frac{-3}{32} + \frac{6}{32}\alpha\right) \left(f_{j+2}+f_{j-2}\right) + \left(\frac{1}{64} - \frac{2}{64}\alpha\right) \left(f_{j+3}+f_{j-3}\right) \f}
!! where \f$ \widehat{f}_j \f$ represents the filtered values at node \f$ y_{j} \f$. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfy_wz(rhswz)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 rhswz(imax,jmax)

      do i = 2, imax

        rhswz(i,1) = (     15.d0 * wz(i,1) +  4.d0 * wz(i,2)
     &                   -  6.d0 * wz(i,3) +  4.d0 * wz(i,4)
     &                   -         wz(i,5) ) / 16.d0
        rhswz(i,2) = (             wz(i,1) + 12.d0 * wz(i,2)
     &                   +  6.d0 * wz(i,3) -  4.d0 * wz(i,4)
     &                   +         wz(i,5) ) / 16.d0
        rhswz(i,3) = (   -         wz(i,1) +  4.d0 * wz(i,2)
     &                   + 10.d0 * wz(i,3) +  4.d0 * wz(i,4)
     &                   -         wz(i,5) ) / 16.d0

        do j = 4, jmax - 3
          rhswz(i,j)  = af *    wz(i,j)
     &                + bf *  ( wz(i,j+1) + wz(i,j-1) )
     &                + cf *  ( wz(i,j+2) + wz(i,j-2) )
     &                + df *  ( wz(i,j+3) + wz(i,j-3) )
        end do

        rhswz(i,jmax-2) = (  -         wz(i,jmax)
     &                       +  4.d0 * wz(i,jmax-1)
     &                       + 10.d0 * wz(i,jmax-2)
     &                       +  4.d0 * wz(i,jmax-3)
     &                       -         wz(i,jmax-4) ) / 16.d0 
        rhswz(i,jmax-1) = (            wz(i,jmax)
     &                       + 12.d0 * wz(i,jmax-1)
     &                       +  6.d0 * wz(i,jmax-2)
     &                       -  4.d0 * wz(i,jmax-3)
     &                       +         wz(i,jmax-4) ) / 16.d0 
        rhswz(i,jmax)   = (    15.d0 * wz(i,jmax)
     &                       +  4.d0 * wz(i,jmax-1)
     &                       -  6.d0 * wz(i,jmax-2)
     &                       +  4.d0 * wz(i,jmax-3)
     &                       -         wz(i,jmax-4) ) / 16.d0 

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfy_T(rhsTxx, rhsTxy, rhsTyy)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8 rhsTxx(imax,jmax), rhsTxy(imax,jmax), 
     &       rhsTyy(imax,jmax)

      do i = 2, imax

        rhsTxx(i,1) = (     15.d0 * Txx(i,1) +  4.d0 * Txx(i,2)
     &                    -  6.d0 * Txx(i,3) +  4.d0 * Txx(i,4)
     &                    -         Txx(i,5) ) / 16.d0
        rhsTxx(i,2) = (             Txx(i,1) + 12.d0 * Txx(i,2)
     &                    +  6.d0 * Txx(i,3) -  4.d0 * Txx(i,4)
     &                    +         Txx(i,5) ) / 16.d0 
        rhsTxx(i,3) = (   -         Txx(i,1) +  4.d0 * Txx(i,2)
     &                    + 10.d0 * Txx(i,3) +
     &                       4.d0 * Txx(i,4) -         Txx(i,5) )
     &                    / 16.d0 

        rhsTxy(i,1) = (     15.d0 * Txy(i,1) +  4.d0 * Txy(i,2)
     &                    -  6.d0 * Txy(i,3) +  4.d0 * Txy(i,4)
     &                    -         Txy(i,5) ) / 16.d0
        rhsTxy(i,2) = (             Txy(i,1) + 12.d0 * Txy(i,2)
     &                    +  6.d0 * Txy(i,3) -  4.d0 * Txy(i,4)
     &                    +         Txy(i,5) ) / 16.d0 
        rhsTxy(i,3) = (   -         Txy(i,1) +  4.d0 * Txy(i,2)
     &                    + 10.d0 * Txy(i,3) +
     &                       4.d0 * Txy(i,4) -         Txy(i,5) )
     &                    / 16.d0 

        rhsTyy(i,1) = (     15.d0 * Tyy(i,1) +  4.d0 * Tyy(i,2)
     &                    -  6.d0 * Tyy(i,3) +  4.d0 * Tyy(i,4)
     &                    -         Tyy(i,5) ) / 16.d0
        rhsTyy(i,2) = (             Tyy(i,1) + 12.d0 * Tyy(i,2)
     &                    +  6.d0 * Tyy(i,3) -  4.d0 * Tyy(i,4)
     &                    +         Tyy(i,5) ) / 16.d0 
        rhsTyy(i,3) = (   -         Tyy(i,1) +  4.d0 * Tyy(i,2)
     &                    + 10.d0 * Tyy(i,3) +
     &                       4.d0 * Tyy(i,4) -         Tyy(i,5) )
     &                    / 16.d0 


        do j = 4, jmax - 3
          rhsTxx(i,j) = af *   Txx(i,j)
     &                + bf * ( Txx(i,j+1) + Txx(i,j-1) )
     &                + cf * ( Txx(i,j+2) + Txx(i,j-2) )
     &                + df * ( Txx(i,j+3) + Txx(i,j-3) )
          rhsTxy(i,j) = af *   Txy(i,j)
     &                + bf * ( Txy(i,j+1) + Txy(i,j-1) )
     &                + cf * ( Txy(i,j+2) + Txy(i,j-2) )
     &                + df * ( Txy(i,j+3) + Txy(i,j-3) )
          rhsTyy(i,j) = af *   Tyy(i,j)
     &                + bf * ( Tyy(i,j+1) + Tyy(i,j-1) )
     &                + cf * ( Tyy(i,j+2) + Tyy(i,j-2) )
     &                + df * ( Tyy(i,j+3) + Tyy(i,j-3) )
        end do

        rhsTxx(i,jmax-2) = ( -         Txx(i,jmax)
     &                       +  4.d0 * Txx(i,jmax-1)
     &                       + 10.d0 * Txx(i,jmax-2)
     &                       +  4.d0 * Txx(i,jmax-3)
     &                       -         Txx(i,jmax-4) ) / 16.d0 
        rhsTxx(i,jmax-1) = (           Txx(i,jmax)
     &                       + 12.d0 * Txx(i,jmax-1)
     &                       +  6.d0 * Txx(i,jmax-2)
     &                       -  4.d0 * Txx(i,jmax-3)
     &                       +         Txx(i,jmax-4) ) / 16.d0 
        rhsTxx(i,jmax)   = (   15.d0 * Txx(i,jmax)
     &                       +  4.d0 * Txx(i,jmax-1)
     &                       -  6.d0 * Txx(i,jmax-2)
     &                       +  4.d0 * Txx(i,jmax-3)
     &                       -         Txx(i,jmax-4) ) / 16.d0 

        rhsTxy(i,jmax-2) = ( -         Txy(i,jmax)
     &                       +  4.d0 * Txy(i,jmax-1)
     &                       + 10.d0 * Txy(i,jmax-2)
     &                       +  4.d0 * Txy(i,jmax-3)
     &                       -         Txy(i,jmax-4) ) / 16.d0 
        rhsTxy(i,jmax-1) = (           Txy(i,jmax)
     &                       + 12.d0 * Txy(i,jmax-1)
     &                       +  6.d0 * Txy(i,jmax-2)
     &                       -  4.d0 * Txy(i,jmax-3)
     &                       +         Txy(i,jmax-4) ) / 16.d0 
        rhsTxy(i,jmax)   = (   15.d0 * Txy(i,jmax)
     &                       +  4.d0 * Txy(i,jmax-1)
     &                       -  6.d0 * Txy(i,jmax-2)
     &                       +  4.d0 * Txy(i,jmax-3)
     &                       -         Txy(i,jmax-4) ) / 16.d0 

        rhsTyy(i,jmax-2) = ( -         Tyy(i,jmax)
     &                       +  4.d0 * Tyy(i,jmax-1)
     &                       + 10.d0 * Tyy(i,jmax-2)
     &                       +  4.d0 * Tyy(i,jmax-3)
     &                       -         Tyy(i,jmax-4) ) / 16.d0 
        rhsTyy(i,jmax-1) = (           Tyy(i,jmax)
     &                       + 12.d0 * Tyy(i,jmax-1)
     &                       +  6.d0 * Tyy(i,jmax-2)
     &                       -  4.d0 * Tyy(i,jmax-3)
     &                       +         Tyy(i,jmax-4) ) / 16.d0 
        rhsTyy(i,jmax)   = (   15.d0 * Tyy(i,jmax)
     &                       +  4.d0 * Tyy(i,jmax-1)
     &                       -  6.d0 * Tyy(i,jmax-2)
     &                       +  4.d0 * Tyy(i,jmax-3)
     &                       -         Tyy(i,jmax-4) ) / 16.d0 

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_tridfy_Psi(rhsPsixx, rhsPsixy, rhsPsiyy)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer i, j
      real*8  rhsPsixx(imax,jmax), rhsPsixy(imax,jmax), 
     &        rhsPsiyy(imax,jmax)

      do i = 2, imax

        rhsPsixx(i,1) = (   15.d0 * Psixx(i,1) +  4.d0 * Psixx(i,2)
     &                    -  6.d0 * Psixx(i,3) +  4.d0 * Psixx(i,4)
     &                    -         Psixx(i,5) ) / 16.d0
        rhsPsixx(i,2) = (           Psixx(i,1) + 12.d0 * Psixx(i,2)
     &                    +  6.d0 * Psixx(i,3) -  4.d0 * Psixx(i,4)
     &                    +         Psixx(i,5) ) / 16.d0 
        rhsPsixx(i,3) = ( -         Psixx(i,1) +  4.d0 * Psixx(i,2)
     &                    + 10.d0 * Psixx(i,3) +
     &                       4.d0 * Psixx(i,4) -         Psixx(i,5) )
     &                    / 16.d0 

        rhsPsixy(i,1) = (   15.d0 * Psixy(i,1) +  4.d0 * Psixy(i,2)
     &                    -  6.d0 * Psixy(i,3) +  4.d0 * Psixy(i,4)
     &                    -         Psixy(i,5) ) / 16.d0
        rhsPsixy(i,2) = (           Psixy(i,1) + 12.d0 * Psixy(i,2)
     &                    +  6.d0 * Psixy(i,3) -  4.d0 * Psixy(i,4)
     &                    +         Psixy(i,5) ) / 16.d0 
        rhsPsixy(i,3) = ( -         Psixy(i,1) +  4.d0 * Psixy(i,2)
     &                    + 10.d0 * Psixy(i,3) +
     &                       4.d0 * Psixy(i,4) -         Psixy(i,5) )
     &                    / 16.d0 

        rhsPsiyy(i,1) = (   15.d0 * Psiyy(i,1) +  4.d0 * Psiyy(i,2)
     &                    -  6.d0 * Psiyy(i,3) +  4.d0 * Psiyy(i,4)
     &                    -         Psiyy(i,5) ) / 16.d0
        rhsPsiyy(i,2) = (           Psiyy(i,1) + 12.d0 * Psiyy(i,2)
     &                    +  6.d0 * Psiyy(i,3) -  4.d0 * Psiyy(i,4)
     &                    +         Psiyy(i,5) ) / 16.d0 
        rhsPsiyy(i,3) = ( -         Psiyy(i,1) +  4.d0 * Psiyy(i,2)
     &                    + 10.d0 * Psiyy(i,3) +
     &                       4.d0 * Psiyy(i,4) -         Psiyy(i,5) )
     &                    / 16.d0 


        do j = 4, jmax - 3
          rhsPsixx(i,j) = af * Psixx(i,j)
     &                + bf * ( Psixx(i,j+1) + Psixx(i,j-1) )
     &                + cf * ( Psixx(i,j+2) + Psixx(i,j-2) )
     &                + df * ( Psixx(i,j+3) + Psixx(i,j-3) )
          rhsPsixy(i,j) = af * Psixy(i,j)
     &                + bf * ( Psixy(i,j+1) + Psixy(i,j-1) )
     &                + cf * ( Psixy(i,j+2) + Psixy(i,j-2) )
     &                + df * ( Psixy(i,j+3) + Psixy(i,j-3) )
          rhsPsiyy(i,j) = af * Psiyy(i,j)
     &                + bf * ( Psiyy(i,j+1) + Psiyy(i,j-1) )
     &                + cf * ( Psiyy(i,j+2) + Psiyy(i,j-2) )
     &                + df * ( Psiyy(i,j+3) + Psiyy(i,j-3) )
        end do

        rhsPsixx(i,jmax-2) = ( -       Psixx(i,jmax)
     &                       +  4.d0 * Psixx(i,jmax-1)
     &                       + 10.d0 * Psixx(i,jmax-2)
     &                       +  4.d0 * Psixx(i,jmax-3)
     &                       -         Psixx(i,jmax-4) ) / 16.d0 
        rhsPsixx(i,jmax-1) = (         Psixx(i,jmax)
     &                       + 12.d0 * Psixx(i,jmax-1)
     &                       +  6.d0 * Psixx(i,jmax-2)
     &                       -  4.d0 * Psixx(i,jmax-3)
     &                       +         Psixx(i,jmax-4) ) / 16.d0 
        rhsPsixx(i,jmax)   = ( 15.d0 * Psixx(i,jmax)
     &                       +  4.d0 * Psixx(i,jmax-1)
     &                       -  6.d0 * Psixx(i,jmax-2)
     &                       +  4.d0 * Psixx(i,jmax-3)
     &                       -         Psixx(i,jmax-4) ) / 16.d0 

        rhsPsixy(i,jmax-2) = ( -       Psixy(i,jmax)
     &                       +  4.d0 * Psixy(i,jmax-1)
     &                       + 10.d0 * Psixy(i,jmax-2)
     &                       +  4.d0 * Psixy(i,jmax-3)
     &                       -         Psixy(i,jmax-4) ) / 16.d0 
        rhsPsixy(i,jmax-1) = (         Psixy(i,jmax)
     &                       + 12.d0 * Psixy(i,jmax-1)
     &                       +  6.d0 * Psixy(i,jmax-2)
     &                       -  4.d0 * Psixy(i,jmax-3)
     &                       +         Psixy(i,jmax-4) ) / 16.d0 
        rhsPsixy(i,jmax)   = ( 15.d0 * Psixy(i,jmax)
     &                       +  4.d0 * Psixy(i,jmax-1)
     &                       -  6.d0 * Psixy(i,jmax-2)
     &                       +  4.d0 * Psixy(i,jmax-3)
     &                       -         Psixy(i,jmax-4) ) / 16.d0 

        rhsPsiyy(i,jmax-2) = ( -       Psiyy(i,jmax)
     &                       +  4.d0 * Psiyy(i,jmax-1)
     &                       + 10.d0 * Psiyy(i,jmax-2)
     &                       +  4.d0 * Psiyy(i,jmax-3)
     &                       -         Psiyy(i,jmax-4) ) / 16.d0 
        rhsPsiyy(i,jmax-1) = (         Psiyy(i,jmax)
     &                       + 12.d0 * Psiyy(i,jmax-1)
     &                       +  6.d0 * Psiyy(i,jmax-2)
     &                       -  4.d0 * Psiyy(i,jmax-3)
     &                       +         Psiyy(i,jmax-4) ) / 16.d0 
        rhsPsiyy(i,jmax)   = ( 15.d0 * Psiyy(i,jmax)
     &                       +  4.d0 * Psiyy(i,jmax-1)
     &                       -  6.d0 * Psiyy(i,jmax-2)
     &                       +  4.d0 * Psiyy(i,jmax-3)
     &                       -         Psiyy(i,jmax-4) ) / 16.d0 

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Mounts the left hand side of the matrix for the second derivative in x-direction
!!
!!\f{eqnarray*}{\left[ \begin{array}{ccccccccccc} 1 &\quad &\quad &\quad &\quad &\quad &\quad &\quad &\quad \\ \quad & 1 & \quad &\quad &\quad &\quad &\quad &\quad &\quad\\ \quad & \quad & 1 &\quad &\quad & \quad &\quad &\quad &\quad \\ \quad&\quad  & \alpha_{f} & 1 & \alpha_{f} & \quad & \quad &\quad &\quad \\ \quad &\quad & \quad & \quad &\ddots & \quad & \quad &\quad &\quad \\ \quad & \quad &\quad &\quad &\alpha_{f} & 1 & \alpha_{f}  &\quad & \quad \\ \quad & \quad & \quad &\quad &\quad &\quad & 1 & \quad \\ \quad & \quad & \quad &\quad &\quad &\quad &\quad & 1 & \quad \\ \quad & \quad & \quad  &\quad &\quad &\quad &\quad &\quad & 1 \end{array} \right]\f}
      subroutine lhs_tridfx(a, b, c)

      implicit none
      include 'par.nn'
      integer i
      real*8 a(imax), b(imax), c(imax)

      a = 0.d0
      b = 1.d0
      c = 0.d0

      do i = 4, imax - 3
        a(i)    = alphaf
        c(i)    = alphaf
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Mounts the left hand side of the matrix for the second derivative in y-direction
!!
!!\f{eqnarray*}{\left[ \begin{array}{ccccccccccc} 1 &\quad &\quad &\quad &\quad &\quad &\quad &\quad &\quad \\ \quad & 1 & \quad &\quad &\quad &\quad &\quad &\quad &\quad\\ \quad & \quad & 1 &\quad &\quad & \quad &\quad &\quad &\quad \\ \quad&\quad  & \alpha_{f} & 1 & \alpha_{f} & \quad & \quad &\quad &\quad \\ \quad &\quad & \quad & \quad &\ddots & \quad & \quad &\quad &\quad \\ \quad & \quad &\quad &\quad &\alpha_{f} & 1 & \alpha_{f}  &\quad & \quad \\ \quad & \quad & \quad &\quad &\quad &\quad & 1 & \quad \\ \quad & \quad & \quad &\quad &\quad &\quad &\quad & 1 & \quad \\ \quad & \quad & \quad  &\quad &\quad &\quad &\quad &\quad & 1 \end{array} \right]\f}
      subroutine lhs_tridfy(a, b, c)

      implicit none
      include 'par.nn'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a = 0.d0
      b = 1.d0
      c = 0.d0

      do j = 4, jmax - 3
        a(j)    = alphaf
        c(j)    = alphaf
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            end of loop calculations                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
