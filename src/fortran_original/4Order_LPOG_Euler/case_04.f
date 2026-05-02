cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_04(t)
      
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
         uxe(i,j) = (8*Cos(pii*x)*Sin(pii*y))/dexp(at*t)

         uye(i,j) = (-8*Cos(pii*y)*Sin(pii*x))/dexp(at*t)

         wze(i,j) = (16*pii*Cos(pii*x)*Cos(pii*y))/dexp(at*t)

        psie(i,j) = (-8*Cos(pii*x)*Cos(pii*y))/
     -  (dexp(at*t)*pii)

         !!! Tensors
        Txxe(i,j) =  -(((-1 + betann)*c1)/dexp(at*t))

        Txye(i,j) = -(((-1 + betann)*c2)/dexp(at*t))

        Tyye(i,j) = -(((-1 + betann)*c3)/dexp(at*t))
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_04(t)

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
        tfwz(i,j) = (16*pii*(-2*betann*pii**2 + 
     -      at*Rey)*Cos(pii*x)*Cos(pii*y))
     -   /(dexp(at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) = ((-1 + betann)*
     -    (Rey*(-((alphaG*
     -              (c1**2 + c2**2) + 
     -              c1*(c1 + c3)*epsylon)*
     -            Rey*Wi) + 
     -         c1*dexp(at*t)*(-1 + at*Wi)) 
     -       + 16*c2*pii*Rey*Wi*
     -       Cos(pii*x)*Cos(pii*y) - 
     -      16*pii*
     -       (dexp(at*t) - 
     -         c1*Rey*Wi*(-1 + xi))*
     -       Sin(pii*x)*Sin(pii*y)))/
     -  (dexp(2*at*t)*Rey)

       tfTxy(i,j) = ((-1 + betann)*
     -    (-(c2*(c1 + c3)*
     -         (alphaG + epsylon)*Rey*Wi) 
     -       + c2*dexp(at*t)*
     -       (-1 + at*Wi) + 
     -      8*(-c1 + c3)*pii*Wi*
     -       Cos(pii*x)*Cos(pii*y)))/
     -  dexp(2*at*t)

       tfTyy(i,j) =  ((-1 + betann)*
     -    (Rey*(-((alphaG*
     -              (c2**2 + c3**2) + 
     -              c3*(c1 + c3)*epsylon)*
     -            Rey*Wi) + 
     -         c3*dexp(at*t)*(-1 + at*Wi)) 
     -       - 16*c2*pii*Rey*Wi*
     -       Cos(pii*x)*Cos(pii*y) + 
     -      16*pii*
     -       (dexp(at*t) - 
     -         c3*Rey*Wi*(-1 + xi))*
     -       Sin(pii*x)*Sin(pii*y)))/
     -  (dexp(2*at*t)*Rey)

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_04(t)

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
         ux(i,j) = (8*Cos(pii*x)*Sin(pii*y))/dexp(at*t)

         uy(i,j) = (-8*Cos(pii*y)*Sin(pii*x))/dexp(at*t)

         wz(i,j) = (16*pii*Cos(pii*x)*Cos(pii*y))/dexp(at*t)

        psi(i,j) = (-8*Cos(pii*x)*Cos(pii*y))/
     -  (dexp(at*t)*pii)

         !!! Tensors
        Txx(i,j) =  -(((-1 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1 + betann)*c3)/dexp(at*t)) 
         end do
      end do
      !!! Boundary: top and botton
      do j = 1, jmax, jmax - 1
         y = y0 + dble(j-1)*dy
         do i = 1, imax
            x     = x0 + dble(i-1)*dx
         !!! Velocitys 
         ux(i,j) = (8*Cos(pii*x)*Sin(pii*y))/dexp(at*t)

         uy(i,j) = (-8*Cos(pii*y)*Sin(pii*x))/dexp(at*t)

         wz(i,j) = (16*pii*Cos(pii*x)*Cos(pii*y))/dexp(at*t)

        psi(i,j) = (-8*Cos(pii*x)*Cos(pii*y))/
     -  (dexp(at*t)*pii)

         !!! Tensors
        Txx(i,j) =  -(((-1 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1 + betann)*c3)/dexp(at*t)) 
         end do
      end do
      return
      end
