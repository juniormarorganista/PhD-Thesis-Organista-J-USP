cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_05(t)
      
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
     &  (dexp(at*t)*pii)

         !!! Tensors
        Txxe(i,j) =  ((1 - betann)*
     &    (a1xx + a2xx*x + a3xx*x**2)*
     &    (b1xx + b2xx*y + b3xx*y**2))/
     &  dexp(at*t)

        Txye(i,j) = ((1 - betann)*
     &    (a1xy + a2xy*x + a3xy*x**2)*
     &    (b1xy + b2xy*y + b3xy*y**2))/
     &  dexp(at*t)

        Tyye(i,j) = ((1 - betann)*
     &    (a1yy + a2yy*x + a3yy*x**2)*
     &    (b1yy + b2yy*y + b3yy*y**2))/
     &  dexp(at*t)
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_05(t)

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
        tfwz(i,j) = (2*b3xy*(1 - betann)*
     &     (a1xy + a2xy*x + a3xy*x**2))/
     &   dexp(at*t) + 
     &  ((1 - betann)*(a2xx + 2*a3xx*x)*
     &     (b2xx + 2*b3xx*y))/dexp(at*t) - 
     &  ((1 - betann)*(a2yy + 2*a3yy*x)*
     &     (b2yy + 2*b3yy*y))/dexp(at*t) - 
     &  (2*a3xy*(1 - betann)*
     &     (b1xy + b2xy*y + b3xy*y**2))/
     &   dexp(at*t) + 
     &  (16*at*pii*Cos(pii*x)*Cos(pii*y))/
     &   dexp(at*t) - 
     &  (32*betann*pii**3*Cos(pii*x)*
     &     Cos(pii*y))/(dexp(at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) = (alphaG*Rey*Wi*
     &     (((1 - betann)**2*
     &          (a1xx + a2xx*x + 
     &             a3xx*x**2)**2*
     &          (b1xx + b2xx*y + 
     &             b3xx*y**2)**2)/
     &        dexp(2*at*t) + 
     &       ((1 - betann)**2*
     &          (a1xy + a2xy*x + 
     &             a3xy*x**2)**2*
     &          (b1xy + b2xy*y + 
     &             b3xy*y**2)**2)/
     &        dexp(2*at*t)))/(1 - betann) 
     &   + ((1 - betann)*
     &     (a1xx + a2xx*x + a3xx*x**2)*
     &     (b1xx + b2xx*y + b3xx*y**2)*
     &     (1 + (epsylon*Rey*Wi*
     &          (((1 - betann)*
     &               (a1xx + a2xx*x + 
     &               a3xx*x**2)*
     &               (b1xx + b2xx*y + 
     &               b3xx*y**2))/dexp(at*t)
     &              + ((1 - betann)*
     &               (a1yy + a2yy*x + 
     &               a3yy*x**2)*
     &               (b1yy + b2yy*y + 
     &               b3yy*y**2))/dexp(at*t)
     &            ))/(1 - betann)))/
     &   dexp(at*t) + 
     &  (16*(1 - betann)*pii*Sin(pii*x)*
     &     Sin(pii*y))/(dexp(at*t)*Rey) + 
     &  Wi*(-((at*(1 - betann)*
     &          (a1xx + a2xx*x + 
     &            a3xx*x**2)*
     &          (b1xx + b2xx*y + 
     &            b3xx*y**2))/dexp(at*t)) 
     &      - (16*(1 - betann)*pii*
     &        (a1xy + a2xy*x + a3xy*x**2)*
     &        (b1xy + b2xy*y + b3xy*y**2)*
     &        Cos(pii*x)*Cos(pii*y))/
     &      dexp(2*at*t) - 
     &     (8*(1 - betann)*
     &        (a1xx + a2xx*x + a3xx*x**2)*
     &        (b2xx + 2*b3xx*y)*
     &        Cos(pii*y)*Sin(pii*x))/
     &      dexp(2*at*t) + 
     &     (8*(1 - betann)*
     &        (a2xx + 2*a3xx*x)*
     &        (b1xx + b2xx*y + b3xx*y**2)*
     &        Cos(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t) + 
     &     (16*(1 - betann)*pii*
     &        (a1xx + a2xx*x + a3xx*x**2)*
     &        (b1xx + b2xx*y + b3xx*y**2)*
     &        Sin(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t) - 
     &     (16*(1 - betann)*pii*
     &        (a1xx + a2xx*x + a3xx*x**2)*
     &        xi*
     &        (b1xx + b2xx*y + b3xx*y**2)*
     &        Sin(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t))

       tfTxy(i,j) = (alphaG*Rey*Wi*
     &     (((1 - betann)**2*
     &          (a1xx + a2xx*x + 
     &            a3xx*x**2)*
     &          (a1xy + a2xy*x + 
     &            a3xy*x**2)*
     &          (b1xx + b2xx*y + 
     &            b3xx*y**2)*
     &          (b1xy + b2xy*y + 
     &            b3xy*y**2))/dexp(2*at*t) 
     &        + ((1 - betann)**2*
     &          (a1xy + a2xy*x + 
     &            a3xy*x**2)*
     &          (a1yy + a2yy*x + 
     &            a3yy*x**2)*
     &          (b1xy + b2xy*y + 
     &            b3xy*y**2)*
     &          (b1yy + b2yy*y + 
     &            b3yy*y**2))/dexp(2*at*t))
     &     )/(1 - betann) + 
     &  ((1 - betann)*
     &     (a1xy + a2xy*x + a3xy*x**2)*
     &     (b1xy + b2xy*y + b3xy*y**2)*
     &     (1 + (epsylon*Rey*Wi*
     &          (((1 - betann)*
     &               (a1xx + a2xx*x + 
     &               a3xx*x**2)*
     &               (b1xx + b2xx*y + 
     &               b3xx*y**2))/dexp(at*t)
     &              + ((1 - betann)*
     &               (a1yy + a2yy*x + 
     &               a3yy*x**2)*
     &               (b1yy + b2yy*y + 
     &               b3yy*y**2))/dexp(at*t)
     &            ))/(1 - betann)))/
     &   dexp(at*t) + 
     &  Wi*(-((at*(1 - betann)*
     &          (a1xy + a2xy*x + 
     &            a3xy*x**2)*
     &          (b1xy + b2xy*y + 
     &            b3xy*y**2))/dexp(at*t)) 
     &      + (8*(1 - betann)*pii*
     &        (a1xx + a2xx*x + a3xx*x**2)*
     &        (b1xx + b2xx*y + b3xx*y**2)*
     &        Cos(pii*x)*Cos(pii*y))/
     &      dexp(2*at*t) - 
     &     (8*(1 - betann)*pii*
     &        (a1yy + a2yy*x + a3yy*x**2)*
     &        (b1yy + b2yy*y + b3yy*y**2)*
     &        Cos(pii*x)*Cos(pii*y))/
     &      dexp(2*at*t) - 
     &     (8*(1 - betann)*
     &        (a1xy + a2xy*x + a3xy*x**2)*
     &        (b2xy + 2*b3xy*y)*
     &        Cos(pii*y)*Sin(pii*x))/
     &      dexp(2*at*t) + 
     &     (8*(1 - betann)*
     &        (a2xy + 2*a3xy*x)*
     &        (b1xy + b2xy*y + b3xy*y**2)*
     &        Cos(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t))

       tfTyy(i,j) =  (alphaG*Rey*Wi*(((1 - betann)**2*
     & (a1xy + a2xy*x + a3xy*x**2)**2*(b1xy + b2xy*y + b3xy*y**2)**2)/
     & dexp(2*at*t) + ((1 - betann)**2*(a1yy + a2yy*x + a3yy*x**2)**2*
     & (b1yy + b2yy*y + b3yy*y**2)**2)/dexp(2*at*t)))/(1 - betann) 
     & + ((1 - betann)*(a1yy + a2yy*x + a3yy*x**2)*
     & (b1yy + b2yy*y + b3yy*y**2)*(1 + (epsylon*Rey*Wi*(((1 - betann)*
     & (a1xx + a2xx*x + a3xx*x**2)*(b1xx + b2xx*y + 
     &               b3xx*y**2))/dexp(at*t)
     &              + ((1 - betann)*
     &               (a1yy + a2yy*x + 
     &               a3yy*x**2)*
     &               (b1yy + b2yy*y + 
     &               b3yy*y**2))/dexp(at*t)
     &            ))/(1 - betann)))/
     &   dexp(at*t) - 
     &  (16*(1 - betann)*pii*Sin(pii*x)*
     &     Sin(pii*y))/(dexp(at*t)*Rey) + 
     &  Wi*(-((at*(1 - betann)*
     &          (a1yy + a2yy*x + 
     &            a3yy*x**2)*
     &          (b1yy + b2yy*y + 
     &            b3yy*y**2))/dexp(at*t)) 
     &      + (16*(1 - betann)*pii*
     &        (a1xy + a2xy*x + a3xy*x**2)*
     &        (b1xy + b2xy*y + b3xy*y**2)*
     &        Cos(pii*x)*Cos(pii*y))/
     &      dexp(2*at*t) - 
     &     (8*(1 - betann)*
     &        (a1yy + a2yy*x + a3yy*x**2)*
     &        (b2yy + 2*b3yy*y)*
     &        Cos(pii*y)*Sin(pii*x))/
     &      dexp(2*at*t) + 
     &     (8*(1 - betann)*
     &        (a2yy + 2*a3yy*x)*
     &        (b1yy + b2yy*y + b3yy*y**2)*
     &        Cos(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t) - 
     &     (16*(1 - betann)*pii*
     &        (a1yy + a2yy*x + a3yy*x**2)*
     &        (b1yy + b2yy*y + b3yy*y**2)*
     &        Sin(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t) + 
     &     (16*(1 - betann)*pii*
     &        (a1yy + a2yy*x + a3yy*x**2)*
     &        xi*
     &        (b1yy + b2yy*y + b3yy*y**2)*
     &        Sin(pii*x)*Sin(pii*y))/
     &      dexp(2*at*t))

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_05(t)

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
     &  (dexp(at*t)*pii)

         !!! Tensors
        Txx(i,j) =  ((1 - betann)*(a1xx + a2xx*x + a3xx*x**2)*
     &    (b1xx + b2xx*y + b3xx*y**2))/dexp(at*t)

        Txy(i,j) = ((1 - betann)*(a1xy + a2xy*x + a3xy*x**2)*
     &    (b1xy + b2xy*y + b3xy*y**2))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*(a1yy + a2yy*x + a3yy*x**2)*
     &    (b1yy + b2yy*y + b3yy*y**2))/dexp(at*t) 
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
     &  (dexp(at*t)*pii)

         !!! Tensors
        Txx(i,j) =  ((1 - betann)*(a1xx + a2xx*x + a3xx*x**2)*
     &    (b1xx + b2xx*y + b3xx*y**2))/dexp(at*t)

        Txy(i,j) = ((1 - betann)*(a1xy + a2xy*x + a3xy*x**2)*
     &    (b1xy + b2xy*y + b3xy*y**2))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*(a1yy + a2yy*x + a3yy*x**2)*
     &    (b1yy + b2yy*y + b3yy*y**2))/dexp(at*t) 
         end do
      end do
      return
      end
