cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_06(t)
      
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
        Txxe(i,j) = ((1 - betann)*
     -    (a1xx + a3xx*Cos(pii*x) + 
     -      a2xx*Sin(pii*x))*
     -    (b1xx + b3xx*Cos(pii*y) + 
     -      b2xx*Sin(pii*y)))/dexp(at*t) 

        Txye(i,j) = ((1 - betann)*
     -    (a1xy + a3xy*Cos(pii*x) + 
     -      a2xy*Sin(pii*x))*
     -    (b1xy + b3xy*Cos(pii*y) + 
     -      b2xy*Sin(pii*y)))/dexp(at*t)

        Tyye(i,j) = ((1 - betann)*
     -    (a1yy + a3yy*Cos(pii*x) + 
     -      a2yy*Sin(pii*x))*
     -    (b1yy + b3yy*Cos(pii*y) + 
     -      b2yy*Sin(pii*y)))/dexp(at*t)
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_06(t)

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
        tfwz(i,j) = (16*at*pii*Cos(pii*x)*Cos(pii*y))/
     -   dexp(at*t) - 
     -  (32*betann*pii**3*Cos(pii*x)*
     -     Cos(pii*y))/(dexp(at*t)*Rey) - 
     -  ((1 - betann)*
     -     (-(a3xy*pii**2*Cos(pii*x)) - 
     -       a2xy*pii**2*Sin(pii*x))*
     -     (b1xy + b3xy*Cos(pii*y) + 
     -       b2xy*Sin(pii*y)))/dexp(at*t) 
     -   + ((1 - betann)*
     -     (a2xx*pii*Cos(pii*x) - 
     -       a3xx*pii*Sin(pii*x))*
     -     (b2xx*pii*Cos(pii*y) - 
     -       b3xx*pii*Sin(pii*y)))/
     -   dexp(at*t) - 
     -  ((1 - betann)*
     -     (a2yy*pii*Cos(pii*x) - 
     -       a3yy*pii*Sin(pii*x))*
     -     (b2yy*pii*Cos(pii*y) - 
     -       b3yy*pii*Sin(pii*y)))/
     -   dexp(at*t) + 
     -  ((1 - betann)*
     -     (a1xy + a3xy*Cos(pii*x) + 
     -       a2xy*Sin(pii*x))*
     -     (-(b3xy*pii**2*Cos(pii*y)) - 
     -       b2xy*pii**2*Sin(pii*y)))/
     -   dexp(at*t)

       !!! Tensors
       tfTxx(i,j) = (16*(1 - betann)*pii*Sin(pii*x)*
     -     Sin(pii*y))/(dexp(at*t)*Rey) + 
     -  (alphaG*Rey*Wi*
     -     (((1 - betann)**2*
     -          (a1xx + a3xx*Cos(pii*x) + 
     -             a2xx*Sin(pii*x))**2*
     -          (b1xx + b3xx*Cos(pii*y) + 
     -             b2xx*Sin(pii*y))**2)/
     -        dexp(2*at*t) + 
     -       ((1 - betann)**2*
     -          (a1xy + a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))**2*
     -          (b1xy + b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))**2)/
     -        dexp(2*at*t)))/(1 - betann) 
     -   + Wi*(-((at*(1 - betann)*
     -          (a1xx + a3xx*Cos(pii*x) + 
     -            a2xx*Sin(pii*x))*
     -          (b1xx + b3xx*Cos(pii*y) + 
     -            b2xx*Sin(pii*y)))/
     -        dexp(at*t)) + 
     -     (16*(1 - betann)*pii*
     -        Sin(pii*x)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (16*(1 - betann)*pii*xi*
     -        Sin(pii*x)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y)))/
     -      dexp(2*at*t) + 
     -     (8*(1 - betann)*Cos(pii*x)*
     -        (a2xx*pii*Cos(pii*x) - 
     -          a3xx*pii*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (16*(1 - betann)*pii*
     -        Cos(pii*x)*Cos(pii*y)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (8*(1 - betann)*Cos(pii*y)*
     -        Sin(pii*x)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b2xx*pii*Cos(pii*y) - 
     -          b3xx*pii*Sin(pii*y)))/
     -      dexp(2*at*t)) + 
     -  ((1 - betann)*
     -     (a1xx + a3xx*Cos(pii*x) + 
     -       a2xx*Sin(pii*x))*
     -     (b1xx + b3xx*Cos(pii*y) + 
     -       b2xx*Sin(pii*y))*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + 
     -               a3xx*Cos(pii*x) + 
     -               a2xx*Sin(pii*x))*
     -               (b1xx + 
     -               b3xx*Cos(pii*y) + 
     -               b2xx*Sin(pii*y)))/
     -             dexp(at*t) + 
     -            ((1 - betann)*
     -               (a1yy + 
     -               a3yy*Cos(pii*x) + 
     -               a2yy*Sin(pii*x))*
     -               (b1yy + 
     -               b3yy*Cos(pii*y) + 
     -               b2yy*Sin(pii*y)))/
     -             dexp(at*t)))/
     -        (1 - betann)))/dexp(at*t)

       tfTxy(i,j) = (alphaG*Rey*Wi*
     -     (((1 - betann)**2*
     -          (a1xx + a3xx*Cos(pii*x) + 
     -            a2xx*Sin(pii*x))*
     -          (a1xy + a3xy*Cos(pii*x) + 
     -            a2xy*Sin(pii*x))*
     -          (b1xx + b3xx*Cos(pii*y) + 
     -            b2xx*Sin(pii*y))*
     -          (b1xy + b3xy*Cos(pii*y) + 
     -            b2xy*Sin(pii*y)))/
     -        dexp(2*at*t) + 
     -       ((1 - betann)**2*
     -          (a1xy + a3xy*Cos(pii*x) + 
     -            a2xy*Sin(pii*x))*
     -          (a1yy + a3yy*Cos(pii*x) + 
     -            a2yy*Sin(pii*x))*
     -          (b1xy + b3xy*Cos(pii*y) + 
     -            b2xy*Sin(pii*y))*
     -          (b1yy + b3yy*Cos(pii*y) + 
     -            b2yy*Sin(pii*y)))/
     -        dexp(2*at*t)))/(1 - betann) 
     -   + Wi*((8*(1 - betann)*pii*
     -        Cos(pii*x)*Cos(pii*y)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (at*(1 - betann)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y)))/
     -      dexp(at*t) + 
     -     (8*(1 - betann)*Cos(pii*x)*
     -        (a2xy*pii*Cos(pii*x) - 
     -          a3xy*pii*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (8*(1 - betann)*pii*Cos(pii*x)*
     -        Cos(pii*y)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (8*(1 - betann)*Cos(pii*y)*
     -        Sin(pii*x)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b2xy*pii*Cos(pii*y) - 
     -          b3xy*pii*Sin(pii*y)))/
     -      dexp(2*at*t)) + 
     -  ((1 - betann)*
     -     (a1xy + a3xy*Cos(pii*x) + 
     -       a2xy*Sin(pii*x))*
     -     (b1xy + b3xy*Cos(pii*y) + 
     -       b2xy*Sin(pii*y))*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + 
     -               a3xx*Cos(pii*x) + 
     -               a2xx*Sin(pii*x))*
     -               (b1xx + 
     -               b3xx*Cos(pii*y) + 
     -               b2xx*Sin(pii*y)))/
     -             dexp(at*t) + 
     -            ((1 - betann)*
     -               (a1yy + 
     -               a3yy*Cos(pii*x) + 
     -               a2yy*Sin(pii*x))*
     -               (b1yy + 
     -               b3yy*Cos(pii*y) + 
     -               b2yy*Sin(pii*y)))/
     -             dexp(at*t)))/
     -        (1 - betann)))/dexp(at*t)

       tfTyy(i,j) =  (-16*(1 - betann)*pii*Sin(pii*x)*
     -     Sin(pii*y))/(dexp(at*t)*Rey) + 
     -  (alphaG*Rey*Wi*
     -     (((1 - betann)**2*
     -          (a1xy + a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))**2*
     -          (b1xy + b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))**2)/
     -        dexp(2*at*t) + 
     -       ((1 - betann)**2*
     -          (a1yy + a3yy*Cos(pii*x) + 
     -             a2yy*Sin(pii*x))**2*
     -          (b1yy + b3yy*Cos(pii*y) + 
     -             b2yy*Sin(pii*y))**2)/
     -        dexp(2*at*t)))/(1 - betann) 
     -   + Wi*((16*(1 - betann)*pii*
     -        Cos(pii*x)*Cos(pii*y)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (at*(1 - betann)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*pii*
     -        Sin(pii*x)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y)))/
     -      dexp(2*at*t) + 
     -     (16*(1 - betann)*pii*xi*
     -        Sin(pii*x)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y)))/
     -      dexp(2*at*t) + 
     -     (8*(1 - betann)*Cos(pii*x)*
     -        (a2yy*pii*Cos(pii*x) - 
     -          a3yy*pii*Sin(pii*x))*
     -        Sin(pii*y)*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y)))/
     -      dexp(2*at*t) - 
     -     (8*(1 - betann)*Cos(pii*y)*
     -        Sin(pii*x)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b2yy*pii*Cos(pii*y) - 
     -          b3yy*pii*Sin(pii*y)))/
     -      dexp(2*at*t)) + 
     -  ((1 - betann)*
     -     (a1yy + a3yy*Cos(pii*x) + 
     -       a2yy*Sin(pii*x))*
     -     (b1yy + b3yy*Cos(pii*y) + 
     -       b2yy*Sin(pii*y))*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + 
     -               a3xx*Cos(pii*x) + 
     -               a2xx*Sin(pii*x))*
     -               (b1xx + 
     -               b3xx*Cos(pii*y) + 
     -               b2xx*Sin(pii*y)))/
     -             dexp(at*t) + 
     -            ((1 - betann)*
     -               (a1yy + 
     -               a3yy*Cos(pii*x) + 
     -               a2yy*Sin(pii*x))*
     -               (b1yy + 
     -               b3yy*Cos(pii*y) + 
     -               b2yy*Sin(pii*y)))/
     -             dexp(at*t)))/
     -        (1 - betann)))/dexp(at*t)

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_06(t)

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
        Txx(i,j) = ((1 - betann)*
     -    (a1xx + a3xx*Cos(pii*x) + 
     -      a2xx*Sin(pii*x))*
     -    (b1xx + b3xx*Cos(pii*y) + 
     -      b2xx*Sin(pii*y)))/dexp(at*t) 

        Txy(i,j) = ((1 - betann)*
     -    (a1xy + a3xy*Cos(pii*x) + 
     -      a2xy*Sin(pii*x))*
     -    (b1xy + b3xy*Cos(pii*y) + 
     -      b2xy*Sin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*
     -    (a1yy + a3yy*Cos(pii*x) + 
     -      a2yy*Sin(pii*x))*
     -    (b1yy + b3yy*Cos(pii*y) + 
     -      b2yy*Sin(pii*y)))/dexp(at*t) 
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
        Txx(i,j) = ((1 - betann)*
     -    (a1xx + a3xx*Cos(pii*x) + 
     -      a2xx*Sin(pii*x))*
     -    (b1xx + b3xx*Cos(pii*y) + 
     -      b2xx*Sin(pii*y)))/dexp(at*t) 

        Txy(i,j) = ((1 - betann)*
     -    (a1xy + a3xy*Cos(pii*x) + 
     -      a2xy*Sin(pii*x))*
     -    (b1xy + b3xy*Cos(pii*y) + 
     -      b2xy*Sin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*
     -    (a1yy + a3yy*Cos(pii*x) + 
     -      a2yy*Sin(pii*x))*
     -    (b1yy + b3yy*Cos(pii*y) + 
     -      b2yy*Sin(pii*y)))/dexp(at*t) 
         end do
      end do
      return
      end
