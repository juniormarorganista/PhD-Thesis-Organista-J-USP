cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_12(t)
      
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
        Txxe(i,j) =  ((1.0d0 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*dsin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*dsin(pii*y)))/dexp(at*t)

        Txye(i,j) = ((1.0d0 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*dsin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*dsin(pii*y)))/dexp(at*t)

        Tyye(i,j) = ((1.0d0 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*dsin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*dsin(pii*y)))/dexp(at*t)

         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_12(t)

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
        tfwz(i,j) = (-(Rey*y) + 2.0d0*Rey*x*y + 
     & Rey*y*dcos(2.0d0*pii*x) - 
     & 2.0d0*Rey*x*y*dcos(2.0d0*pii*x) + 
     & 2.0d0*pii**2.0d0*Rey*x*y*dcos(2.0d0*pii*x) - 
     & 6*pii**2.0d0*Rey*x**2.0d0*y*dcos(2.0d0*pii*x) + 
     & 4*pii**2.0d0*Rey*x**3*y*dcos(2.0d0*pii*x) + 
     & pii*Rey*y*dsin(2.0d0*pii*x) - 
     & 2.0d0*pii*Rey*x*y*dsin(2.0d0*pii*x) + 
     & 2.0d0*pii*Rey*x**2.0d0*y*dsin(2.0d0*pii*x) - 
     & 2.0d0*pii**3*Rey*x**2.0d0*y*dsin(2.0d0*pii*x) + 
     & 4*pii**3*Rey*x**3*y*dsin(2.0d0*pii*x) - 
     & 2.0d0*pii**3*Rey*x**4*y*dsin(2.0d0*pii*x) + 
     & Rey*y*dcos(2.0d0*pii*y)*
     & (-((-1.0d0 + 2.0d0*x)*
     & (1.0d0 + 2.0d0*pii**2.0d0*(-1.0d0 + x)*x - dcos(2.0d0*pii*x))) + 
     & pii*(-1.0d0 - 2.0d0*(-1.0d0 + x)*x)*dsin(2.0d0*pii*x)) + 
     & 2.0d0*dexp(at*t)*pii*dcos(pii*x)*
     & (-(a3xy*b1xy*(-1.0d0 + betann)*pii*Rey) + 
     & pii*((a2xx*b2xx - a2yy*b2yy)*Rey + 
     & betann*(-8 - a2xx*b2xx*Rey + 
     & a2yy*b2yy*Rey + 16*x))*dcos(pii*y) + 
     & ((a2xx*b3xx - a2yy*b3yy)*(-1.0d0 + betann)*pii*
     & Rey - 2.0d0*(4*betann*pii**2.0d0 - at*Rey)*
     &        (-1.0d0 + 2.0d0*x)*y)*dsin(pii*y)) + 
     & dexp(at*t)*(-2.0d0*a2xy*b1xy*(-1.0d0 + betann)*pii**2*Rey*
     &     dsin(pii*x) + 
     &    2.0d0*pii*dcos(pii*y)*
     &     (a1xy*b3xy*(-1.0d0 + betann)*pii*Rey + 
     &       (Rey*(-(a3xx*b2xx*pii) + a3yy*b2yy*pii + 
     &             2.0d0*at*(-1.0d0 + x)*x) + 
     &          betann*
     &           (8 + a3xx*b2xx*pii*Rey - 
     &             pii*(a3yy*b2yy*Rey + 8*pii*(-1.0d0 + x)*x)
     &             ))*dsin(pii*x)) + 
     &    2.0d0*(a1xy*b2xy*(-1.0d0 + betann)*pii**2*Rey + 
     &       (-((a3xx*b3xx - a3yy*b3yy)*(-1.0d0 + betann)*
     &             pii**2.0d0*Rey) + 
     &          2.0d0*(at*Rey + 
     &             pii**2*
     &              (-(at*Rey*(-1.0d0 + x)*x) + 
     &                2*betann*(-4 + pii**2*(-1.0d0 + x)*x)))
     &            *y)*dsin(pii*x))*dsin(pii*y)) + 
     & pii*Rey*((-1.0d0 + 2*x)*
     &     (y**2 + (-1.0d0 + x)*x*(-1.0d0 + 2*pii**2*y**2)) + 
     &    (-1.0d0 + 2*x)*((-1.0d0 + x)*x - y**2)*dcos(2*pii*x) - 
     &    pii*((-1.0d0 + x)**2*x**2 + 
     &       (-1.0d0 - 2*(-1.0d0 + x)*x)*y**2)*dsin(2*pii*x))*
     &  dsin(2*pii*y))/(2.*dexp(2*at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) =  ((2*(-1.0d0 + betann)*dexp(at*t)*
     &    (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &      (-1.0d0 + 2*x)*dsin(pii*x))*
     &    (pii*y*dcos(pii*y) + dsin(pii*y)))/Rey + 
     & alphaG*(1.0d0 - betann)*Rey*Wi*
     &  ((a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))**2*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y))**2
     &     + (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))**
     &      2*(b1xy + b3xy*dcos(pii*y) + 
     &        b2xy*dsin(pii*y))**2) + 
     & (1.0d0 - betann)*dexp(at*t)*
     &  (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &  (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y))*
     &  (1.0d0 + (epsylon*Rey*Wi*
     &       (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &       (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)))/
     &     dexp(at*t) + 
     &    (epsylon*Rey*Wi*
     &       (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &       (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)))/
     &     dexp(at*t)) + 
     & Wi*(at*(-1.0d0 + betann)*dexp(at*t)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*y*dcos(pii*y)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*pii*(-1.0d0 + x)*x*dcos(pii*x)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*(-1.0d0 + x)*dsin(pii*x)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*x*dsin(pii*x)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &     (-(a2xx*dcos(pii*x)) + a3xx*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    2*(-1.0d0 + betann)*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*y*
     &     (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b2xx*dcos(pii*y) - b3xx*dsin(pii*y)) - 
     &    2*(-1.0d0 + betann)*pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))*
     &     (-2*dcos(pii*y) + pii*y*dsin(pii*y)) + 
     &    (1.0d0 - betann)*xi*
     &     (2*(a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &        (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &          (-1.0d0 + 2*x)*dsin(pii*x))*
     &        (pii*y*dcos(pii*y) + dsin(pii*y))*
     &        (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y))
     &        + y*(a1xy + a3xy*dcos(pii*x) + 
     &          a2xy*dsin(pii*x))*
     &        (2*pii*(1.0d0 - 2*x)*dcos(pii*x) + 
     &          (-2 + pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
     &        dsin(pii*y)*
     &        (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))
     &        - pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &        (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &        (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))*
     &        (-2*dcos(pii*y) + pii*y*dsin(pii*y)))))/
     &  dexp(2*at*t)
     
       tfTxy(i,j) = ((2*(-1.0d0 + betann)*dexp(at*t)*
     &    (pii*(-1.0d0 + x)*x*dcos(pii*y)*dsin(pii*x) - 
     &      y*(pii*(-1.0d0 + 2*x)*dcos(pii*x) + dsin(pii*x))*
     &       dsin(pii*y)))/Rey + 
     & alphaG*(1.0d0 - betann)*Rey*Wi*
     &  (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &  (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))*
     &  ((a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y))) + 
     & (1.0d0 - betann)*dexp(at*t)*
     &  (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &  (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))*
     &  (1.0d0 + (epsylon*Rey*Wi*
     &       (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &       (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)))/
     &     dexp(at*t) + 
     &    (epsylon*Rey*Wi*
     &       (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &       (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)))/
     &     dexp(at*t)) + 
     & (-1.0d0 + betann)*Wi*
     &  (y*(a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &     (2*pii*(1.0d0 - 2*x)*dcos(pii*x) + 
     &       (-2 + pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
     &     dsin(pii*y)*
     &     (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)) + 
     &    at*dexp(at*t)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) + 
     &    pii*y*dcos(pii*y)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) + 
     &    (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) - 
     &    pii*(-1.0d0 + x)*x*dcos(pii*x)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) - 
     &    (-1.0d0 + x)*dsin(pii*x)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) - 
     &    x*dsin(pii*x)*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) + 
     &    pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &     (-(a2xy*dcos(pii*x)) + a3xy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) + 
     &    pii*y*(a1xy + a3xy*dcos(pii*x) + 
     &       a2xy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b2xy*dcos(pii*y) - b3xy*dsin(pii*y)) + 
     &    pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y))*
     &     (2*dcos(pii*y) - pii*y*dsin(pii*y)) + 
     &    xi*(-(pii*(-1.0d0 + x)*x*dcos(pii*y)*dsin(pii*x)) + 
     &       y*(pii*(-1.0d0 + 2*x)*dcos(pii*x) + dsin(pii*x))*
     &        dsin(pii*y))*
     &     (a1xx*b1xx + a1yy*b1yy + 
     &       dcos(pii*y)*
     &        (a1xx*b3xx + a1yy*b3yy + 
     &          (a2xx*b3xx + a2yy*b3yy)*dsin(pii*x)) + 
     &       (a1xx*b2xx + a1yy*b2yy)*dsin(pii*y) + 
     &       dsin(pii*x)*
     &        (a2xx*b1xx + a2yy*b1yy + 
     &          (a2xx*b2xx + a2yy*b2yy)*dsin(pii*y)) + 
     &       dcos(pii*x)*
     &        (a3xx*b1xx + a3yy*b1yy + 
     &          (a3xx*b3xx + a3yy*b3yy)*dcos(pii*y) + 
     &          (a3xx*b2xx + a3yy*b2yy)*dsin(pii*y)))))/
     &  dexp(2*at*t)
     
       tfTyy(i,j) = ((-2*(-1.0d0 + betann)*dexp(at*t)*
     &    (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &      (-1.0d0 + 2*x)*dsin(pii*x))*
     &    (pii*y*dcos(pii*y) + dsin(pii*y)))/Rey + 
     & (1.0d0 - betann)*dexp(at*t)*
     &  (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &  (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y))*
     &  (1.0d0 + (epsylon*Rey*Wi*
     &       (a1xx + a3xx*dcos(pii*x) + a2xx*dsin(pii*x))*
     &       (b1xx + b3xx*dcos(pii*y) + b2xx*dsin(pii*y)))/
     &     dexp(at*t) + 
     &    (epsylon*Rey*Wi*
     &       (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &       (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)))/
     &     dexp(at*t)) + 
     & alphaG*(1.0d0 - betann)*Rey*Wi*
     &  ((a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))**2*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))**2 
     &     + (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))**
     &      2*(b1yy + b3yy*dcos(pii*y) + 
     &        b2yy*dsin(pii*y))**2) + 
     & Wi*(2*(-1.0d0 + betann)*y*
     &     (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &     (2*pii*(1.0d0 - 2*x)*dcos(pii*x) + 
     &       (-2 + pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
     &     dsin(pii*y)*
     &     (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) + 
     &    at*(-1.0d0 + betann)*dexp(at*t)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*y*dcos(pii*y)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*pii*(-1.0d0 + x)*x*dcos(pii*x)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*(-1.0d0 + x)*dsin(pii*x)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) - 
     &    (-1.0d0 + betann)*x*dsin(pii*x)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &     (-(a2yy*dcos(pii*x)) + a3yy*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) - 
     &    2*(-1.0d0 + betann)*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*
     &     (pii*y*dcos(pii*y) + dsin(pii*y))*
     &     (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) + 
     &    (-1.0d0 + betann)*pii*y*
     &     (a1yy + a3yy*dcos(pii*x) + a2yy*dsin(pii*x))*
     &     (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &       (-1.0d0 + 2*x)*dsin(pii*x))*dsin(pii*y)*
     &     (b2yy*dcos(pii*y) - b3yy*dsin(pii*y)) + 
     &    (1.0d0 - betann)*xi*
     &     (y*(a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &        (2*pii*(1.0d0 - 2*x)*dcos(pii*x) + 
     &          (-2 + pii**2*(-1.0d0 + x)*x)*dsin(pii*x))*
     &        dsin(pii*y)*
     &        (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y)) 
     &        - 2*(a1yy + a3yy*dcos(pii*x) + 
     &          a2yy*dsin(pii*x))*
     &        (pii*(-1.0d0 + x)*x*dcos(pii*x) + 
     &          (-1.0d0 + 2*x)*dsin(pii*x))*
     &        (pii*y*dcos(pii*y) + dsin(pii*y))*
     &        (b1yy + b3yy*dcos(pii*y) + b2yy*dsin(pii*y)) 
     &        - pii*(-1.0d0 + x)*x*dsin(pii*x)*
     &        (a1xy + a3xy*dcos(pii*x) + a2xy*dsin(pii*x))*
     &        (b1xy + b3xy*dcos(pii*y) + b2xy*dsin(pii*y))*
     &        (-2*dcos(pii*y) + pii*y*dsin(pii*y)))))/
     &  dexp(2*at*t)
     

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_12(t)

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
        Txx(i,j) =  ((1.0d0 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*dsin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*dsin(pii*y)))/dexp(at*t)

        Txy(i,j) = ((1.0d0 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*dsin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*dsin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1.0d0 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*dsin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*dsin(pii*y)))/dexp(at*t) 

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
        Txx(i,j) =  ((1.0d0 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*dsin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*dsin(pii*y)))/dexp(at*t)

        Txy(i,j) = ((1.0d0 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*dsin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*dsin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1.0d0 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*dsin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*dsin(pii*y)))/dexp(at*t) 

         end do
      end do

      call boundarie_vorticity_cavity_flow 

      return
      end
