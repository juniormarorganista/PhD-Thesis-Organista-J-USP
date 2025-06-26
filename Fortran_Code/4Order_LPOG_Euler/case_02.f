cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_02(t)
      
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
         uxe(i,j) = (16*(-1 + x)**2*x**2*y*
     -    (-1 + 2*y**2))/dexp(at*t)

         uye(i,j) = (-16*(-1 + x)*x*(-1 + 2*x)*y**2*
     -    (-1 + y**2))/dexp(at*t)

         wze(i,j) = (16*(-((-1 + x)**2*x**2) + 
     -      (-1 + 6*(x - 2*x**3 + x**4))*
     -       y**2 + 
     -      (1 + 6*(-1 + x)*x)*y**4))/
     -  dexp(at*t)

        psie(i,j) = (8*(-1 + x)**2*x**2*y**2*
     -    (-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txxe(i,j) =  -(((-1 + betann)*
     -      (a1xx + x*(a2xx + a3xx*x))*
     -      (b1xx + y*(b2xx + b3xx*y)))/
     -    dexp(at*t))

        Txye(i,j) = -(((-1 + betann)*
     -      (a1xy + x*(a2xy + a3xy*x))*
     -      (b1xy + y*(b2xy + b3xy*y)))/
     -    dexp(at*t))

        Tyye(i,j) = -(((-1 + betann)*
     -      (a1yy + x*(a2yy + a3yy*x))*
     -      (b1yy + y*(b2yy + b3yy*y)))/
     -    dexp(at*t))
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_02(t)

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
        tfwz(i,j) = (betann*dexp(at*t)*
     -     (-64 + a2yy*b2yy*Rey - 
     -       2*a1xy*b3xy*Rey + 384*x - 
     -       2*a3xx*b2xx*Rey*x + 
     -       2*a3yy*b2yy*Rey*x - 
     -       2*a2xy*b3xy*Rey*x - 
     -       192*x**2 - 384*x**3 + 
     -       192*x**4 + 
     -       2*Rey*
     -        (a2yy*b3yy - 
     -          2*a3xx*b3xx*x + 
     -          2*a3yy*b3yy*x)*y + 
     -       192*(1 + 12*(-1 + x)*x)*
     -        y**2 + 192*y**4 - 
     -       a2xx*Rey*(b2xx + 2*b3xx*y) + 
     -       2*a3xy*Rey*
     -        (b1xy - b3xy*x**2 + 
     -          b2xy*y + b3xy*y**2)) + 
     -    Rey*(-512*(-1 + x)*x*(-1 + 2*x)*
     -        y*((-1 + x)**2*x**2 - 
     -          (1 + (-1 + x)*x)*
     -           (1 + 2*(-1 + x)*x)*y**2 
     -           + 3*(1 + (-1 + x)*x)*
     -           (1 + 2*(-1 + x)*x)*y**4 
     -           - 2*(1 + 3*(-1 + x)*x)*
     -           y**6) + 
     -       dexp(at*t)*
     -        (-(a2yy*b2yy) + 
     -          2*a1xy*b3xy + 
     -          a2xx*(b2xx + 2*b3xx*y) - 
     -          2*a3xy*
     -           (b1xy - b3xy*x**2 + 
     -             b2xy*y + b3xy*y**2) + 
     -          2*
     -           (x*
     -              (-(a3yy*b2yy) + 
     -               a2xy*b3xy - 
     -               8*at*(-1 + x)**2*x) 
     -              - b3yy*
     -              (a2yy + 2*a3yy*x)*y + 
     -             8*at*
     -              (-1 + 
     -               6*(x - 2*x**3 + x**4)
     -               )*y**2 + 
     -             8*at*
     -              (1 + 6*(-1 + x)*x)*
     -              y**4 + 
     -             a3xx*x*
     -              (b2xx + 2*b3xx*y)))))/
     -  (dexp(2*at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) = ((64*(-1 + betann)*dexp(at*t)*
     -       (-1 + x)*x*(-1 + 2*x)*y*
     -       (-1 + 2*y**2))/Rey + 
     -    alphaG*(1 - betann)*Rey*Wi*
     -     ((a1xx + x*(a2xx + a3xx*x))**2*
     -        (b1xx + y*(b2xx + b3xx*y))**
     -         2 + 
     -       (a1xy + x*(a2xy + a3xy*x))**
     -         2*
     -        (b1xy + y*(b2xy + b3xy*y))**
     -         2) + 
     -    Wi*(16*(-1 + betann)*(-1 + x)*x*
     -        (-1 + 2*x)*
     -        (a1xx + x*(a2xx + a3xx*x))*
     -        y**2*(b2xx + 2*b3xx*y)*
     -        (-1 + y**2) + 
     -       at*(-1 + betann)*dexp(at*t)*
     -        (a1xx + x*(a2xx + a3xx*x))*
     -        (b1xx + y*(b2xx + b3xx*y)) 
     -        - 16*(-1 + betann)*
     -        (-1 + x)**2*x**2*
     -        (a2xx + 2*a3xx*x)*y*
     -        (-1 + 2*y**2)*
     -        (b1xx + y*(b2xx + b3xx*y)) 
     -        + 64*(-1 + betann)*(-1 + x)*
     -        x*(-1 + 2*x)*
     -        (a1xx + x*(a2xx + a3xx*x))*
     -        y*(-1 + 2*y**2)*
     -        (b1xx + y*(b2xx + b3xx*y)) 
     -        + 32*(-1 + betann)*
     -        (-1 + x)**2*x**2*
     -        (a1xy + x*(a2xy + a3xy*x))*
     -        (-1 + 6*y**2)*
     -        (b1xy + y*(b2xy + b3xy*y)) 
     -        + 16*(1 - betann)*xi*
     -        (4*(-1 + x)*x*(-1 + 2*x)*
     -           (a1xx + 
     -             x*(a2xx + a3xx*x))*y*
     -           (-1 + 2*y**2)*
     -           (b1xx + 
     -             y*(b2xx + b3xx*y)) - 
     -          (1 + 6*(-1 + x)*x)*
     -           (a1xy + 
     -             x*(a2xy + a3xy*x))*
     -           y**2*(-1 + y**2)*
     -           (b1xy + 
     -             y*(b2xy + b3xy*y)) + 
     -          (-1 + x)**2*x**2*
     -           (a1xy + 
     -             x*(a2xy + a3xy*x))*
     -           (-1 + 6*y**2)*
     -           (b1xy + 
     -             y*(b2xy + b3xy*y)))) + 
     -    (1 - betann)*
     -     (a1xx + x*(a2xx + a3xx*x))*
     -     (b1xx + y*(b2xx + b3xx*y))*
     -     (dexp(at*t) + 
     -       epsylon*Rey*Wi*
     -        (a1xx*
     -           (b1xx + 
     -             y*(b2xx + b3xx*y)) + 
     -          a1yy*
     -           (b1yy + 
     -             y*(b2yy + b3yy*y)) + 
     -          x*
     -           (a2yy*b1yy + 
     -             a2yy*y*
     -              (b2yy + b3yy*y) + 
     -             a2xx*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3xx*x*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3yy*x*
     -              (b1yy + 
     -               y*(b2yy + b3yy*y)))))
     -    )/dexp(2*at*t)

       tfTxy(i,j) = ((-16*(-1 + betann)*dexp(at*t)*
     -       ((-1 + x)**2*x**2 + 
     -         (-1 - 
     -            6*(-1 + x)*x*
     -             (1 + (-1 + x)*x))*y**2 
     -          + (1 + 6*(-1 + x)*x)*y**4)
     -       )/Rey + 
     -    alphaG*(1 - betann)*Rey*Wi*
     -     (a1xy + x*(a2xy + a3xy*x))*
     -     (b1xy + y*(b2xy + b3xy*y))*
     -     ((a1xx + x*(a2xx + a3xx*x))*
     -        (b1xx + y*(b2xx + b3xx*y)) 
     -        + (a1yy + 
     -          x*(a2yy + a3yy*x))*
     -        (b1yy + y*(b2yy + b3yy*y))) 
     -     - (-1 + betann)*
     -     (a1xy + x*(a2xy + a3xy*x))*
     -     (b1xy + y*(b2xy + b3xy*y))*
     -     (dexp(at*t) + 
     -       epsylon*Rey*Wi*
     -        (a1xx*
     -           (b1xx + 
     -             y*(b2xx + b3xx*y)) + 
     -          a1yy*
     -           (b1yy + 
     -             y*(b2yy + b3yy*y)) + 
     -          x*
     -           (a2yy*b1yy + 
     -             a2yy*y*
     -              (b2yy + b3yy*y) + 
     -             a2xx*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3xx*x*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3yy*x*
     -              (b1yy + 
     -               y*(b2yy + b3yy*y)))))
     -      + (-1 + betann)*Wi*
     -     (16*(-1 + x)*x*(-1 + 2*x)*
     -        (a1xy + x*(a2xy + a3xy*x))*
     -        y**2*(b2xy + 2*b3xy*y)*
     -        (-1 + y**2) - 
     -       16*(1 + 6*(-1 + x)*x)*
     -        (a1xx + x*(a2xx + a3xx*x))*
     -        y**2*(-1 + y**2)*
     -        (b1xx + y*(b2xx + b3xx*y)) 
     -        + at*dexp(at*t)*
     -        (a1xy + x*(a2xy + a3xy*x))*
     -        (b1xy + y*(b2xy + b3xy*y)) 
     -        - 16*(-1 + x)**2*x**2*
     -        (a2xy + 2*a3xy*x)*y*
     -        (-1 + 2*y**2)*
     -        (b1xy + y*(b2xy + b3xy*y)) 
     -        + 16*(-1 + x)**2*x**2*
     -        (a1yy + x*(a2yy + a3yy*x))*
     -        (-1 + 6*y**2)*
     -        (b1yy + y*(b2yy + b3yy*y)) 
     -        + 8*xi*
     -        ((-1 + x)**2*x**2 + 
     -          (-1 - 
     -             6*(-1 + x)*x*
     -              (1 + (-1 + x)*x))*y**2
     -            + (1 + 6*(-1 + x)*x)*
     -           y**4)*
     -        (a1xx*
     -           (b1xx + 
     -             y*(b2xx + b3xx*y)) + 
     -          a1yy*
     -           (b1yy + 
     -             y*(b2yy + b3yy*y)) + 
     -          x*
     -           (a2yy*b1yy + 
     -             a2yy*y*
     -              (b2yy + b3yy*y) + 
     -             a2xx*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3xx*x*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3yy*x*
     -              (b1yy + 
     -               y*(b2yy + b3yy*y)))))
     -    )/dexp(2*at*t)

       tfTyy(i,j) = ((-64*(-1 + betann)*dexp(at*t)*
     -       (-1 + x)*x*(-1 + 2*x)*y*
     -       (-1 + 2*y**2))/Rey + 
     -    alphaG*(1 - betann)*Rey*Wi*
     -     ((a1xy + x*(a2xy + a3xy*x))**2*
     -        (b1xy + y*(b2xy + b3xy*y))**
     -         2 + 
     -       (a1yy + x*(a2yy + a3yy*x))**
     -         2*
     -        (b1yy + y*(b2yy + b3yy*y))**
     -         2) + 
     -    Wi*(16*(-1 + betann)*(-1 + x)*x*
     -        (-1 + 2*x)*
     -        (a1yy + x*(a2yy + a3yy*x))*
     -        y**2*(b2yy + 2*b3yy*y)*
     -        (-1 + y**2) - 
     -       32*(-1 + betann)*
     -        (1 + 6*(-1 + x)*x)*
     -        (a1xy + x*(a2xy + a3xy*x))*
     -        y**2*(-1 + y**2)*
     -        (b1xy + y*(b2xy + b3xy*y)) 
     -        + at*(-1 + betann)*
     -        dexp(at*t)*
     -        (a1yy + x*(a2yy + a3yy*x))*
     -        (b1yy + y*(b2yy + b3yy*y)) 
     -        - 16*(-1 + betann)*
     -        (-1 + x)**2*x**2*
     -        (a2yy + 2*a3yy*x)*y*
     -        (-1 + 2*y**2)*
     -        (b1yy + y*(b2yy + b3yy*y)) 
     -        - 64*(-1 + betann)*(-1 + x)*
     -        x*(-1 + 2*x)*
     -        (a1yy + x*(a2yy + a3yy*x))*
     -        y*(-1 + 2*y**2)*
     -        (b1yy + y*(b2yy + b3yy*y)) 
     -        + 16*(1 - betann)*xi*
     -        (-((1 + 6*(-1 + x)*x)*
     -             (a1xy + 
     -               x*(a2xy + a3xy*x))*
     -             y**2*(-1 + y**2)*
     -             (b1xy + 
     -               y*(b2xy + b3xy*y))) 
     -           + (-1 + x)**2*x**2*
     -           (a1xy + 
     -             x*(a2xy + a3xy*x))*
     -           (-1 + 6*y**2)*
     -           (b1xy + 
     -             y*(b2xy + b3xy*y)) - 
     -          4*(-1 + x)*x*(-1 + 2*x)*
     -           (a1yy + 
     -             x*(a2yy + a3yy*x))*y*
     -           (-1 + 2*y**2)*
     -           (b1yy + 
     -             y*(b2yy + b3yy*y)))) + 
     -    (1 - betann)*
     -     (a1yy + x*(a2yy + a3yy*x))*
     -     (b1yy + y*(b2yy + b3yy*y))*
     -     (dexp(at*t) + 
     -       epsylon*Rey*Wi*
     -        (a1xx*
     -           (b1xx + 
     -             y*(b2xx + b3xx*y)) + 
     -          a1yy*
     -           (b1yy + 
     -             y*(b2yy + b3yy*y)) + 
     -          x*
     -           (a2yy*b1yy + 
     -             a2yy*y*
     -              (b2yy + b3yy*y) + 
     -             a2xx*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3xx*x*
     -              (b1xx + 
     -               y*(b2xx + b3xx*y)) + 
     -             a3yy*x*
     -              (b1yy + 
     -               y*(b2yy + b3yy*y)))))
     -    )/dexp(2*at*t) 

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_02(t)

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
         ux(i,j) = (16*(-1 + x)**2*x**2*y*
     -    (-1 + 2*y**2))/dexp(at*t)

         uy(i,j) = (-16*(-1 + x)*x*(-1 + 2*x)*y**2*
     -    (-1 + y**2))/dexp(at*t)

         wz(i,j) = (16*(-((-1 + x)**2*x**2) + 
     -      (-1 + 6*(x - 2*x**3 + x**4))*
     -       y**2 + 
     -      (1 + 6*(-1 + x)*x)*y**4))/
     -  dexp(at*t)

        psi(i,j) = (8*(-1 + x)**2*x**2*y**2*
     -    (-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txx(i,j) =  -(((-1 + betann)*
     -      (a1xx + x*(a2xx + a3xx*x))*
     -      (b1xx + y*(b2xx + b3xx*y)))/
     -    dexp(at*t))

        Txy(i,j) = -(((-1 + betann)*
     -      (a1xy + x*(a2xy + a3xy*x))*
     -      (b1xy + y*(b2xy + b3xy*y)))/
     -    dexp(at*t))

        Tyy(i,j) = -(((-1 + betann)*
     -      (a1yy + x*(a2yy + a3yy*x))*
     -      (b1yy + y*(b2yy + b3yy*y)))/
     -    dexp(at*t)) 
         end do
      end do
      !!! Boundary: top and botton
      do j = 1, jmax, jmax - 1
         y = y0 + dble(j-1)*dy
         do i = 1, imax
            x     = x0 + dble(i-1)*dx
         !!! Velocitys 
         ux(i,j) = (16*(-1 + x)**2*x**2*y*
     -    (-1 + 2*y**2))/dexp(at*t)

         uy(i,j) = (-16*(-1 + x)*x*(-1 + 2*x)*y**2*
     -    (-1 + y**2))/dexp(at*t)

         wz(i,j) = (16*(-((-1 + x)**2*x**2) + 
     -      (-1 + 6*(x - 2*x**3 + x**4))*
     -       y**2 + 
     -      (1 + 6*(-1 + x)*x)*y**4))/
     -  dexp(at*t)

        psi(i,j) = (8*(-1 + x)**2*x**2*y**2*
     -    (-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txx(i,j) =  -(((-1 + betann)*
     -      (a1xx + x*(a2xx + a3xx*x))*
     -      (b1xx + y*(b2xx + b3xx*y)))/
     -    dexp(at*t))

        Txy(i,j) = -(((-1 + betann)*
     -      (a1xy + x*(a2xy + a3xy*x))*
     -      (b1xy + y*(b2xy + b3xy*y)))/
     -    dexp(at*t))

        Tyy(i,j) = -(((-1 + betann)*
     -      (a1yy + x*(a2yy + a3yy*x))*
     -      (b1yy + y*(b2yy + b3yy*y)))/
     -    dexp(at*t)) 
         end do
      end do
      return
      end
