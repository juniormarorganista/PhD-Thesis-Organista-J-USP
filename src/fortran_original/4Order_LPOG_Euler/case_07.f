cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_07(t)
      
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
         uxe(i,j) = -8*(-1 + x)**2*x**2*y**2*
     -  (-1 + Tanh(4 - 8*at*t))

         uye(i,j) = (16*(-1 + x)*x*(-1 + 2*x)*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         wze(i,j) = (-16*y*(y**2 + 
     -      3*(-1 + x)*x*
     -       ((-1 + x)*x + 2*y**2))*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

        psie(i,j) = (-8*(-1 + x)**2*x**2*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         !!! Tensors
        Txxe(i,j) = -(((-1 + betann)*c1)/dexp(at*t))

        Txye(i,j) = -(((-1 + betann)*c2)/dexp(at*t))

        Tyye(i,j) = -(((-1 + betann)*c3)/dexp(at*t))
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_07(t)

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
        tfwz(i,j) = (-64*y*(2*at*Rey*
     -       (y**2 + 
     -         3*(-1 + x)*x*
     -          ((-1 + x)*x + 2*y**2))*
     -       (1.d0/Cosh(4 - 8*at*t))**2 + 
     -      (-1 + Tanh(4 - 8*at*t))*
     -       (3*betann*
     -          (1 + 6*(-1 + x)*x + y**2) 
     -          + 4*Rey*(-1 + x)*x*
     -          (-1 + 2*x)*y**2*
     -          (y**2 + 
     -            (-1 + x)*x*
     -             (-2*(-1 + x)*x + 
     -              3*y**2)) + 
     -         4*Rey*(-1 + x)*x*
     -          (-1 + 2*x)*y**2*
     -          (2*(-1 + x)**2*x**2 + 
     -            (-1 - 3*(-1 + x)*x)*y**2
     -            )*Tanh(4 - 8*at*t))))/
     -  (3.*Rey)

       !!! Tensors
       tfTxx(i,j) = ((-1 + betann)*
     -    (-3*alphaG*c2**2*
     -       (dexp(8.0d0) + dexp(16*at*t))*Rey**2*
     -       Wi - 
     -      3*c1**2*(dexp(8.0d0) + dexp(16*at*t))*
     -       (alphaG + epsylon)*Rey**2*Wi 
     -       - 3*c1*Rey*
     -       (c3*dexp(8.0d0)*epsylon*Rey*Wi + 
     -         c3*dexp(16*at*t)*epsylon*
     -          Rey*Wi + 
     -         dexp(8 + at*t)*
     -          (1 - at*Wi) + 
     -         dexp(17*at*t)*
     -          (1 - at*Wi + 
     -            64*Wi*(-1 + x)*x*
     -             (-1 + 2*x)*(-1 + xi)*
     -             y**2)) + 
     -      32*dexp(17*at*t)*y*
     -       (6*dexp(at*t)*(-1 + x)*x*
     -          (-1 + 2*x)*y + 
     -         c2*Rey*Wi*
     -          (-3*(-1 + x)**2*x**2*
     -             (-2 + xi) + 
     -            (1 + 6*(-1 + x)*x)*xi*
     -             y**2))))/
     -  (3.*dexp(2*at*t)*
     -    (dexp(8.0d0) + dexp(16*at*t))*Rey)

       tfTxy(i,j) = ((-1 + betann)*
     -    (3*c2*(dexp(8.0d0) + dexp(16*at*t))*Rey*
     -       (-((c1 + c3)*
     -            (alphaG + epsylon)*Rey*
     -            Wi) + 
     -         dexp(at*t)*(-1 + at*Wi)) + 
     -      16*dexp(17*at*t)*y*
     -       (2*dexp(at*t)*
     -          (3*(-1 + x)**2*x**2 + 
     -            (-1 - 6*(-1 + x)*x)*y**2
     -            ) + 
     -         Rey*Wi*
     -          (-3*(-1 + x)**2*x**2*
     -             (-2*c3 + (c1 + c3)*xi) 
     -             + (1 + 6*(-1 + x)*x)*
     -             (-2*c1 + (c1 + c3)*xi)*
     -             y**2))))/
     -  (3.*dexp(2*at*t)*
     -    (dexp(8.0d0) + dexp(16*at*t))*Rey)

       tfTyy(i,j) =  -0.3333333333333333*
     -  ((-1 + betann)*
     -     (3*alphaG*c2**2*
     -        (dexp(8.0d0) + dexp(16*at*t))*
     -        Rey**2*Wi + 
     -       3*c3**2*
     -        (dexp(8.0d0) + dexp(16*at*t))*
     -        (alphaG + epsylon)*Rey**2*Wi
     -         + 32*dexp(17*at*t)*y*
     -        (6*dexp(at*t)*(-1 + x)*x*
     -           (-1 + 2*x)*y + 
     -          c2*Rey*Wi*
     -           (3*(-1 + x)**2*x**2*xi - 
     -             (1 + 6*(-1 + x)*x)*
     -              (-2 + xi)*y**2)) + 
     -       3*c3*Rey*
     -        (c1*dexp(8.0d0)*epsylon*Rey*Wi + 
     -          c1*dexp(16*at*t)*epsylon*
     -           Rey*Wi + 
     -          dexp(8 + at*t)*
     -           (1 - at*Wi) + 
     -          dexp(17*at*t)*
     -           (1 - at*Wi - 
     -             64*Wi*(-1 + x)*x*
     -              (-1 + 2*x)*(-1 + xi)*
     -              y**2))))/
     -   (dexp(2*at*t)*
     -     (dexp(8.0d0) + dexp(16*at*t))*Rey)

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_07(t)

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
         ux(i,j) = -8*(-1 + x)**2*x**2*y**2*
     -  (-1 + Tanh(4 - 8*at*t))

         uy(i,j) = (16*(-1 + x)*x*(-1 + 2*x)*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         wz(i,j) = (-16*y*(y**2 + 
     -      3*(-1 + x)*x*
     -       ((-1 + x)*x + 2*y**2))*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

        psi(i,j) = (-8*(-1 + x)**2*x**2*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         !!! Tensors
        Txx(i,j) = -(((-1 + betann)*c1)/dexp(at*t))

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
         ux(i,j) = -8*(-1 + x)**2*x**2*y**2*
     -  (-1 + Tanh(4 - 8*at*t))

         uy(i,j) = (16*(-1 + x)*x*(-1 + 2*x)*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         wz(i,j) = (-16*y*(y**2 + 
     -      3*(-1 + x)*x*
     -       ((-1 + x)*x + 2*y**2))*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

        psi(i,j) = (-8*(-1 + x)**2*x**2*y**3*
     -    (-1 + Tanh(4 - 8*at*t)))/3.

         !!! Tensors
        Txx(i,j) = -(((-1 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1 + betann)*c3)/dexp(at*t)) 
         end do
      end do
      return
      end
