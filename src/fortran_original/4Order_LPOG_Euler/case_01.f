cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_01(t)
      
      implicit none
      include 'par.nn'
      include 'comm.vare'
      integer i, j
      real*8 x, y, t

      !write(*,*) 'Sol. anal.: ', t
      !write(*,*) '********************'
      do j = 1, jmax
         y = y0 + dble(j-1)*dy
         do i = 1, imax
             x    = x0 + dble(i-1)*dx

         !!! Velocitys 
         uxe(i,j) = (16.0d0*(-1.0d0 + x)**2*x**2*y*
     & (-1.0d0 + 2.0d0*y**2))/dexp(at*t)

         uye(i,j) = (-16.0d0*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y**2*(-1.0d0 + y**2))/dexp(at*t)

         wze(i,j) = (16.0d0*(-((-1.0d0 + x)**2*x**2) + 
     & (-1.0d0 + 6.0d0*(x - 2.0d0*x**3 + x**4))*
     & y**2 + (1.0d0 + 6.0d0*(-1 + x)*x)*y**4))/
     & dexp(at*t)

        psie(i,j) = (8.0d0*(-1.0d0 + x)**2*x**2*y**2*
     & (-1.0d0 + y**2))/dexp(at*t)

         !!! Tensors
        Txxe(i,j) = -(((-1.0d0 + betann)*c1)/dexp(at*t))

        Txye(i,j) = -(((-1.0d0 + betann)*c2)/dexp(at*t))

        Tyye(i,j) = -(((-1.0d0 + betann)*c3)/dexp(at*t))

         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_01(t)

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
        tfwz(i,j) = (16.0d0*(4.0d0*betann*dexp(at*t)*
     & (-1.0d0 + 3.0d0*y**2 + 3.0d0*((-2.0d0 + x)*(-1.0d0 + x)*x*
     & (1.0d0 + x) + 12.0d0*(-1.0d0 + x)*x*y**2 + y**4)) + 
     & Rey*(at*dexp(at*t)*(-((-1.0d0 + x)**2*x**2) + 
     & (-1.0d0 + 6.0d0*(x - 2.0d0*x**3 + x**4))*y**2 + 
     & (1.0d0 + 6.0d0*(-1.0d0 + x)*x)*y**4) + 32.0d0*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y*(-((-1.0d0 + x)**2*x**2) + 
     & (1.0d0 + (-1.0d0 + x)*x)*(1.0d0 + 2.0d0*(-1.0d0 + x)*x)*y**2
     & - 3.0d0*(1.0d0 + (-1.0d0 + x)*x)*(1.0d0 + 
     & 2.0d0*(-1.0d0 + x)*x)*y**4 + 2.0d0*
     & (1.0d0 + 3.0d0*(-1.0d0 + x)*x)*y**6))))/
     & (dexp(2.0d0*at*t)*Rey)

       !!! Tensors
       tfTxx(i,j) = ((-1.0d0 + betann)*(-(alphaG*(c1**2 + c2**2)*
     & Rey*Wi) - c1*(dexp(at*t) + (c1 + c3)*epsylon*Rey*Wi)
     & + (64.0d0*dexp(at*t)*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y*(-1.0d0 + 2.0d0*y**2))
     & /Rey + Wi*(at*c1*dexp(at*t) + 64.0d0*c1*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*(-1.0d0 + xi)*y*
     & (1.0d0 - 2.0d0*y**2) + 16.0d0*c2*((-1.0d0 + x)**2*x**2*
     & (-2.0d0 + xi) - (xi + 6.0d0*(-1.0d0 + x)*x*((-1.0d0 + x)*
     & x*(-2.0d0 + xi) + xi))*y**2
     & + (1.0d0 + 6.0d0*(-1.0d0 + x)*x)*
     & xi*y**4))))/dexp(2*at*t)

       tfTxy(i,j) = ((-1.0d0 + betann)*(c2*Rey*(-((c1 + c3)*
     & (alphaG + epsylon)*Rey*Wi) + dexp(at*t)*(-1.0d0 + at*Wi)) +
     & 8.0d0*(2.0d0*dexp(at*t)*(-((-1.0d0 + x)**2*x**2) + (1.0d0 + 
     & 6.0d0*(-1.0d0 + x)*x*(1.0d0 + (-1.0d0 + x)*x))*y**2 + 
     & (-1.0d0 - 6.0d0*(-1.0d0 + x)*x)*y**4) + Rey*Wi*
     & ((-1.0d0 + x)**2*x**2*(c3*(-2.0d0 + xi) + c1*xi)
     & - (2.0d0*(-6.0d0*c3*(-1.0d0 + x)**2*x**2 + 
     & c1*(-1.0d0 - 6.0d0*(-1.0d0 + x)*x)) + (c1 + c3)*
     & (1.0d0 + 6.0d0*(-1.0d0 + x)*x*(1.0d0 + (-1.0d0 + x)*x))*xi)
     & *y**2 + (1.0d0 + 6.0d0*(-1.0d0 + x)*x)*
     & (c1*(-2.0d0 + xi) + c3*xi)*y**4))))/(dexp(2.0d0*at*t)*Rey)

       tfTyy(i,j) =  ((-1.0d0 + betann)*(-(alphaG*c2**2*Rey**2*Wi)-
     & c3**2*(alphaG + epsylon)*Rey**2*Wi - 
     & 64.0d0*dexp(at*t)*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y*(-1.0d0 + 2.0d0*y**2) + 
     & 16.0d0*c2*Rey*Wi*((-1.0d0 + x)**2*x**2*xi - 
     & (-2.0d0 + xi + 6.0d0*(-1.0d0 + x)*x*
     & (-2.0d0 + xi + (-1.0d0 + x)*x*xi))*y**2
     & + (1.0d0 + 6.0d0*(-1.0d0 + x)*x)*
     & (-2.0d0 + xi)*y**4) + c3*Rey*
     & (-(c1*epsylon*Rey*Wi) + dexp(at*t)*(-1.0d0 + at*Wi) + 
     & 64.0d0*Wi*(-1.0d0 + x)*x*(-1.0d0 + 2.0d0*x)*(-1.0d0 + xi)*y*
     & (-1.0d0 + 2.0d0*y**2))))/(dexp(2.0d0*at*t)*Rey)

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_01(t)

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
         ux(i,j) = (16.0d0*(-1.0d0 + x)**2*x**2*y*
     & (-1.0d0 + 2.0d0*y**2))/dexp(at*t)

         uy(i,j) = (-16.0d0*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y**2*(-1.0d0 + y**2))/dexp(at*t)

         wz(i,j) = (16.0d0*(-((-1.0d0 + x)**2*x**2) + 
     & (-1.0d0 + 6.0d0*(x - 2.0d0*x**3 + x**4))*
     & y**2 + (1.0d0 + 6.0d0*(-1 + x)*x)*y**4))/
     & dexp(at*t)

        psi(i,j) = (8.0d0*(-1.0d0 + x)**2*x**2*y**2*
     & (-1.0d0 + y**2))/dexp(at*t)

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
         ux(i,j) = (16.0d0*(-1.0d0 + x)**2*x**2*y*
     & (-1.0d0 + 2.0d0*y**2))/dexp(at*t)

         uy(i,j) = (-16.0d0*(-1.0d0 + x)*x*
     & (-1.0d0 + 2.0d0*x)*y**2*(-1.0d0 + y**2))/dexp(at*t)

         wz(i,j) = (16.0d0*(-((-1.0d0 + x)**2*x**2) + 
     & (-1.0d0 + 6.0d0*(x - 2.0d0*x**3 + x**4))*
     & y**2 + (1.0d0 + 6.0d0*(-1 + x)*x)*y**4))/
     & dexp(at*t)

        psi(i,j) = (8.0d0*(-1.0d0 + x)**2*x**2*y**2*
     & (-1.0d0 + y**2))/dexp(at*t)

         !!! Tensors
        Txx(i,j) = -(((-1.0d0 + betann)*c1)/dexp(at*t))

        Txy(i,j) = -(((-1.0d0 + betann)*c2)/dexp(at*t))

        Tyy(i,j) = -(((-1.0d0 + betann)*c3)/dexp(at*t))

         end do
      end do
      return
      end
