cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_03(t)
      
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
         uxe(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2*y + 4*y**3))/dexp(at*t)

         uye(i,j) = (-32*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t)

         wze(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) + (32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t)

        psie(i,j) = (16*(x**2/2. - x**3 + x**4/2.)*
     & y**2*(-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txxe(i,j) =  ((1 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t)

        Txye(i,j) = ((1 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(at*t)

        Tyye(i,j) = ((1 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(at*t)

         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_03(t)

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
        tfwz(i,j) = (8*at*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) + (32*at*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*at*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/
     & dexp(at*t) + (32*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.)*((192*(x**2 - 2*x**3 + x**4)*
     & y)/dexp(at*t) + (32*x*(-3 + 4*x)*(-y + 2*y**3))/dexp(at*t) 
     & + (32*(1 - 3*x + 2*x**2)*(-y + 2*y**3))/dexp(at*t)))
     & /dexp(at*t) + (betann*((192*(x**2 - 2*x**3 + x**4))/
     & dexp(at*t) + (32*x*(-3 + 4*x)*(-1 + 6*y**2))/dexp(at*t) 
     & + (32*(1 - 3*x + 2*x**2)*(-1 + 6*y**2))/dexp(at*t) 
     & + (8*(2 - 12*x + 12*x**2)*(-2 + 12*y**2))/dexp(at*t) 
     & + (384*(-0.5*y**2 + y**4/2.))/dexp(at*t)))/Rey - 
     & (8*(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)*
     & ((8*(2*x - 6*x**2 + 4*x**3)*(-2 + 12*y**2))/dexp(at*t) 
     & + (128*x*(-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (64*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/dexp(at*t)))/dexp(at*t) + 
     & (32*x*(1 - 3*x + 2*x**2)*(-y + 2*y**3)*
     & ((8*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2))/dexp(at*t) 
     & + (32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/
     & dexp(at*t)))/dexp(at*t) - (8*(2*x - 6*x**2 + 4*x**3)*
     & (-2*y + 4*y**3)*((8*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) + (32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/
     & dexp(at*t)))/dexp(at*t) - ((1 - betann)*
     & (-(a3xy*pii**2*dcos(pii*x)) - a2xy*pii**2*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/dexp(at*t) 
     & + ((1 - betann)*(a2xx*pii*dcos(pii*x) - a3xx*pii*Sin(pii*x))*
     & (b2xx*pii*dcos(pii*y) - b3xx*pii*Sin(pii*y)))/dexp(at*t) - 
     & ((1 - betann)*(a2yy*pii*dcos(pii*x) - a3yy*pii*Sin(pii*x))*
     & (b2yy*pii*dcos(pii*y) - b3yy*pii*Sin(pii*y)))/dexp(at*t) + 
     & ((1 - betann)*(a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (-(b3xy*pii**2*dcos(pii*y)) - b2xy*pii**2*Sin(pii*y)))/
     & dexp(at*t)

       !!! Tensors
       tfTxx(i,j) = (-16*(1 - betann)*(2*x - 6*x**2 + 4*x**3)*
     & (-2*y + 4*y**3))/(dexp(at*t)*Rey) + (alphaG*Rey*Wi*
     & (((1 - betann)**2*(a1xx+a3xx*dcos(pii*x)+a2xx*Sin(pii*x))**2*
     & (b1xx + b3xx*dcos(pii*y) +b2xx*Sin(pii*y))**2)/dexp(2*at*t) + 
     & ((1 - betann)**2*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))**2*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y))**2)/dexp(2*at*t)))/(1 - betann) 
     & + Wi*(-((at*(1 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t)) - (32*(1 - betann)*x*
     & (1 - 3*x + 2*x**2)*(-y + 2*y**3)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/
     & dexp(2*at*t) - (8*(1 - betann)*(2*x - 6*x**2 + 4*x**3)*
     & (-2*y + 4*y**3)*(a1xx + a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*
     & (b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/dexp(2*at*t) + 
     & (8*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)*
     & (a2xx*pii*dcos(pii*x) - a3xx*pii*Sin(pii*x))*
     & (b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/dexp(2*at*t) - 
     & (16*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*
     & (a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/
     & dexp(2*at*t) - (32*(1 - betann)*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b2xx*pii*dcos(pii*y) - 
     & b3xx*pii*Sin(pii*y)))/dexp(2*at*t) + xi*((16*(1 - betann)*
     & (2*x - 6*x**2 + 4*x**3)*(-2*y + 4*y**3)*
     & (a1xx + a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*
     & (b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/dexp(2*at*t) + 
     & (8*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*
     & (a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*(b1xy + 
     & b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/dexp(2*at*t) + 
     & ((1 - betann)*((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.)
     & )/dexp(at*t) - (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.)
     & )/dexp(at*t))*(a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/dexp(at*t))) + 
     & ((1 - betann)*(a1xx + a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*
     & (b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y))*
     & (1 + (epsylon*Rey*Wi*(((1 - betann)*(a1xx + 
     & a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*(b1xx + 
     & b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/dexp(at*t) + 
     & ((1 - betann)*(a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/dexp(at*t)))/
     & (1 - betann)))/dexp(at*t)

       tfTxy(i,j) = -(((1 - betann)*((8*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) - (32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) - (32*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t)))/Rey) + (alphaG*Rey*Wi*
     & (((1 - betann)**2*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(2*at*t) + 
     & ((1 - betann)**2*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(2*at*t)))/(1 - betann) 
     & + ((1 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y))*(1 + (epsylon*Rey*Wi*
     & (((1 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t) + ((1 - betann)*
     & (a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/
     & dexp(at*t)))/(1 - betann)))/dexp(at*t) + 
     & Wi*(-(((1 - betann)*((-32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) - 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.)
     & )/dexp(at*t))*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t)) - 
     & (at*(1 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(at*t) + (8*(1 - betann)*
     & (x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)*
     & (a2xy*pii*dcos(pii*x) - a3xy*pii*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/
     & dexp(2*at*t) - (8*(1 - betann)*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(2*at*t) - (32*(1 - betann)*x*
     & (1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.)*
     & (a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (b2xy*pii*dcos(pii*y) - b3xy*pii*Sin(pii*y)))/
     & dexp(2*at*t) + xi*((4*(1 - betann)*
     & (x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*
     & (a1xx + a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*
     & (b1xx + b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/
     & dexp(2*at*t) + ((1 - betann)*((-32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) - (32*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t))*(a1xx+a3xx*dcos(pii*x)+
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/(2.*dexp(at*t)) - (32*(1 - betann)*x*
     & (1 - 3*x + 2*x**2)*(-y + 2*y**3)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(2*at*t) + (8*(1 - betann)*
     & (2*x - 6*x**2 + 4*x**3)*(-2*y + 4*y**3)*(a1xy + 
     & a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*(b1xy + 
     & b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/dexp(2*at*t) + 
     & (4*(1 - betann)*(x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*
     & (a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*(b1yy + 
     & b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/dexp(2*at*t) + 
     & ((1 - betann)*((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.)
     & )/dexp(at*t) - (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.)
     & )/dexp(at*t))*(a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/
     & (2.*dexp(at*t))))

       tfTyy(i,j) =  (64*(1 - betann)*x*(1 - 3*x + 2*x**2)*
     & (-y + 2*y**3))/(dexp(at*t)*Rey) + (alphaG*Rey*Wi*
     & (((1 - betann)**2*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))**2*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y))**2)/dexp(2*at*t) + ((1 - betann)**2*
     & (a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))**2*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y))**2)/
     & dexp(2*at*t)))/(1 - betann) + ((1 - betann)*
     & (a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y))*
     & (1 + (epsylon*Rey*Wi*(((1 - betann)*
     & (a1xx + a3xx*dcos(pii*x) + a2xx*Sin(pii*x))*(b1xx + 
     & b3xx*dcos(pii*y) + b2xx*Sin(pii*y)))/dexp(at*t) + 
     & ((1 - betann)*(a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/dexp(at*t)))/
     & (1 - betann)))/dexp(at*t) + Wi*((-2*(1 - betann)*
     & ((-32*x*(-3 + 4*x)*(-0.5*y**2 + y**4/2.))/dexp(at*t) - 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/dexp(at*t))*
     & (a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/dexp(at*t) - 
     & (at*(1 - betann)*(a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/dexp(at*t) + 
     & (32*(1 - betann)*x*(1 - 3*x + 2*x**2)*(-y + 2*y**3)*
     & (a1yy + a3yy*dcos(pii*x) + a2yy*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/
     & dexp(2*at*t) + (8*(1 - betann)*(2*x - 6*x**2 + 4*x**3)*
     & (-2*y + 4*y**3)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(2*at*t) + (8*(1 - betann)*
     & (x**2 - 2*x**3 + x**4)*(-2*y + 4*y**3)*
     & (a2yy*pii*dcos(pii*x) - a3yy*pii*Sin(pii*x))*
     & (b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/
     & dexp(2*at*t) - (32*(1 - betann)*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b2yy*pii*dcos(pii*y) - 
     & b3yy*pii*Sin(pii*y)))/dexp(2*at*t) + xi*((8*(1 - betann)*
     & (x**2 - 2*x**3 + x**4)*(-2 + 12*y**2)*(a1xy + 
     & a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*(b1xy + 
     & b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/
     & dexp(2*at*t) + ((1 - betann)*
     & ((-32*x*(-3 + 4*x)* (-0.5*y**2 + y**4/2.)
     & )/dexp(at*t) - (32*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t))*
     & (a1xy + a3xy*dcos(pii*x) + a2xy*Sin(pii*x))*
     & (b1xy + b3xy*dcos(pii*y) + b2xy*Sin(pii*y)))/
     & dexp(at*t) - (64*(1 - betann)*x*(1 - 3*x + 2*x**2)*
     & (-y + 2*y**3)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + b2yy*Sin(pii*y)))/
     & dexp(2*at*t)))

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_03(t)

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
         ux(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2*y + 4*y**3))/dexp(at*t)

         uy(i,j) = (-32*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t)

         wz(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) +(32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/
     & dexp(at*t)

        psi(i,j) = (16*(x**2/2. - x**3 + x**4/2.)*
     & y**2*(-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txx(i,j) =  ((1 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t)

        Txy(i,j) = ((1 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(at*t) 

         end do
      end do
      !!! Boundary: top and botton
      do j = 1, jmax, jmax - 1
         y = y0 + dble(j-1)*dy
         do i = 1, imax
            x     = x0 + dble(i-1)*dx

         !!! Velocitys 
         ux(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2*y + 4*y**3))/dexp(at*t)

         uy(i,j) = (-32*x*(1 - 3*x + 2*x**2)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t)

         wz(i,j) = (8*(x**2 - 2*x**3 + x**4)*
     & (-2 + 12*y**2))/dexp(at*t) +(32*x*(-3 + 4*x)*
     & (-0.5*y**2 + y**4/2.))/dexp(at*t) + 
     & (32*(1 - 3*x + 2*x**2)*(-0.5*y**2 + y**4/2.))/
     & dexp(at*t)

        psi(i,j) = (16*(x**2/2. - x**3 + x**4/2.)*
     & y**2*(-1 + y**2))/dexp(at*t)

         !!! Tensors
        Txx(i,j) =  ((1 - betann)*(a1xx + a3xx*dcos(pii*x) + 
     & a2xx*Sin(pii*x))*(b1xx + b3xx*dcos(pii*y) + 
     & b2xx*Sin(pii*y)))/dexp(at*t)

        Txy(i,j) = ((1 - betann)*(a1xy + a3xy*dcos(pii*x) + 
     & a2xy*Sin(pii*x))*(b1xy + b3xy*dcos(pii*y) + 
     & b2xy*Sin(pii*y)))/dexp(at*t)

        Tyy(i,j) = ((1 - betann)*(a1yy + a3yy*dcos(pii*x) + 
     & a2yy*Sin(pii*x))*(b1yy + b3yy*dcos(pii*y) + 
     & b2yy*Sin(pii*y)))/dexp(at*t) 

         end do
      end do
      return
      end
