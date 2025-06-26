ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                begin of solvers                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Solve tridiagonal matrix for the derivatives in x direction using Thomas Algorithm (TDMA)
!!\f{eqnarray*}{\left[ \begin{array}{cccccccc} b_{1} & c_{1} & \quad &\quad\\a_{2} & b_{2} & c_{2} &\quad &\quad\\ \quad & a_{3} & b_{3} & \ddots & \quad \\\quad & \quad & \ddots & \ddots & c_{n-1}\\ \quad & \quad & \quad & a_{n} & b_{n} \end{array} \right] \left[ \begin{array}{c} x_{1}\\ x_{2} \\ x_{3} \\ \vdots\\ x_{n} \end{array}\right] = \left[ \begin{array}{c} d_{1}\\ d_{2} \\ d_{3} \\ \vdots\\ d_{n} \end{array} \right]\ \ \ \Rightarrow\ \ \ \begin{array}{l}c'_{i} = \left\{\begin{array}{ll}\frac{c_{i}}{b_{i}}, & i = 1\\ \frac{d_{i}-a_{i}d'_{i-1}}{b_{i}-a_{i}c'_{i-1}}, & i = 2,3,\cdots,n-1\end{array} \right. \\ \mbox{and} \\ d'_{i} = \left\{\begin{array}{ll}\frac{d_{i}}{b_{i}}, & i = 1\\ \frac{c_{i}}{b_{i}-a_{i}c'_{i-1}}, & i = 2,3,\cdots,n\end{array} \right.\end{array}\ \ \ \Rightarrow\ \ \ \begin{array}{l} x_{n} = d'_{n} \\ x_{i} = d'_{i} - c'_{i}x_{i+1},\ \ i = n-1, n-2, \cdots, 1\end{array} \f}
      subroutine tridx(a,b,c,rhs)

      implicit none
      include 'par.nn'
      integer i, j
      real*8 a(imax), b(imax), c(imax), gam(imax), bet
      real*8 rhs(imax,jmax), u(imax)

      do j = 1, jmax
        bet  = b(1)
        u(1) = rhs(1,j) / bet
        do i = 2, imax
          gam(i) = c(i-1) / bet
          bet    = b(i) - a(i) * gam(i)
          u(i)   = ( rhs(i,j) - a(i) * u(i-1) ) / bet
        end do
        do i = imax - 1, 1, -1
          u(i) = u(i) - gam(i+1) * u(i+1)
        end do
        do i = 1, imax
          rhs(i,j) = u(i)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Solves the tridiagonal matrix for the derivatives in y direction using Thomas Algorithm (TDMA)
!!\f{eqnarray*}{\left[ \begin{array}{cccccccc} b_{1} & c_{1} & \quad &\quad\\a_{2} & b_{2} & c_{2} &\quad &\quad\\ \quad & a_{3} & b_{3} & \ddots & \quad \\\quad & \quad & \ddots & \ddots & c_{n-1}\\ \quad & \quad & \quad & a_{n} & b_{n} \end{array} \right] \left[ \begin{array}{c} x_{1}\\ x_{2} \\ x_{3} \\ \vdots\\ x_{n} \end{array}\right] = \left[ \begin{array}{c} d_{1}\\ d_{2} \\ d_{3} \\ \vdots\\ d_{n} \end{array} \right]\ \ \ \Rightarrow\ \ \ \begin{array}{l}c'_{i} = \left\{\begin{array}{ll}\frac{c_{i}}{b_{i}}, & i = 1\\ \frac{d_{i}-a_{i}d'_{i-1}}{b_{i}-a_{i}c'_{i-1}}, & i = 2,3,\cdots,n-1\end{array} \right. \\ \mbox{and} \\ d'_{i} = \left\{\begin{array}{ll}\frac{d_{i}}{b_{i}}, & i = 1\\ \frac{c_{i}}{b_{i}-a_{i}c'_{i-1}}, & i = 2,3,\cdots,n\end{array} \right.\end{array}\ \ \ \Rightarrow\ \ \ \begin{array}{l} x_{n} = d'_{n} \\ x_{i} = d'_{i} - c'_{i}x_{i+1},\ \ i = n-1, n-2, \cdots, 1\end{array} \f}
      subroutine tridy(a, b, c, rhs)

      implicit none
      include 'par.nn'
      integer i, j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      real*8 rhs(imax,jmax), u(jmax)

      do i = 1, imax
        bet  = b(1)
        u(1) = rhs(i,1) / bet
        do j = 2, jmax
          gam(j) = c(j-1) / bet
          bet    = b(j) - a(j) * gam(j)
          u(j)   = ( rhs(i,j) - a(j) * u(j-1) ) / bet
        end do
        do j = jmax - 1, 1, -1
          u(j) = u(j) - gam(j+1) * u(j+1)
        end do
        do j = 1, jmax
          rhs(i,j) = u(j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Mounts the left hand side of the pentadiagonal matrix
!!\f{eqnarray*}{\left( \begin{array}{cccccc} 1 & 0 & 0 & 0 & 0 & 0 \\ 254 & - 432 - 360 \frac{\triangle y^{2}}{\triangle x^{2}} & 162 - 90 \frac{\triangle y^{2}}{\triangle x^{2}} & 16 & 0& 0 \\ 3 & 48 - 20 \frac{\triangle y^{2}}{\triangle x^{2}} & -102 - 110 \frac{\triangle y^{2}}{\triangle x^{2}} & 48 - 20 \frac{\triangle y^{2}}{\triangle x^{2}} & 3 & 0 \\ \quad & \ddots & \ddots & \ddots & \ddots & \quad \\ 0 & 3 & 48 - 20 \frac{\triangle y^{2}}{\triangle x^{2}} & -102 - 110 \frac{\triangle y^{2}}{\triangle x^{2}} & 48 - 20 \frac{\triangle y^{2}}{\triangle x^{2}} & 3 \\ 0 & 0 & 16 & 162 - 90 \frac{\triangle y^{2}}{\triangle x^{2}} & - 432 - 360 \frac{\triangle y^{2}}{\triangle x^{2}} & 254 \\ 0 & 0 & 0 & 0 & 0 & 1\end{array} \right)\ \ \ \mbox{or}\ \ \ \left( \begin{array}{cccccc} 1 & 0 & 0 & 0 & 0 & 0 \\ 254 & - 432 - 180 \frac{\triangle y^{2}}{\triangle x^{2}} & 162 - 45 \frac{\triangle y^{2}}{\triangle x^{2}} & 16 & 0& 0 \\ 3 & 48 - 10 \frac{\triangle y^{2}}{\triangle x^{2}} & -102 - 55 \frac{\triangle y^{2}}{\triangle x^{2}} & 48 - 10 \frac{\triangle y^{2}}{\triangle x^{2}} & 3 & 0 \\ \quad & \ddots & \ddots & \ddots & \ddots & \quad \\ 0 & 3 & 48 - 10 \frac{\triangle y^{2}}{\triangle x^{2}} & -102 - 55 \frac{\triangle y^{2}}{\triangle x^{2}} & 48 - 10 \frac{\triangle y^{2}}{\triangle x^{2}} & 3 \\ 0 & 0 & 16 & 162 - 45 \frac{\triangle y^{2}}{\triangle x^{2}} & - 432 - 180 \frac{\triangle y^{2}}{\triangle x^{2}} & 254 \\ 0 & 0 & 0 & 0 & 0 & 1\end{array} \right)\f}
!! that corresponds to the central or border points, respectively, resulting of the discretization on y direction
!! - \f$j = 2\f$ (and similar to \f$ j=jmax-1\f$) 
!! \f{eqnarray*}{\frac{4}{5}\frac{\partial^{2}v}{\partial y^{2}}\Big|_{i,2} + \frac{1}{5}\frac{\partial^{2}v}{\partial y^{2}}\Big|_{i,3} = \frac{254v_{i,1}-432v_{i,2}+162v_{i,3}+16v_{i,4}}{180\triangle y^{2}}\f}
!! - \f$3 \leq j \leq jmax-2\f$
!! \f{eqnarray*}{\frac{2}{15}\frac{\partial^{2}v}{\partial y^{2}}\Big|_{i,j-1} + \frac{11}{15}\frac{\partial^{2}v}{\partial y^{2}}\Big|_{i,j} + \frac{2}{15}\frac{\partial^{2}v}{\partial y^{2}}\Big|_{i,j+1} = \frac{3v_{i,j-2}+48v_{i,j-1}-102v_{i,j}+48v_{i,j+1}+3v_{i,j+2}}{60\triangle y^{2}}\f}i
      subroutine band5_poi(a, n, al, indy)

      ! solve the LHS of the pentadiagonal matrix
      implicit none
      include 'par.nn'
      integer m1, n, indy(jmax), i, j, k, l, mm
      real*8 d, a(jmax,5), al(jmax,5), TINY, dum
      parameter (TINY = 1.e-20)

      m1 = 2
      mm = 5
      l  = m1
      do i = 1, m1
        do j = m1 + 2 - i, mm
          a(i,j-l) = a(i,j)
        end do
        l = l - 1
        do j = mm - l, mm
          a(i,j) = 0.d0
        end do
      end do
      d = 1.d0
      l = m1
      do k = 1, n
        dum = a(k,1)
        i   = k
        if (l.lt.n) l = l + 1
        do j = k + 1, l
          if (dabs(a(j,1)) .gt. dabs(dum)) then
            dum = a(j,1)
            i   = j
          end if
        end do
        indy(k) = i
        if (dum .eq. 0.d0) a(k,1) = TINY
        if (i .ne. k) then
          d = - d
          do j = 1, mm
            dum    = a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = dum
          end do
        endif
        do i = k + 1, l
          dum       = a(i,1) / a(k,1)
          al(k,i-k) = dum
          do j = 2, mm
            a(i,j-1) = a(i,j) - dum * a(k,j)
          end do
          a(i,mm) = 0.d0
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> solve the the pentadiagonal matrix the term a and al comes from
!! subroutine band5_poi, the rhs variable is the input and at the
!! end is the result of matrix
      subroutine banbk5_poi(a, n, al, indy, rhs)

      implicit none
      include 'par.nn'
      integer n, i, k, l, indy(jmax)
      real*8 a(jmax,5), al(jmax,5), rhs(jmax), dum

      l = 2
      do k = 1, n
        i = indy(k)
        if (i.ne.k) then
          dum    = rhs(k)
          rhs(k) = rhs(i)
          rhs(i) = dum
        endif
        if (l.lt.n) l = l + 1
        do i = k + 1, l
          rhs(i) = rhs(i) - al(k,i-k) * rhs(k)
        end do
      end do
      l = 1
      do i = n, 1, -1
        dum = rhs(i)
        do k = 2, l
          dum = dum - a(i,k) * rhs(k+i-1)
        end do
        rhs(i) = dum / a(i,1)
        if (l.lt.5) l = l + 1
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  end of solvers                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
