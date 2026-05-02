cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_09(t)
      
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
         uxe(i,j) = 8*(1 - x)**2*x**2*y**2*
     -  (1 - Tanh(4 - 8*at*t))

         uye(i,j) = (-16*(1 - x)**2*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*(1 - x)*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

         wze(i,j) =  16*(1 - x)**2*x**2*y*
     -   (1 - Tanh(4 - 8*at*t)) + 
     -  (16*(1 - x)**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. - 
     -  (64*(1 - x)*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

        psie(i,j) = (-16*(x**2/2. - x**3 + x**4/2.)*
     -    y**3*(-1 + Tanh(4 - 8*at*t)))/3.

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
      subroutine source_term_case_09(t)

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
        tfwz(i,j) = -128*at*(1 - x)**2*x**2*y*
     -   (1.0d0/Cosh(4 - 8*at*t))**2 - 
     -  (128*at*(1 - x)**2*y**3*
     -     (1.0d0/Cosh(4 - 8*at*t))**2)/3. + 
     -  (512*at*(1 - x)*x*y**3*
     -     (1.0d0/Cosh(4 - 8*at*t))**2)/3. - 
     -  (128*at*x**2*y**3*
     -     (1.0d0/Cosh(4 - 8*at*t))**2)/3. - 
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
     -   dexp(at*t) + 
     -  (betann*(32*(1 - x)**2*y*
     -        (1 - Tanh(4 - 8*at*t)) - 
     -       128*(1 - x)*x*y*
     -        (1 - Tanh(4 - 8*at*t)) + 
     -       32*x**2*y*
     -        (1 - Tanh(4 - 8*at*t)) + 
     -       (32*(1 - x)**2 - 
     -          128*(1 - x)*x + 32*x**2)*
     -        y*(1 - Tanh(4 - 8*at*t)) + 
     -       64*y**3*
     -        (1 - Tanh(4 - 8*at*t))))/Rey
     -    - (-16*(1 - x)**2*x*y**2*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     16*(1 - x)*x**2*y**2*
     -      (1 - Tanh(4 - 8*at*t)))*
     -   (16*(1 - x)**2*x**2*y*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     (16*(1 - x)**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      - (64*(1 - x)*x*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      + (16*x**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3.) 
     -   - (16*(1 - x)**2*x**2*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     16*(1 - x)**2*y**2*
     -      (1 - Tanh(4 - 8*at*t)) - 
     -     64*(1 - x)*x*y**2*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     16*x**2*y**2*
     -      (1 - Tanh(4 - 8*at*t)))*
     -   ((-16*(1 - x)**2*x*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      + (16*(1 - x)*x**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3.) 
     -   - 8*(1 - x)**2*x**2*y**2*
     -   (32*(1 - x)**2*x*y*
     -      (1 - Tanh(4 - 8*at*t)) - 
     -     32*(1 - x)*x**2*y*
     -      (1 - Tanh(4 - 8*at*t)) - 
     -     32*(1 - x)*y**3*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     32*x*y**3*
     -      (1 - Tanh(4 - 8*at*t)))*
     -   (1 - Tanh(4 - 8*at*t)) - 
     -  16*(1 - x)**2*x*y**2*
     -   (16*(1 - x)**2*x**2*y*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     (16*(1 - x)**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      - (64*(1 - x)*x*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      + (16*x**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3.)*
     -   (1 - Tanh(4 - 8*at*t)) + 
     -  16*(1 - x)*x**2*y**2*
     -   (16*(1 - x)**2*x**2*y*
     -      (1 - Tanh(4 - 8*at*t)) + 
     -     (16*(1 - x)**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      - (64*(1 - x)*x*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3. 
     -      + (16*x**2*y**3*
     -        (1 - Tanh(4 - 8*at*t)))/3.)*
     -   (1 - Tanh(4 - 8*at*t))

       !!! Tensors
       tfTxx(i,j) = (alphaG*Rey*Wi*
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
     -   + ((1 - betann)*
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
     -        (1 - betann)))/dexp(at*t) - 
     -  (2*(1 - betann)*
     -     (16*(1 - x)**2*x*y**2*
     -        (1 - Tanh(4 - 8*at*t)) - 
     -       16*(1 - x)*x**2*y**2*
     -        (1 - Tanh(4 - 8*at*t))))/Rey
     -    + Wi*(-((at*(1 - betann)*
     -          (a1xx + a3xx*Cos(pii*x) + 
     -            a2xx*Sin(pii*x))*
     -          (b1xx + b3xx*Cos(pii*y) + 
     -            b2xx*Sin(pii*y)))/
     -        dexp(at*t)) - 
     -     (2*(1 - betann)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        (16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) - 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        (-16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) + 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b2xx*pii*Cos(pii*y) - 
     -          b3xx*pii*Sin(pii*y))*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*((2*(1 - betann)*
     -           (a1xx + 
     -             a3xx*Cos(pii*x) + 
     -             a2xx*Sin(pii*x))*
     -           (b1xx + 
     -             b3xx*Cos(pii*y) + 
     -             b2xx*Sin(pii*y))*
     -           (16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               - 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           ((-16*(1 - x)**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              + (64*(1 - x)*x*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              - (16*x**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3.)
     -           )/dexp(at*t) + 
     -        (16*(1 - betann)*(1 - x)**2*
     -           x**2*y*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        y**2*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        y**2*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*y**2*
     -        (a2xx*pii*Cos(pii*x) - 
     -          a3xx*pii*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (32*(1 - betann)*(1 - x)**2*
     -        x**2*y*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t))

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
     -   + ((1 - betann)*
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
     -        (1 - betann)))/dexp(at*t) - 
     -  ((1 - betann)*
     -     (16*(1 - x)**2*x**2*y*
     -        (1 - Tanh(4 - 8*at*t)) - 
     -       (16*(1 - x)**2*y**3*
     -          (1 - Tanh(4 - 8*at*t)))/3.
     -         + (64*(1 - x)*x*y**3*
     -          (1 - Tanh(4 - 8*at*t)))/3.
     -         - (16*x**2*y**3*
     -          (1 - Tanh(4 - 8*at*t)))/3.
     -       ))/Rey + 
     -  Wi*(-((at*(1 - betann)*
     -          (a1xy + a3xy*Cos(pii*x) + 
     -            a2xy*Sin(pii*x))*
     -          (b1xy + b3xy*Cos(pii*y) + 
     -            b2xy*Sin(pii*y)))/
     -        dexp(at*t)) - 
     -     ((1 - betann)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        (16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) - 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) - 
     -     ((1 - betann)*
     -        (a1xx + a3xx*Cos(pii*x) + 
     -          a2xx*Sin(pii*x))*
     -        (b1xx + b3xx*Cos(pii*y) + 
     -          b2xx*Sin(pii*y))*
     -        ((-16*(1 - x)**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (64*(1 - x)*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. - 
     -          (16*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b2xy*pii*Cos(pii*y) - 
     -          b3xy*pii*Sin(pii*y))*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*(((1 - betann)*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           (16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               - 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           (-16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               + 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xx + 
     -             a3xx*Cos(pii*x) + 
     -             a2xx*Sin(pii*x))*
     -           (b1xx + 
     -             b3xx*Cos(pii*y) + 
     -             b2xx*Sin(pii*y))*
     -           ((-16*(1 - x)**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              + (64*(1 - x)*x*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              - (16*x**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3.)
     -           )/(2.*dexp(at*t)) + 
     -        ((1 - betann)*
     -           (a1yy + 
     -             a3yy*Cos(pii*x) + 
     -             a2yy*Sin(pii*x))*
     -           (b1yy + 
     -             b3yy*Cos(pii*y) + 
     -             b2yy*Sin(pii*y))*
     -           ((-16*(1 - x)**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              + (64*(1 - x)*x*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              - (16*x**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3.)
     -           )/(2.*dexp(at*t)) + 
     -        (8*(1 - betann)*(1 - x)**2*
     -           x**2*y*
     -           (a1xx + 
     -             a3xx*Cos(pii*x) + 
     -             a2xx*Sin(pii*x))*
     -           (b1xx + 
     -             b3xx*Cos(pii*y) + 
     -             b2xx*Sin(pii*y))*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t) + 
     -        (8*(1 - betann)*(1 - x)**2*
     -           x**2*y*
     -           (a1yy + 
     -             a3yy*Cos(pii*x) + 
     -             a2yy*Sin(pii*x))*
     -           (b1yy + 
     -             b3yy*Cos(pii*y) + 
     -             b2yy*Sin(pii*y))*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        y**2*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        y**2*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*y**2*
     -        (a2xy*pii*Cos(pii*x) - 
     -          a3xy*pii*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)**2*
     -        x**2*y*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t))

       tfTyy(i,j) = (alphaG*Rey*Wi*
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
     -   + ((1 - betann)*
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
     -        (1 - betann)))/dexp(at*t) - 
     -  (2*(1 - betann)*
     -     (-16*(1 - x)**2*x*y**2*
     -        (1 - Tanh(4 - 8*at*t)) + 
     -       16*(1 - x)*x**2*y**2*
     -        (1 - Tanh(4 - 8*at*t))))/Rey
     -    + Wi*(-((at*(1 - betann)*
     -          (a1yy + a3yy*Cos(pii*x) + 
     -            a2yy*Sin(pii*x))*
     -          (b1yy + b3yy*Cos(pii*y) + 
     -            b2yy*Sin(pii*y)))/
     -        dexp(at*t)) - 
     -     ((1 - betann)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y))*
     -        (-16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) + 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) - 
     -     (2*(1 - betann)*
     -        (a1xy + a3xy*Cos(pii*x) + 
     -          a2xy*Sin(pii*x))*
     -        (b1xy + b3xy*Cos(pii*y) + 
     -          b2xy*Sin(pii*y))*
     -        ((-16*(1 - x)**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (64*(1 - x)*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. - 
     -          (16*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b2yy*pii*Cos(pii*y) - 
     -          b3yy*pii*Sin(pii*y))*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*((2*(1 - betann)*
     -           (a1yy + 
     -             a3yy*Cos(pii*x) + 
     -             a2yy*Sin(pii*x))*
     -           (b1yy + 
     -             b3yy*Cos(pii*y) + 
     -             b2yy*Sin(pii*y))*
     -           (-16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               + 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           ((-16*(1 - x)**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              + (64*(1 - x)*x*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3. 
     -              - (16*x**2*y**3*
     -               (1 - 
     -              Tanh(4 - 8*at*t)))/3.)
     -           )/dexp(at*t) + 
     -        (16*(1 - betann)*(1 - x)**2*
     -           x**2*y*
     -           (a1xy + 
     -             a3xy*Cos(pii*x) + 
     -             a2xy*Sin(pii*x))*
     -           (b1xy + 
     -             b3xy*Cos(pii*y) + 
     -             b2xy*Sin(pii*y))*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        y**2*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        y**2*
     -        (a1yy + a3yy*Cos(pii*x) + 
     -          a2yy*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*y**2*
     -        (a2yy*pii*Cos(pii*x) - 
     -          a3yy*pii*Sin(pii*x))*
     -        (b1yy + b3yy*Cos(pii*y) + 
     -          b2yy*Sin(pii*y))*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t)) 

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_09(t)

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
         ux(i,j) = 8*(1 - x)**2*x**2*y**2*
     -  (1 - Tanh(4 - 8*at*t))

         uy(i,j) = (-16*(1 - x)**2*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*(1 - x)*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

         wz(i,j) =  16*(1 - x)**2*x**2*y*
     -   (1 - Tanh(4 - 8*at*t)) + 
     -  (16*(1 - x)**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. - 
     -  (64*(1 - x)*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

        psi(i,j) = (-16*(x**2/2. - x**3 + x**4/2.)*
     -    y**3*(-1 + Tanh(4 - 8*at*t)))/3.

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
         ux(i,j) = 8*(1 - x)**2*x**2*y**2*
     -  (1 - Tanh(4 - 8*at*t))

         uy(i,j) = (-16*(1 - x)**2*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*(1 - x)*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

         wz(i,j) =  16*(1 - x)**2*x**2*y*
     -   (1 - Tanh(4 - 8*at*t)) + 
     -  (16*(1 - x)**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. - 
     -  (64*(1 - x)*x*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3. + 
     -  (16*x**2*y**3*
     -     (1 - Tanh(4 - 8*at*t)))/3.

        psi(i,j) = (-16*(x**2/2. - x**3 + x**4/2.)*
     -    y**3*(-1 + Tanh(4 - 8*at*t)))/3.

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
