cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions_case_08(t)
      
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

         wze(i,j) = 16*(1 - x)**2*x**2*y*
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
     -    (a1xx + a2xx*x + a3xx*x**2)*
     -    (b1xx + b2xx*y + b3xx*y**2))/
     -  dexp(at*t) 

        Txye(i,j) = ((1 - betann)*
     -    (a1xy + a2xy*x + a3xy*x**2)*
     -    (b1xy + b2xy*y + b3xy*y**2))/
     -  dexp(at*t)

        Tyye(i,j) = ((1 - betann)*
     -    (a1yy + a2yy*x + a3yy*x**2)*
     -    (b1yy + b2yy*y + b3yy*y**2))/
     -  dexp(at*t)
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term_case_08(t)

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
     -     (a1xy + a2xy*x + a3xy*x**2))/
     -   dexp(at*t) + 
     -  ((1 - betann)*(a2xx + 2*a3xx*x)*
     -     (b2xx + 2*b3xx*y))/dexp(at*t) - 
     -  ((1 - betann)*(a2yy + 2*a3yy*x)*
     -     (b2yy + 2*b3yy*y))/dexp(at*t) - 
     -  (2*a3xy*(1 - betann)*
     -     (b1xy + b2xy*y + b3xy*y**2))/
     -   dexp(at*t) - 
     -  128*at*(1 - x)**2*x**2*y*
     -   (1.0d0/Cosh(4 - 8*at*t))**2 - 
     -  (128*at*(1 - x)**2*y**3*
     -     (1.0d0/Cosh(4 - 8*at*t))**2)/3. + 
     -  (512*at*(1 - x)*x*y**3*
     -     (1.0d0/Cosh(4 - 8*at*t))**2)/3. - 
     -  (128*at*x**2*y**3*
     -     (1.d0/Cosh(4 - 8*at*t))**2)/3. + 
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
     -          (a1xx + a2xx*x + 
     -             a3xx*x**2)**2*
     -          (b1xx + b2xx*y + 
     -             b3xx*y**2)**2)/
     -        dexp(2*at*t) + 
     -       ((1 - betann)**2*
     -          (a1xy + a2xy*x + 
     -             a3xy*x**2)**2*
     -          (b1xy + b2xy*y + 
     -             b3xy*y**2)**2)/
     -        dexp(2*at*t)))/(1 - betann) 
     -   + ((1 - betann)*
     -     (a1xx + a2xx*x + a3xx*x**2)*
     -     (b1xx + b2xx*y + b3xx*y**2)*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + a2xx*x + 
     -               a3xx*x**2)*
     -               (b1xx + b2xx*y + 
     -               b3xx*y**2))/dexp(at*t)
     -              + ((1 - betann)*
     -               (a1yy + a2yy*x + 
     -               a3yy*x**2)*
     -               (b1yy + b2yy*y + 
     -               b3yy*y**2))/dexp(at*t)
     -            ))/(1 - betann)))/
     -   dexp(at*t) - 
     -  (2*(1 - betann)*
     -     (16*(1 - x)**2*x*y**2*
     -        (1 - Tanh(4 - 8*at*t)) - 
     -       16*(1 - x)*x**2*y**2*
     -        (1 - Tanh(4 - 8*at*t))))/Rey
     -    + Wi*(-((at*(1 - betann)*
     -          (a1xx + a2xx*x + 
     -            a3xx*x**2)*
     -          (b1xx + b2xx*y + 
     -            b3xx*y**2))/dexp(at*t)) 
     -      - (2*(1 - betann)*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        (b1xx + b2xx*y + b3xx*y**2)*
     -        (16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) - 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        (b1xx + b2xx*y + b3xx*y**2)*
     -        (-16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) + 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) + 
     -     ((1 - betann)*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        (b2xx + 2*b3xx*y)*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*((2*(1 - betann)*
     -           (a1xx + a2xx*x + 
     -             a3xx*x**2)*
     -           (b1xx + b2xx*y + 
     -             b3xx*y**2)*
     -           (16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               - 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
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
     -           x**2*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*y*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*(a2xx + 2*a3xx*x)*y**2*
     -        (b1xx + b2xx*y + b3xx*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        y**2*
     -        (b1xx + b2xx*y + b3xx*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        y**2*
     -        (b1xx + b2xx*y + b3xx*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (32*(1 - betann)*(1 - x)**2*
     -        x**2*
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        y*(b1xy + b2xy*y + 
     -          b3xy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t))

       tfTxy(i,j) = (alphaG*Rey*Wi*
     -     (((1 - betann)**2*
     -          (a1xx + a2xx*x + 
     -            a3xx*x**2)*
     -          (a1xy + a2xy*x + 
     -            a3xy*x**2)*
     -          (b1xx + b2xx*y + 
     -            b3xx*y**2)*
     -          (b1xy + b2xy*y + 
     -            b3xy*y**2))/dexp(2*at*t) 
     -        + ((1 - betann)**2*
     -          (a1xy + a2xy*x + 
     -            a3xy*x**2)*
     -          (a1yy + a2yy*x + 
     -            a3yy*x**2)*
     -          (b1xy + b2xy*y + 
     -            b3xy*y**2)*
     -          (b1yy + b2yy*y + 
     -            b3yy*y**2))/dexp(2*at*t))
     -     )/(1 - betann) + 
     -  ((1 - betann)*
     -     (a1xy + a2xy*x + a3xy*x**2)*
     -     (b1xy + b2xy*y + b3xy*y**2)*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + a2xx*x + 
     -               a3xx*x**2)*
     -               (b1xx + b2xx*y + 
     -               b3xx*y**2))/dexp(at*t)
     -              + ((1 - betann)*
     -               (a1yy + a2yy*x + 
     -               a3yy*x**2)*
     -               (b1yy + b2yy*y + 
     -               b3yy*y**2))/dexp(at*t)
     -            ))/(1 - betann)))/
     -   dexp(at*t) - 
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
     -          (a1xy + a2xy*x + 
     -            a3xy*x**2)*
     -          (b1xy + b2xy*y + 
     -            b3xy*y**2))/dexp(at*t)) 
     -      - ((1 - betann)*
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        (b1xy + b2xy*y + b3xy*y**2)*
     -        (16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) - 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) - 
     -     ((1 - betann)*
     -        (a1xx + a2xx*x + a3xx*x**2)*
     -        (b1xx + b2xx*y + b3xx*y**2)*
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
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        (b2xy + 2*b3xy*y)*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*(((1 - betann)*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
     -           (16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               - 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
     -           (-16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               + 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xx + a2xx*x + 
     -             a3xx*x**2)*
     -           (b1xx + b2xx*y + 
     -             b3xx*y**2)*
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
     -           (a1yy + a2yy*x + 
     -             a3yy*x**2)*
     -           (b1yy + b2yy*y + 
     -             b3yy*y**2)*
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
     -           x**2*
     -           (a1xx + a2xx*x + 
     -             a3xx*x**2)*y*
     -           (b1xx + b2xx*y + 
     -             b3xx*y**2)*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t) + 
     -        (8*(1 - betann)*(1 - x)**2*
     -           x**2*
     -           (a1yy + a2yy*x + 
     -             a3yy*x**2)*y*
     -           (b1yy + b2yy*y + 
     -             b3yy*y**2)*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*(a2xy + 2*a3xy*x)*y**2*
     -        (b1xy + b2xy*y + b3xy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        y**2*
     -        (b1xy + b2xy*y + b3xy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        y**2*
     -        (b1xy + b2xy*y + b3xy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)**2*
     -        x**2*
     -        (a1yy + a2yy*x + a3yy*x**2)*
     -        y*(b1yy + b2yy*y + 
     -          b3yy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t))

       tfTyy(i,j) = (alphaG*Rey*Wi*
     -     (((1 - betann)**2*
     -          (a1xy + a2xy*x + 
     -             a3xy*x**2)**2*
     -          (b1xy + b2xy*y + 
     -             b3xy*y**2)**2)/
     -        dexp(2*at*t) + 
     -       ((1 - betann)**2*
     -          (a1yy + a2yy*x + 
     -             a3yy*x**2)**2*
     -          (b1yy + b2yy*y + 
     -             b3yy*y**2)**2)/
     -        dexp(2*at*t)))/(1 - betann) 
     -   + ((1 - betann)*
     -     (a1yy + a2yy*x + a3yy*x**2)*
     -     (b1yy + b2yy*y + b3yy*y**2)*
     -     (1 + (epsylon*Rey*Wi*
     -          (((1 - betann)*
     -               (a1xx + a2xx*x + 
     -               a3xx*x**2)*
     -               (b1xx + b2xx*y + 
     -               b3xx*y**2))/dexp(at*t)
     -              + ((1 - betann)*
     -               (a1yy + a2yy*x + 
     -               a3yy*x**2)*
     -               (b1yy + b2yy*y + 
     -               b3yy*y**2))/dexp(at*t)
     -            ))/(1 - betann)))/
     -   dexp(at*t) - 
     -  (2*(1 - betann)*
     -     (-16*(1 - x)**2*x*y**2*
     -        (1 - Tanh(4 - 8*at*t)) + 
     -       16*(1 - x)*x**2*y**2*
     -        (1 - Tanh(4 - 8*at*t))))/Rey
     -    + Wi*(-((at*(1 - betann)*
     -          (a1yy + a2yy*x + 
     -            a3yy*x**2)*
     -          (b1yy + b2yy*y + 
     -            b3yy*y**2))/dexp(at*t)) 
     -      - ((1 - betann)*
     -        (a1yy + a2yy*x + a3yy*x**2)*
     -        (b1yy + b2yy*y + b3yy*y**2)*
     -        (-16*(1 - x)**2*x*y**2*
     -           (1 - Tanh(4 - 8*at*t)) + 
     -          16*(1 - x)*x**2*y**2*
     -           (1 - Tanh(4 - 8*at*t))))/
     -      dexp(at*t) - 
     -     (2*(1 - betann)*
     -        (a1xy + a2xy*x + a3xy*x**2)*
     -        (b1xy + b2xy*y + b3xy*y**2)*
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
     -        (a1yy + a2yy*x + a3yy*x**2)*
     -        (b2yy + 2*b3yy*y)*
     -        ((-16*(1 - x)**2*x*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3. + 
     -          (16*(1 - x)*x**2*y**3*
     -             (1 - Tanh(4 - 8*at*t)))
     -            /3.))/dexp(at*t) + 
     -     xi*((2*(1 - betann)*
     -           (a1yy + a2yy*x + 
     -             a3yy*x**2)*
     -           (b1yy + b2yy*y + 
     -             b3yy*y**2)*
     -           (-16*(1 - x)**2*x*y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -               + 16*(1 - x)*x**2*
     -              y**2*
     -              (1 - Tanh(4 - 8*at*t))
     -             ))/dexp(at*t) + 
     -        ((1 - betann)*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
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
     -           x**2*
     -           (a1xy + a2xy*x + 
     -             a3xy*x**2)*y*
     -           (b1xy + b2xy*y + 
     -             b3xy*y**2)*
     -           (1 - Tanh(4 - 8*at*t)))/
     -         dexp(at*t)) + 
     -     (8*(1 - betann)*(1 - x)**2*
     -        x**2*(a2yy + 2*a3yy*x)*y**2*
     -        (b1yy + b2yy*y + b3yy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) + 
     -     (16*(1 - betann)*(1 - x)**2*x*
     -        (a1yy + a2yy*x + a3yy*x**2)*
     -        y**2*
     -        (b1yy + b2yy*y + b3yy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t) - 
     -     (16*(1 - betann)*(1 - x)*x**2*
     -        (a1yy + a2yy*x + a3yy*x**2)*
     -        y**2*
     -        (b1yy + b2yy*y + b3yy*y**2)*
     -        (1 - Tanh(4 - 8*at*t)))/
     -      dexp(at*t)) 

        end do  
      end do  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions_case_08(t)

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

         wz(i,j) = 16*(1 - x)**2*x**2*y*
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
     -    (a1xx + a2xx*x + a3xx*x**2)*
     -    (b1xx + b2xx*y + b3xx*y**2))/
     -  dexp(at*t) 

        Txy(i,j) = ((1 - betann)*
     -    (a1xy + a2xy*x + a3xy*x**2)*
     -    (b1xy + b2xy*y + b3xy*y**2))/
     -  dexp(at*t)

        Tyy(i,j) = ((1 - betann)*
     -    (a1yy + a2yy*x + a3yy*x**2)*
     -    (b1yy + b2yy*y + b3yy*y**2))/
     -  dexp(at*t) 

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

         wz(i,j) = 16*(1 - x)**2*x**2*y*
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
     -    (a1xx + a2xx*x + a3xx*x**2)*
     -    (b1xx + b2xx*y + b3xx*y**2))/
     -  dexp(at*t) 

        Txy(i,j) = ((1 - betann)*
     -    (a1xy + a2xy*x + a3xy*x**2)*
     -    (b1xy + b2xy*y + b3xy*y**2))/
     -  dexp(at*t)

        Tyy(i,j) = ((1 - betann)*
     -    (a1yy + a2yy*x + a3yy*x**2)*
     -    (b1yy + b2yy*y + b3yy*y**2))/
     -  dexp(at*t) 

         end do
      end do
      return
      end
