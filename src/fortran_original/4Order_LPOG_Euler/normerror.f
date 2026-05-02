ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 subroutines writes the results        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine normerror(erru,errv,errwz,errpsi,errTxx,errTxy,errTyy)
 
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      real*8 erru, errv, errwz, errpsi, errTxx, errTxy, errTyy
      real*8 UXb(imax,jmax),  UYb(imax,jmax),  WZb(imax,jmax), 
     &      PSIb(imax,jmax), TXXb(imax,jmax), TXYb(imax,jmax),
     &      TYYb(imax,jmax)

      UXb  = abs(ux  - uxe)
      UYb  = abs(uy  - uye) 
      WZb  = abs(wz  - wze) 
      PSIb = abs(psi - psie)
      TXXb = abs(Txx - Txxe)
      TXYb = abs(Txy - Txye)
      TYYb = abs(Tyy - Tyye)
     
      call integrals_roles(erru  ,  UXb)
      call integrals_roles(errv  ,  UYb)
      call integrals_roles(errwz ,  WZb)
      call integrals_roles(errpsi, PSIb)
      call integrals_roles(errTxx, TXXb)
      call integrals_roles(errTxy, TXYb)
      call integrals_roles(errTyy, TYYb)
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine normerror_q(eu,ev,ewz,epsi,eTxx,eTxy,eTyy,q)
 
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j, q
      real*8 eu, ev, ewz, epsi, eTxx, eTxy, eTyy
      real*8 UXb(imax,jmax),  UYb(imax,jmax),  WZb(imax,jmax), 
     &      PSIb(imax,jmax), TXXb(imax,jmax), TXYb(imax,jmax),
     &      TYYb(imax,jmax)

      UXb  = abs(ux  - uxe )
      UYb  = abs(uy  - uye ) 
      WZb  = abs(wz  - wze ) 
      PSIb = abs(psi - psie)
      TXXb = abs(Txx - Txxe)
      TXYb = abs(Txy - Txye)
      TYYb = abs(Tyy - Tyye)

      eu   = 0.0d0
      ev   = 0.0d0
      ewz  = 0.0d0
      epsi = 0.0d0
      eTxx = 0.0d0
      eTxy = 0.0d0
      eTyy = 0.0d0

      do j=1,jmax
         do i=1,imax
            eu   = eu   +  UXb(i,j)**q
            ev   = ev   +  UYb(i,j)**q
            ewz  = ewz  +  WZb(i,j)**q
            epsi = epsi + PSIb(i,j)**q
            eTxx = eTxx + TXXb(i,j)**q
            eTxy = eTxy + TXYb(i,j)**q
            eTyy = eTyy + TYYb(i,j)**q
         end do
      end do
      eu   = (eu   * dx * dy)**(1.0d0/q)
      ev   = (ev   * dx * dy)**(1.0d0/q)
      ewz  = (ewz  * dx * dy)**(1.0d0/q)
      epsi = (epsi * dx * dy)**(1.0d0/q)
      eTxx = (eTxx * dx * dy)**(1.0d0/q)
      eTxy = (eTxy * dx * dy)**(1.0d0/q) 
      eTyy = (eTyy * dx * dy)**(1.0d0/q)
 
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrals_roles(JJ,F)
 
      implicit none
      include 'par.nn'
      real*8 F(imax,jmax), JJ 
      
      if ((mod(imax+2,3).eq.0).and.(mod(jmax+2,3).eq.0)) then
         call simpson_role38(JJ,F)
      else if ((mod(imax,2).eq.1).and.(mod(jmax,2).eq.1)) then
         call simpson_role13(JJ,F)
      else 
         call trapezoidal(JJ,F)
      end if
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trapezoidal(JJ,F)
 
      implicit none
      include 'par.nn'
      integer i, j
      real*8 F(imax,jmax), JJ, L, J1, J2, K1, K2

      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      do i = 1, imax
        K1 = F(i,1) + F(i,jmax)
        K2 = 0.0d0 
        do j = 2, jmax-1
           K2 = K2 + F(i,j)
        end do
        L = (K1 + 2.0d0*K2)*dy/2.0d0
        if ((i.eq.1).or.(i.eq.imax)) then
          J1 = J1 + L
        else 
          J2 = J2 + L
        end if 
      end do
      JJ =  (J1 + 2.0d0*J2)*dx/2.0d0
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpson_role13(JJ,F)
 
      implicit none
      include 'par.nn'
      integer i, j
      real*8 F(imax,jmax), JJ, L, J1, J2, J3, K1, K2, K3, Q

      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      J3 = JJ 

      do i = 1, imax
        K1 = F(i,1) + F(i,jmax)
        K2 = 0.0d0 
        K3 = 0.0d0
        do j = 2, jmax-1
           Q  = F(i,j)
           if (mod(j,2).eq.1) then
             K2 = K2 + Q
           else if (mod(j,2).eq.0) then
             K3 = K3 + Q
           end if 
        end do
        L = (K1 + 2.0d0*K2 + 4.0d0*K3) * dy / 3.0d0
        if ((i.eq.1).or.(i.eq.imax)) then
          J1 = J1 + L
        else if (mod(i,2).eq.1) then
          J2 = J2 + L
        else if (mod(i,2).eq.0) then
          J3 = J3 + L
        end if 
      end do
      JJ =  (J1 + 2.0d0*J2 + 4.0d0*J3) * dx / 3.0d0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpson_role38(JJ,F)
 
      implicit none
      include 'par.nn'
      integer i, j
      real*8 F(imax,jmax), JJ, L, J1, J2, J3, J4
      real*8 K1, K2, K3, K4, Q

      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      J3 = JJ 
      J4 = JJ 

      do i = 1, imax
        K1 = F(i,1) + F(i,jmax)
        K2 = 0.0d0 
        K3 = 0.0d0
        K4 = 0.0d0
        do j = 2, jmax-1
           Q  = F(i,j)
           if (mod(j,3).eq.2) then
             K2 = K2 + Q
           else if (mod(j,3).eq.0) then
             K3 = K3 + Q
           else if (mod(j,3).eq.1) then
             K4 = K4 + Q
           end if 
        end do
        L = 3.0d0*(K1+3.0d0*K2+3.0d0*K3+2.0d0*K4)*dy/8.0d0
        if ((i.eq.1).or.(i.eq.imax)) then
          J1 = J1 + L
        else if (mod(i,3).eq.2) then
          J2 = J2 + L
        else if (mod(i,3).eq.0) then
          J3 = J3 + L
        else if (mod(i,3).eq.1) then
          J4 = J4 + L
        end if 
      end do
      JJ = 3.0d0*(J1+3.0d0*J2+3.0d0*J3+2.0d0*J4)*dx/8.0d0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpson_role245(JJ,F)
 
      implicit none
      include 'par.nn'
      integer i, j
      real*8 F(imax,jmax), JJ, L, J1, J2, J3, J4, J5
      real*8 K1, K2, K3, K4, K5, Q

      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      J3 = JJ 
      J4 = JJ 
      J5 = JJ 

      do i = 1, imax
        K1 = F(i,1) + F(i,jmax)
        K2 = 0.0d0 
        K3 = 0.0d0
        K4 = 0.0d0
        K5 = 0.0d0
        do j = 2, jmax-1
           Q  = F(i,j)
           if (mod(j,4).eq.2) then
             K2 = K2 + Q
           else if (mod(j,4).eq.3) then
             K3 = K3 + Q
           else if (mod(j,4).eq.0) then
             K4 = K4 + Q
           else if (mod(j,4).eq.1) then
             K5 = K5 + Q
           end if 
        end do
        L = 2.0d0*(K1+32.0d0*K2+12.0d0*K3+32.0d0*K4+14.d0*K5)*dy/45.0d0
        if ((i.eq.1).or.(i.eq.imax)) then
          J1 = J1 + L
        else if (mod(i,4).eq.2) then
          J2 = J2 + L
        else if (mod(i,4).eq.3) then
          J3 = J3 + L
        else if (mod(i,4).eq.0) then
          J4 = J4 + L
        else if (mod(i,4).eq.1) then
          J5 = J5 + L
        end if 
      end do
      JJ = 2.0d0*(J1+32.0d0*J2+12.0d0*J3+32.0d0*J4+14.d0*J5)*dx/45.0d0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpson_role_13_and_38_rxy2(JJ,F,rx,ry)
 
      implicit none
      include 'par.nn'
      integer i, j, rx, ry
      real*8 F(imax,jmax), JJ, L, J1, J2, J3, J4
      real*8 K1, K2, K3, K4, Q
      
      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      J3 = JJ 
      J4 = JJ 
      do i = 1, imax-rx
        K1 = F(i,1) + F(i,jmax-ry)
        K2 = 0.0d0 
        K3 = 0.0d0
        K4 = 0.0d0
        do j = 2, jmax-ry-1
           Q = F(i,j)
           if (mod(j,3).eq.2) then
             K2 = K2 + Q
           else if (mod(j,3).eq.0) then
             K3 = K3 + Q
           else if (mod(j,3).eq.1) then
             K4 = K4 + Q
           end if 
        end do
        L = 3.0d0*(K1+3.0d0*K2+3.0d0*K3+2.0d0*K4)*dy/8.0d0

        !!! ***  ***  ***  ***  *** !!!
        L = L + (F(i,jmax-2)+2.0d0*F(i,jmax-1)
     &        + 4.0d0*F(i,jmax))*dy/3.0d0
        !!! ***  ***  ***  ***  *** !!!

        if ((i.eq.1).or.(i.eq.(imax-rx))) then
          J1 = J1 + L
        else if (mod(i,3).eq.2) then
          J2 = J2 + L
        else if (mod(i,3).eq.0) then
          J3 = J3 + L
        else if (mod(i,3).eq.1) then
          J4 = J4 + L
        end if
      end do
      JJ = 3.0d0*(J1+3.0d0*J2+3.0d0*J3+2.0d0*J4)*dx/8.0d0

      !!! ***  ***  ***  ***  *** !!!
      J1   = 0.0d0
      J2   = 0.0d0
      J3   = 0.0d0
      do i = imax-rx, imax
        K1 = F(i,1) + F(i,jmax-ry)
        K2 = 0.0d0 
        K3 = 0.0d0
        K4 = 0.0d0
        do j = 2, jmax-ry-1
           Q = F(i,j)
           if (mod(j,3).eq.2) then
             K2 = K2 + Q
           else if (mod(j,3).eq.0) then
             K3 = K3 + Q
           else if (mod(j,3).eq.1) then
             K4 = K4 + Q
           end if 
        end do
        L = 3.0d0*(K1+3.0d0*K2+3.0d0*K3+2.0d0*K4)*dy/8.0d0
        L = L + (F(i,jmax-2)+2.0d0*F(i,jmax-1)
     &        + 4.0d0*F(i,jmax))*dy/3.0d0
        if ((i.eq.(imax-rx)).or.(i.eq.imax)) then
          J1 = J1 + L
        else if (mod(i,2).eq.1) then
          J2 = J2 + L
        else if (mod(i,2).eq.0) then
          J3 = J3 + L
        end if 
      end do
      JJ = JJ+(J1+2.0d0*J2+4.0d0*J3)*dx/3.0d0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpson_role_13_and_38_rxy1(JJ,F,rx,ry)
 
      implicit none
      include 'par.nn'
      integer i, j, rx, ry
      real*8 F(imax,jmax), JJ, L, J1, J2, J3, J4
      real*8 K1, K2, K3, K4, Q
      
      JJ = 0.0d0 
      J1 = JJ 
      J2 = JJ 
      J3 = JJ 
      J4 = JJ 
      do i = 1, imax-rx-3
         K1 = F(i,1) + F(i,jmax-ry-3)
         K2 = 0.0d0 
         K3 = 0.0d0
         K4 = 0.0d0
         do j = 2, jmax-ry-3-1
            Q = F(i,j)
            if (mod(j,3).eq.2) then
              K2 = K2 + Q
            else if (mod(j,3).eq.0) then
              K3 = K3 + Q
            else if (mod(j,3).eq.1) then
              K4 = K4 + Q
            end if 
         end do
         L = 3.0d0*(K1+3.0d0*K2+3.0d0*K3+2.0d0*K4)*dy/8.0d0

         !!! ***  ***  ***  ***  *** !!!
         K1 = F(i,jmax-ry-3) + F(i,jmax)
         K2 = 0.0d0 
         K3 = 0.0d0
         do j = jmax-ry-4, jmax-1
            Q  = F(i,j)
            if (mod(j,2).eq.0) then
              K2 = K2 + Q
            else if (mod(j,2).eq.1) then
              K3 = K3 + Q
            end if 
         end do
         L = L+(K1+2.0d0*K2+4.0d0*K3)*dy/3.0d0
         !!! ***  ***  ***  ***  *** !!!

         if ((i.eq.1).or.(i.eq.(imax-rx-3))) then
           J1 = J1 + L
         else if (mod(i,3).eq.2) then
           J2 = J2 + L
         else if (mod(i,3).eq.0) then
           J3 = J3 + L
         else if (mod(i,3).eq.1) then
           J4 = J4 + L
         end if
       end do
       JJ = 3.0d0*(J1+3.0d0*J2+3.0d0*J3+2.0d0*J4)*dx/8.0d0

       !!! ***  ***  ***  ***  ***  ***  ***  ***  *** !!!
       !!! ***  ***            ***            ***  *** !!!
       J1   = 0.0d0
       J2   = 0.0d0
       J3   = 0.0d0
       do i = imax-rx-3, imax
         K1 = F(i,1) + F(i,jmax-ry-3)
         K2 = 0.0d0 
         K3 = 0.0d0
         K4 = 0.0d0
         do j = 2, jmax-ry-3-1
            Q = F(i,j)
            if (mod(j,3).eq.2) then
              K2 = K2 + Q
            else if (mod(j,3).eq.0) then
              K3 = K3 + Q
            else if (mod(j,3).eq.1) then
              K4 = K4 + Q
            end if 
         end do
         L = 3.0d0*(K1+3.0d0*K2+3.0d0*K3+2.0d0*K4)*dy/8.0d0

         K1 = F(i,jmax-ry-3) + F(i,jmax)
         K2 = 0.0d0 
         K3 = 0.0d0
         do j = jmax-ry-3-1, jmax-1
            Q  = F(i,j)
            if (mod(j,2).eq.0) then
              K2 = K2 + Q
            else if (mod(j,2).eq.1) then
              K3 = K3 + Q
            end if 
         end do
         L = L+(K1+2.0d0*K2+4.0d0*K3)*dy/3.0d0

         if ((i.eq.(imax-rx-3)).or.(i.eq.imax)) then
           J1 = J1 + L
         else if (mod(i,2).eq.1) then
           J2 = J2 + L
         else if (mod(i,2).eq.0) then
           J3 = J3 + L
         end if 
       end do
       !!! ***  ***            ***            ***  *** !!!
       !!! ***  ***  ***  ***  ***  ***  ***  ***  *** !!!
       JJ = JJ+(J1+2.0d0*J2+4.0d0*J3)*dx/3.0d0
      return
      end
