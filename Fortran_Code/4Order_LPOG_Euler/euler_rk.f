cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Starts the simulation process after verifying the number of points in 
!! the x and y directions. After it initializes the variables, the solver 
!! function computes the Runge-Kutta method and verifies the convergence 
!! criteria for vorticity. At the end, this function writes the results.
      subroutine solver_rk4

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j, t, div, auxdiv
      real*8 dt2, dt6, bdfc(imax), dv1z(imax,jmax), dv2z(imax,jmax),
     &       wz1(imax,jmax), Txx1(imax,jmax), Txy1(imax,jmax),
     &       Tyy1(imax,jmax), 
     &       dv1Txx(imax,jmax), dv2Txx(imax,jmax), 
     &       dv1Txy(imax,jmax), dv2Txy(imax,jmax),
     &       dv1Tyy(imax,jmax), dv2Tyy(imax,jmax), time, emv,
     &       errowz , erroTxx, erroTxy, erroTyy
      common/bd/ bdfc
      
      call init(dv1z, dv2z, dv1Txx, dv2Txx, dv1Txy, dv2Txy, 
     &        dv1Tyy, dv2Tyy, wz1, Txx1, Txy1, Tyy1, dt2, dt6, t)

      time   = 0.d0
      call boundaries_condictions(time)
      call source_term(time)
      
      t        = 1
      auxdiv   = 1

      do while (time .lt. abs(timef-dt))

        time = dble(t - 1) * dt
        wz1  = wz
        Txx1 = Txx
        Txy1 = Txy
        Tyy1 = Tyy

        !!! First Runge-Kutta step !!!
        call drv(dv1z, dv1Txx, dv1Txy, dv1Tyy) !!! K1 !!! 
        do j = 2, jmax - 1
          do i = 2, imax - 1
              wz(i,j) =  wz1(i,j) +   dv1z(i,j) * dt2
             Txx(i,j) = Txx1(i,j) + dv1Txx(i,j) * dt2
             Txy(i,j) = Txy1(i,j) + dv1Txy(i,j) * dt2
             Tyy(i,j) = Tyy1(i,j) + dv1Tyy(i,j) * dt2
          end do
        end do
        call source_term(time + dt2)
        call boundaries_condictions(time + dt2)
        call loop(1)

        !!! Second Runge-Kutta step !!!
        call drv(dv2z, dv2Txx, dv2Txy, dv2Tyy) !!! K2 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
              wz(i,j) =  wz1(i,j) +   dv2z(i,j) * dt2
             Txx(i,j) = Txx1(i,j) + dv2Txx(i,j) * dt2
             Txy(i,j) = Txy1(i,j) + dv2Txy(i,j) * dt2
             Tyy(i,j) = Tyy1(i,j) + dv2Tyy(i,j) * dt2
          end do
        end do
        !!! K1 = K1 + 2*K2
        dv1z   = dv1z   + 2.d0 * dv2z
        dv1Txx = dv1Txx + 2.d0 * dv2Txx
        dv1Txy = dv1Txy + 2.d0 * dv2Txy
        dv1Tyy = dv1Tyy + 2.d0 * dv2Tyy
        call boundaries_condictions(time + dt2)
        call loop(1)

        !!! Third Runge-Kutta step !!!
        call drv(dv2z, dv2Txx, dv2Txy, dv2Tyy) !!! K3 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
              wz(i,j) =  wz1(i,j) +   dv2z(i,j) * dt 
             Txx(i,j) = Txx1(i,j) + dv2Txx(i,j) * dt 
             Txy(i,j) = Txy1(i,j) + dv2Txy(i,j) * dt 
             Tyy(i,j) = Tyy1(i,j) + dv2Tyy(i,j) * dt 
          end do
        end do
        !!! K1 = k1 + 2*K3
        dv1z   = dv1z   + 2.d0 * dv2z
        dv1Txx = dv1Txx + 2.d0 * dv2Txx
        dv1Txy = dv1Txy + 2.d0 * dv2Txy
        dv1Tyy = dv1Tyy + 2.d0 * dv2Tyy
        call source_term(time + dt)
        call boundaries_condictions(time + dt)
        call loop(1)

        !!! Fourth Runge-Kutta step !!!
        call drv(dv2z, dv2Txx, dv2Txy, dv2Tyy) !!! K4 !!!
        !!! F_i+1 = F_i +  ( [k1 + 2*K2 + 2*K3] + K4 ) *  (dt/6) 
        do j = 2, jmax - 1
          do i = 2, imax - 1
              wz(i,j) =  wz1(i,j) +(  dv1z(i,j) +   dv2z(i,j)) * dt6
             Txx(i,j) = Txx1(i,j) +(dv1Txx(i,j) + dv2Txx(i,j)) * dt6
             Txy(i,j) = Txy1(i,j) +(dv1Txy(i,j) + dv2Txy(i,j)) * dt6
             Tyy(i,j) = Tyy1(i,j) +(dv1Tyy(i,j) + dv2Tyy(i,j)) * dt6
          end do
        end do
        call boundaries_condictions(time + dt)
        call loop(2)

        !!! Print Error !!!
        div     = floor((time+dt)/0.0001)
        errowz  = maxval(abs(wz  -  wz1))/max(maxval(abs( wz)),1.0d0)
        erroTxx = maxval(abs(Txx - Txx1))/max(maxval(abs(Txx)),1.0d0)
        erroTxy = maxval(abs(Txy - Txy1))/max(maxval(abs(Txy)),1.0d0)
        erroTyy = maxval(abs(Tyy - Tyy1))/max(maxval(abs(Tyy)),1.0d0)

!       write(*,*) t, time+dt, errowz, erroTxx, erroTxy, erroTyy
!       if ((div.eq.auxdiv).or.(time.ge.timef)) then
!          auxdiv = auxdiv + 1
!          write(*,*) 'antes do calculo: ',time+dt
!          call exact_solutions(time + dt)
!          errua   = 0.0d0 
!          errva   = 0.0d0 
!          errwza  = 0.0d0 
!          errpsia = 0.0d0 
!          errTxxa = 0.0d0 
!          errTxya = 0.0d0 
!          errTyya = 0.0d0 
!          call normerror_q(errua,errva,errwza,errpsia,
!    & errTxxa,errTxya,errTyya,2)
!          write(7,*) time+dt, errua, errva, errwza, errpsia,
!    & errTxxa, errTxya, errTyya
!          write(10,*) time+dt,errowz,erroTxx,erroTxy,erroTyy
!       end if

        emv = max(errowz,max(erroTxx,max(erroTxy,erroTyy)))
        if ((emv .lt. (1.0d-5)).and.(time .gt. 0.1)) then
           write(*,*) '----------------------------------------' 
           write(*,*) t, time+dt, errowz, erroTxx, erroTxy, erroTyy
           write(*,*) '----------------------------------------' 
           go to 151
        end if
!       call write_solutions_vtk(t)
        t = t + 1
      end do
151   continue
      !!! ********************************************** !!!
      time = time + dt
      call exact_solutions(time)
      call write_error_solutions(time)
      call write_solutions_bin(t)
      call write_solutions_txt(t)
!     call write_solutions_vtk(t)
      !!! ********************************************** !!!
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solver_log_rk4

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j, t, div, auxdiv
      real*8 dt2, dt6, bdfc(imax), 
     &           dv1z(imax,jmax),      dv2z(imax,jmax),
     &            wz1(imax,jmax),      Txx1(imax,jmax), 
     &           Txy1(imax,jmax),      Tyy1(imax,jmax), 
     &       dv1Psixx(imax,jmax),  dv2Psixx(imax,jmax), 
     &       dv1Psixy(imax,jmax),  dv2Psixy(imax,jmax), 
     &       dv1Psiyy(imax,jmax),  dv2Psiyy(imax,jmax),
     &         Psixx1(imax,jmax),    Psixy1(imax,jmax),
     &         Psiyy1(imax,jmax), time, emv, 
     &         errowz, erroTxx, erroTxy, erroTyy 

      common/bd/ bdfc
      
      call init_log(dv1z, dv2z, wz1, Txx1, Txy1, Tyy1, dt2, dt6, 
     &              t, Psixx1, Psixy1, Psiyy1, dv1Psixx, dv2Psixx, 
     &              dv1Psixy, dv2Psixy, dv1Psiyy, dv2Psiyy)

      t      = 1
      auxdiv = 1
      call Psi_to_T

      do while (time .le. timef) 
        time   = dble(t-1) * dt
        wz1    = wz
        Txx1   = Txx
        Txy1   = Txy
        Tyy1   = Tyy
        !!! Psi
        Psixx1 = Psixx
        Psixy1 = Psixy
        Psiyy1 = Psiyy

        !!! First Runge-Kutta step !!!
        call drv_Psi(dv1z, dv1Psixx, dv1Psixy, dv1Psiyy) !!! K1 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
                wz(i,j) =    wz1(i,j) +     dv1z(i,j) * dt2
             Psixx(i,j) = Psixx1(i,j) + dv1Psixx(i,j) * dt2
             Psixy(i,j) = Psixy1(i,j) + dv1Psixy(i,j) * dt2
             Psiyy(i,j) = Psiyy1(i,j) + dv1Psiyy(i,j) * dt2
          end do
        end do
        call source_term(time + dt2)
        call tf_conformation
        call boundaries_condictions(time + dt2)
        call cc_Psi
        call loop(1)

        !!! Second Runge-Kutta step !!!
        call drv_Psi(dv2z, dv2Psixx, dv2Psixy, dv2Psiyy) !!! K2 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
               wz(i,j) =     wz1(i,j) +     dv2z(i,j) * dt2
            Psixx(i,j) =  Psixx1(i,j) + dv2Psixx(i,j) * dt2
            Psixy(i,j) =  Psixy1(i,j) + dv2Psixy(i,j) * dt2
            Psiyy(i,j) =  Psiyy1(i,j) + dv2Psiyy(i,j) * dt2
          end do
        end do

        !!! K1 + 2*K2
        dv1z     =     dv1z + 2.d0 *     dv2z
        dv1Psixx = dv1Psixx + 2.d0 * dv2Psixx
        dv1Psixy = dv1Psixy + 2.d0 * dv2Psixy
        dv1Psiyy = dv1Psiyy + 2.d0 * dv2Psiyy
        call loop(1)

        !!! Third Runge-Kutta step !!! 
        call drv_Psi(dv2z, dv2Psixx, dv2Psixy, dv2Psiyy) !!! K3 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
                  wz(i,j) =      wz1(i,j) +     dv2z(i,j) * dt
               Psixx(i,j) =   Psixx1(i,j) + dv2Psixx(i,j) * dt
               Psixy(i,j) =   Psixy1(i,j) + dv2Psixy(i,j) * dt
               Psiyy(i,j) =   Psiyy1(i,j) + dv2Psiyy(i,j) * dt
          end do
        end do

        !!! K1 + 2*K2 + 2*K3 
            dv1z =     dv1z + 2.d0 *     dv2z
        dv1Psixx = dv1Psixx + 2.d0 * dv2Psixx
        dv1Psixy = dv1Psixy + 2.d0 * dv2Psixy
        dv1Psiyy = dv1Psiyy + 2.d0 * dv2Psiyy
        call source_term(time + dt)
        call tf_conformation
        call boundaries_condictions(time + dt)
        call cc_Psi
        call loop(1)

        !!! Fourth Runge-Kutta step !!!
        call drv_Psi(dv2z, dv2Psixx, dv2Psixy, dv2Psiyy) !!! K4 !!!
        !!! F_i+1 = F_i +  ( [k1 + 2*K2 + 2*K3] + K4 ) *  (dt/6) 
        do j = 2, jmax - 1
          do i = 2, imax - 1
               wz(i,j) =    wz1(i,j) + dt6 * (dv1z(i,j) + dv2z(i,j))
            Psixx(i,j) = Psixx1(i,j) + dt6 * (dv1Psixx(i,j) +
     &                 dv2Psixx(i,j))
            Psixy(i,j) = Psixy1(i,j) + dt6 * (dv1Psixy(i,j) +
     &                 dv2Psixy(i,j))
            Psiyy(i,j) = Psiyy1(i,j) + dt6 * (dv1Psiyy(i,j) +
     &                 dv2Psiyy(i,j))
          end do
        end do
        call loop(2)

        call Psi_to_T

        !!! Print Error !!!
        div = floor((time+dt)/0.0001)
        errowz  = maxval(abs(wz-wz1))/max(maxval(abs(wz)),1.0d0)
        erroTxx = maxval(abs(Txx-Txx1))/max(maxval(abs(Txx)),1.0d0)
        erroTxy = maxval(abs(Txy-Txy1))/max(maxval(abs(Txy)),1.0d0)
        erroTyy = maxval(abs(Tyy-Tyy1))/max(maxval(abs(Tyy)),1.0d0)
        
        div = floor((time+dt)/0.001)
!       if ((div.eq.auxdiv).or.(time.ge.timef)) then
!          auxdiv = auxdiv + 1
!          call Psi_to_T
!          errou   = maxval(abs(ux-uxe))
!          errov   = maxval(abs(uy-uye))
!          errowz  = maxval(abs(wz-wz1))
!          erropsi = maxval(abs(psi-psie))
!          erroTxx = maxval(abs(Txx-Txx1))
!          erroTxy = maxval(abs(Txy-Txy1))
!          erroTyy = maxval(abs(Tyy-Tyy1))
!          call normerror(erru,errv,errwz,errpsi,errTxx,errTxy,errTyy)
!          write(1,*) time+dt, errou, errov, errowz, erropsi,
!    & erroTxx, erroTxy, erroTyy
!          write(1,*) time+dt, errowz, erroTxx, erroTxy, erroTyy
!          write(*,*) 'Max =>', time+dt, errou, errov, errowz, erropsi,
!    & erroTxx, erroTxy, erroTyy
!          write(*,*) 'Max =>', time+dt,errowz,erroTxx,erroTxy,erroTyy
!          write(2,*) time+dt, erru, errv, errwz, errpsi,
!    & errTxx, errTxy, errTyy
!          write(*,*) 'Int =>',time+dt, erru, errv, errwz, errpsi,
!    & errTxx, errTxy, errTyy
!          write(*,*) '*****************************************',
!    &'*********************************************************',
!    &'*********************************************************',
!    &'*********************************************************'
!       end if

        emv = max(errowz,max(erroTxx,max(erroTxy,erroTyy)))
        if ((emv .lt. 1.0d-7).and.(time .gt. 0.01)) then
!          write(1,*) '----------------------------------------' 
!          write(1,*) time+dt, errowz, erroTxx, erroTxy, erroTyy
!          write(1,*) '----------------------------------------' 
           go to 151
        end if
        t = t + 1
      end do
151   continue
      !!! ********************************************** !!!
      time = time + dt
      call exact_solutions(time)
      call write_error_solutions(time)
      !call write_solutions_bin(t)
      call write_solutions_vtk(t)
      return
      end
      
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>!Starts the simulation process after verifying the number of points in 
!!!the x and y directions. After it initializes the variables, the solver 
!!!function computes the Runge-Kutta method and verifies the convergence 
!!!criteria for vorticity. At the end, this function writes the results.
      subroutine solver_euler

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer t, div, auxdiv
      real*8 dt2, dt6, bdfc(imax), dv1z(imax,jmax), dv2z(imax,jmax),
     &       wz1(imax,jmax), Txx1(imax,jmax), Txy1(imax,jmax),
     &       Tyy1(imax,jmax), emv,
     &       dv1Txx(imax,jmax), dv2Txx(imax,jmax), 
     &       dv1Txy(imax,jmax), dv2Txy(imax,jmax),
     &       dv1Tyy(imax,jmax), dv2Tyy(imax,jmax), time,
     &       errowz , erroTxx, erroTxy, erroTyy
      common/bd/ bdfc
      
      time   = 0.d0

      call init(dv1z, dv2z, dv1Txx, dv2Txx, dv1Txy, dv2Txy, 
     & dv1Tyy, dv2Tyy, wz1, Txx1, Txy1, Tyy1, dt2, dt6, t)
      
      call boundaries_condictions(time)
      call source_term(time)
      
      t        = 1
      auxdiv   = 1
      do while (time .lt. timef-dt)
        time = dble(t - 1) * dt

         wz1 = wz
        Txx1 = Txx
        Txy1 = Txy
        Tyy1 = Tyy

        !!! Euler step !!!
        call drv(dv1z, dv1Txx, dv1Txy, dv1Tyy) 

         wz =  wz1 +   dv1z * dt
        Txx = Txx1 + dv1Txx * dt 
        Txy = Txy1 + dv1Txy * dt 
        Tyy = Tyy1 + dv1Tyy * dt 

        call source_term(time + dt)
        call boundaries_condictions(time + dt)
        call loop(2)

        !!! Print Error !!!
        div = floor((time+dt)/0.0001)
        errowz  = maxval(abs(wz-wz1))/max(maxval(abs(wz)),1.0d0)
        erroTxx = maxval(abs(Txx-Txx1))/max(maxval(abs(Txx)),1.0d0)
        erroTxy = maxval(abs(Txy-Txy1))/max(maxval(abs(Txy)),1.0d0)
        erroTyy = maxval(abs(Tyy-Tyy1))/max(maxval(abs(Tyy)),1.0d0)

c       if ((div.eq.auxdiv).or.(time.ge.timef)) then
c          auxdiv = auxdiv + 1
c          errou   = maxval(abs(ux-uxe))
c          errov   = maxval(abs(uy-uye))
c          errowz  = maxval(abs(wz-wz1))
c          erropsi = maxval(abs(psi-psie))
c          erroTxx = maxval(abs(Txx-Txx1))
c          erroTxy = maxval(abs(Txy-Txy1))
c          erroTyy = maxval(abs(Tyy-Tyy1))
c          call normerror(erru,errv,errwz,errpsi,errTxx,errTxy,errTyy)
c          write(*,*) 'Max =>', time+dt, errou, errov, errowz, erropsi,
c    & erroTxx, erroTxy, erroTyy
           write(*,*) 'Max =>', time+dt,errowz,erroTxx,erroTxy,erroTyy
c          write(*,*) 'Int =>',time+dt, erru, errv, errwz, errpsi,
c    & errTxx, errTxy, errTyy
c          write(*,*) '*****************************************',
c    &'*********************************************************',
c    &'*********************************************************',
c    &'*********************************************************'
c       end if

        emv = max(errowz,max(erroTxx,max(erroTxy,erroTyy)))
        if ((emv .lt. 1.0d-7).and.(time .gt. 0.01)) then
!          write(1,*) '----------------------------------------' 
!          write(1,*) time+dt, errowz, erroTxx, erroTxy, erroTyy
!          write(1,*) '----------------------------------------' 
           go to 151
        end if
        t = t + 1
      end do
151   continue
      !!! ********************************************** !!!
      time = time + dt
      call exact_solutions(time)
      call write_error_solutions(time)
      !call write_solutions_bin(t)
      call write_solutions_vtk(t)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solver_log_euler

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer i, j, t, auxdiv
      real*8 dt2, dt6, bdfc(imax), 
     &           dv1z(imax,jmax),      dv2z(imax,jmax),
     &            wz1(imax,jmax),      Txx1(imax,jmax), 
     &           Txy1(imax,jmax),      Tyy1(imax,jmax), 
     &       dv1Psixx(imax,jmax),  dv2Psixx(imax,jmax), 
     &       dv1Psixy(imax,jmax),  dv2Psixy(imax,jmax), 
     &       dv1Psiyy(imax,jmax),  dv2Psiyy(imax,jmax),
     &         Psixx1(imax,jmax),    Psixy1(imax,jmax),
     &         Psiyy1(imax,jmax), time, emv, 
     &         errowz, erroTxx, erroTxy, erroTyy 



      common/bd/ bdfc
      
      call init_log(dv1z, dv2z, wz1, Txx1, Txy1, Tyy1, dt2, dt6, 
     &              t, Psixx1, Psixy1, Psiyy1, dv1Psixx, dv2Psixx, 
     &              dv1Psixy, dv2Psixy, dv1Psiyy, dv2Psiyy)

      t      = 1
      auxdiv = 1
      call Psi_to_T

      do while (time .le. timef) 
        time   = dble(t-1) * dt
        wz1    = wz
        Txx1   = Txx
        Txy1   = Txy
        Tyy1   = Tyy
        !!! Psi
        Psixx1 = Psixx
        Psixy1 = Psixy
        Psiyy1 = Psiyy

        !!! First Runge-Kutta step !!!
        call drv_Psi(dv1z, dv1Psixx, dv1Psixy, dv1Psiyy) !!! K1 !!!
        do j = 2, jmax - 1
          do i = 2, imax - 1
                wz(i,j) =    wz1(i,j) +     dv1z(i,j) * dt
             Psixx(i,j) = Psixx1(i,j) + dv1Psixx(i,j) * dt
             Psixy(i,j) = Psixy1(i,j) + dv1Psixy(i,j) * dt
             Psiyy(i,j) = Psiyy1(i,j) + dv1Psiyy(i,j) * dt
          end do
        end do
        call source_term(time + dt)
        call tf_conformation
        call boundaries_condictions(time + dt)
        call cc_Psi
        call loop(2)

        call Psi_to_T
        !!! Print Error !!!
        errowz  = maxval(abs(wz-wz1))/max(maxval(abs(wz)),1.0d0)
        erroTxx = maxval(abs(Txx-Txx1))/max(maxval(abs(Txx)),1.0d0)
        erroTxy = maxval(abs(Txy-Txy1))/max(maxval(abs(Txy)),1.0d0)
        erroTyy = maxval(abs(Tyy-Tyy1))/max(maxval(abs(Tyy)),1.0d0)
        
!       div = floor((time+dt)/0.001)
!       if ((div.eq.auxdiv).or.(time.ge.timef)) then
!          auxdiv = auxdiv + 1
!          call Psi_to_T
!          errou   = maxval(abs(ux-uxe))
!          errov   = maxval(abs(uy-uye))
!          errowz  = maxval(abs(wz-wz1))
!          erropsi = maxval(abs(psi-psie))
!          erroTxx = maxval(abs(Txx-Txx1))
!          erroTxy = maxval(abs(Txy-Txy1))
!          erroTyy = maxval(abs(Tyy-Tyy1))
!          call normerror(erru,errv,errwz,errpsi,errTxx,errTxy,errTyy)
!          write(1,*) time+dt, errou, errov, errowz, erropsi,
!    & erroTxx, erroTxy, erroTyy
!          write(1,*) time+dt, errowz, erroTxx, erroTxy, erroTyy
!          write(*,*) 'Max =>', time+dt, errou, errov, errowz, erropsi,
!    & erroTxx, erroTxy, erroTyy
!          write(*,*) 'Max =>', time+dt,errowz,erroTxx,erroTxy,erroTyy
!          write(2,*) time+dt, erru, errv, errwz, errpsi,
!    & errTxx, errTxy, errTyy
!          write(*,*) 'Int =>',time+dt, erru, errv, errwz, errpsi,
!    & errTxx, errTxy, errTyy
!          write(*,*) '*****************************************',
!    &'*********************************************************',
!    &'*********************************************************',
!    &'*********************************************************'
!       end if

        emv = max(errowz,max(erroTxx,max(erroTxy,erroTyy)))
        if ((emv .lt. 1.0d-7).and.(time .gt. 0.01)) then
!          write(1,*) '----------------------------------------' 
!          write(1,*) time+dt, errowz, erroTxx, erroTxy, erroTyy
!          write(1,*) '----------------------------------------' 
           go to 151
        end if
        t = t + 1
      end do
151   continue
      !!! ********************************************** !!!
      time = time + dt
      call exact_solutions(time)
      !call write_error_solutions(time)
      !call write_solutions_bin(t)
      !call write_solutions_vtk(t)

      return
      end

