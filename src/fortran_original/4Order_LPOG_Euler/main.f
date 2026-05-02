cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                main subroutines                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Verifies the number of adopted points in x and y directions and 
!! proceeds the simulation calling the solver function.
      program main

      implicit none
      include 'par.nn'
      
      if ((cfl.lt.dt).or.(visc.lt.dt)) then
         write(*,*) 'Delta t deve ser refinado!!!'
         stop
      endif

      ! verification of adopted points in x and y directions
      if ( ty_sim .eq. 1  ) call solver_rk4
      if ( ty_sim .eq. 12 ) call solver_euler
      if ( ty_sim .eq. 2  ) call solver_log_rk4
      if ( ty_sim .eq. 22 ) call solver_log_euler

      end program

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Initializes the program variables, gives the values of the boundary 
!! layer profile for each node, mounts the LU matrix for the computational 
!! filter and the left hand side of the derivative calculation and calls 
!! the x function of the disturbance strip for the Tollmien-Schlichting disturbance.
      subroutine init(dv1z, dv2z, dv1Txx, dv2Txx, dv1Txy, dv2Txy, 
     & dv1Tyy, dv2Tyy, wz1, Txx1, Txy1, Tyy1, dt2, dt6, t)
      
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      integer t
      real*8 dt2, dt6 , 
     &       dv1z(imax,jmax),   dv2z(imax,jmax),  wz1(imax,jmax), 
     &       Txx1(imax,jmax),   Txy1(imax,jmax), Tyy1(imax,jmax), 
     &     dv1Txx(imax,jmax), dv2Txx(imax,jmax), 
     &     dv1Txy(imax,jmax), dv2Txy(imax,jmax), 
     &     dv1Tyy(imax,jmax), dv2Tyy(imax,jmax),
     &      afilx(imax), bfilx(imax), cfilx(imax),
     &      afily(jmax), bfily(jmax), cfily(jmax)
      common/filtx/ afilx, bfilx, cfilx
      common/filty/ afily, bfily, cfily

      t   = 1
      dt2 = 0.5d0 * dt
      dt6 = dt / 6.d0
      
      ! mounts the lu matrix for the filter
      call lhs_tridfx(afilx, bfilx, cfilx)
      call lhs_tridfy(afily, bfily, cfily)
      
      ! mounts the lhs for the derivative calculation
      call derivs_k
      
      ! initialize all variables
      wz1    = 0.d0
      Txx1   = 0.d0
      Txy1   = 0.d0
      Tyy1   = 0.d0
      dv1z   = 0.d0
      dv2z   = 0.d0
      dv1Txx = 0.d0
      dv2Txx = 0.d0
      dv1Txy = 0.d0
      dv2Txy = 0.d0
      dv1Tyy = 0.d0
      dv2Tyy = 0.d0
       wz    = 0.d0
       ux    = 0.d0
       uy    = 0.d0
      psi    = 0.d0
      Txx    = 0.d0
      Txy    = 0.d0
      Tyy    = 0.d0
     
      call exact_solutions(0.d0)

c      wz = 0.25 *  wze
c      ux = 0.25 *  uxe
c      uy = 0.25 *  uye
c     psi = 0.25 * psie
c     Txx = 0.25 * Txxe
c     Txy = 0.25 * Txye
c     Tyy = 0.25 * Tyye

       wz =  wze
       ux =  uxe
       uy =  uye
      psi = psie
      Txx = Txxe
      Txy = Txye
      Tyy = Tyye     

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_log(dv1z, dv2z, wz1, Txx1, Txy1, Tyy1, dt2, dt6, 
     & t, Psixx1, Psixy1, Psiyy1, dv1Psixx, dv2Psixx, 
     & dv1Psixy, dv2Psixy, dv1Psiyy, dv2Psiyy)

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      real*8 dt2, dt6 , t,
     &       dv1z(imax,jmax), dv2z(imax,jmax),  wz1(imax,jmax), 
     &       Txx1(imax,jmax), Txy1(imax,jmax), Tyy1(imax,jmax), 
     &       afilx(imax), bfilx(imax), cfilx(imax),
     &       afily(jmax), bfily(jmax), cfily(jmax),
     &       Psixx1(imax,jmax), Psixy1(imax,jmax),
     &       Psiyy1(imax,jmax), 
     &     dv1Psixx(imax,jmax), dv2Psixx(imax,jmax), 
     &     dv1Psixy(imax,jmax), dv2Psixy(imax,jmax), 
     &     dv1Psiyy(imax,jmax), dv2Psiyy(imax,jmax)
      common/filtx/ afilx, bfilx, cfilx
      common/filty/ afily, bfily, cfily

      t   = 1
      dt2 = 0.5d0 * dt
      dt6 = dt / 6.d0

      ! mounts the lu matrix for the filter
      call lhs_tridfx(afilx, bfilx, cfilx)
      call lhs_tridfy(afily, bfily, cfily)

      ! mounts the lhs for the derivative calculation
      call derivs_k

      ! initialize all variables
      wz1      = 0.d0
      Txx1     = 0.d0
      Txy1     = 0.d0
      Tyy1     = 0.d0
      dv1z     = 0.d0
      dv2z     = 0.d0
      Psixx1   = 0.d0
      Psixy1   = 0.d0
      Psiyy1   = 0.d0
      dv1Psixx = 0.d0
      dv2Psixx = 0.d0
      dv1Psixy = 0.d0
      dv2Psixy = 0.d0
      dv1Psiyy = 0.d0
      dv2Psiyy = 0.d0

      call source_term(0.d0)
      call tf_conformation

      call exact_solutions(0.d0)
       ux =  uxe
       uy =  uye
      psi = psie
       wz =  wze
      Txx = Txxe
      Txy = Txye
      Tyy = Tyye
      call T_to_Psi

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Actualizes the variables for each Runge-Kutta step
      subroutine loop(num)

      ! actualize the variables for each RK-step
      implicit none
      include 'par.nn'
      include 'comm.var'
      integer num, r
      
      r = num
      call filter_trid_x
      call filter_trid_y
      call poi_psi4a
      call calc_ux
      call calc_uy
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   end of main subroutines                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
