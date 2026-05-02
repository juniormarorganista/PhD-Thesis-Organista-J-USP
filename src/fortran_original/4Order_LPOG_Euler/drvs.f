cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv(dvz, dvTxx, dvTxy, dvTyy)

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 d2wzdx2(imax,jmax),   d2wzdy2(imax,jmax),
     &           uwz(imax,jmax),       vwz(imax,jmax),
     &        duwzdx(imax,jmax),    dvwzdy(imax,jmax),  
     &           dvz(imax,jmax),     dvTxx(imax,jmax),
     &         dvTxy(imax,jmax),     dvTyy(imax,jmax),
     &          uTxx(imax,jmax),      vTxx(imax,jmax),
     &          uTxy(imax,jmax),      vTxy(imax,jmax),
     &          uTyy(imax,jmax),      vTyy(imax,jmax), 
     &     d2Txxdxdy(imax,jmax), d2Tyydxdy(imax,jmax),
     &      d2Txydy2(imax,jmax),  d2Txydx2(imax,jmax),
     &       duTxxdx(imax,jmax),   dvTxxdy(imax,jmax),
     &       duTxydx(imax,jmax),   dvTxydy(imax,jmax),
     &       duTyydx(imax,jmax),   dvTyydy(imax,jmax),
     &          dudx(imax,jmax),      dudy(imax,jmax),
     &          dvdx(imax,jmax),      dvdy(imax,jmax),
     &           aux(imax,jmax),      aux2(imax,jmax),
     &          FTij(imax,jmax)

      call nlterms(uwz, vwz, uTxx, vTxx, uTxy, vTxy, uTyy, vTyy)

!!!!!!#################################
!!!!!!drv_wz
!!!!!!#################################
      ! derivative calculations
      call derx(duwzdx,uwz)
      call derxx(d2wzdx2,wz)
      call dery(dvwzdy,vwz)
      call deryy(d2wzdy2,wz)

      call derx(aux, Txx)
      call dery(d2Txxdxdy, aux)
      call deryy(d2Txydy2, Txy)
      call derxx(d2Txydx2, Txy)
      call derx(aux2, Tyy)
      call dery(d2Tyydxdy, aux2)

      dvz = - duwzdx - dvwzdy + betann*(d2wzdx2 + d2wzdy2) / Rey
     & + d2Txxdxdy +  d2Txydy2 - d2Txydx2 - d2Tyydxdy - tfwz

!!!!!!#################################
!!!!!!drv_T
!!!!!!#################################
      call derx(duTxxdx, uTxx)
      call dery(dvTxxdy, vTxx)
      call derx(duTxydx, uTxy)
      call dery(dvTxydy, vTxy)
      call derx(duTyydx, uTyy)
      call dery(dvTyydy, vTyy)

      call derx(dudx, ux)
      call dery(dudy, ux)
      call derx(dvdx, uy)
      call dery(dvdy, uy)

      FTij  = 1.0d0 + (epsylon*Rey*Wi/(1.0d0 - betann)) * 
     & (Txx + Tyy)

      dvTxx = - FTij*Txx/Wi - duTxxdx - dvTxxdy + 
     & 2.d0*Txx*dudx + 2.d0*Txy*dudy - 
     & 2.0d0*xi*Txx*dudx - xi*Txy*(dudy + dvdx) -
     & (alphaG*Rey/(1.0d0 - betann))*(Txx*Txx + Txy*Txy) +
     & (2.0d0*(1.0d0 - betann)/(Rey*Wi))*dudx + tfTxx/Wi

      dvTxy = - FTij*Txy/Wi - duTxydx - dvTxydy + 
     & Txx*dvdx + Tyy*dudy - 
     & xi*Txy*dudx - 0.5d0*xi*Tyy*(dudy + dvdx) -
     & 0.5d0*xi*Txx*(dudy + dvdx) - xi*Txy*dvdy -
     & (alphaG*Rey/(1.0d0 - betann))*(Txx*Txy + Txy*Tyy) +
     & ((1.0d0 - betann)/(Rey*Wi))*(dudy + dvdx) + tfTxy/Wi

      dvTyy = - FTij*Tyy/Wi - duTyydx - dvTyydy + 
     & 2.d0*Txy*dvdx + 2.d0*Tyy*dvdy - 
     & xi*Txy*(dvdx + dudy) - 2.0d0*xi*Tyy*dvdy -
     & (alphaG*Rey/(1.0d0 - betann))*(Txy*Txy + Tyy*Tyy) +
     & (2.0d0*(1.0d0 - betann)/(Rey*Wi))*dvdy + tfTyy/Wi

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nlterms(uwz, vwz, uTxx, vTxx, uTxy, vTxy, uTyy, vTyy)

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8  uwz(imax,jmax),   vwz(imax,jmax),
     &       uTxx(imax,jmax),  vTxx(imax,jmax),
     &       uTxy(imax,jmax),  vTxy(imax,jmax),
     &       uTyy(imax,jmax),  vTyy(imax,jmax)

      uwz  = ux * wz
      vwz  = uy * wz
      uTxx = ux * Txx
      vTxx = uy * Txx
      uTxy = ux * Txy
      vTxy = uy * Txy
      uTyy = ux * Tyy
      vTyy = uy * Tyy

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv_Psi(dvz, dvPsixx, dvPsixy, dvPsiyy)

      implicit none
      include 'par.nn'
      include 'comm.var'

      real*8 Omegalog(imax,jmax),
     &      duPsixxdx(imax,jmax), duPsixydx(imax,jmax), 
     &      duPsiyydx(imax,jmax), dvPsixxdy(imax,jmax),
     &      dvPsixydy(imax,jmax), dvPsiyydy(imax,jmax),
     &            Bxx(imax,jmax),       Byy(imax,jmax),
     &            Bxy(imax,jmax),       Mxx(imax,jmax), 
     &            Mxy(imax,jmax),       Myy(imax,jmax),
     &        dvPsixx(imax,jmax),   dvPsixy(imax,jmax),
     &        dvPsiyy(imax,jmax),       dvz(imax,jmax),
     &         uPsixx(imax,jmax),    uPsixy(imax,jmax),
     &         uPsiyy(imax,jmax),    vPsixx(imax,jmax),
     &         vPsixy(imax,jmax),    vPsiyy(imax,jmax),
     &           dudx(imax,jmax),      dudy(imax,jmax),
     &        d2wzdx2(imax,jmax),   d2wzdy2(imax,jmax),
     &            uwz(imax,jmax),       vwz(imax,jmax),
     &         duwzdx(imax,jmax),    dvwzdy(imax,jmax),  
     &      d2Txxdxdy(imax,jmax), d2Tyydxdy(imax,jmax),
     &       d2Txydy2(imax,jmax),  d2Txydx2(imax,jmax),
     &            aux(imax,jmax)

      call Psi_to_T

      call kernel_conformation(Omegalog,Bxx,Bxy,Byy,Mxx,Mxy,Myy)

      call nlterms_logPsi(uwz, vwz, uPsixx, vPsixx, uPsixy, 
     &                    vPsixy, uPsiyy, vPsiyy)

!!!!!!#################################
!!!!!!drv_wz
!!!!!!#################################
      !!! Derivative calculations
      call derx(duwzdx,uwz)
      call derxx(d2wzdx2,wz)
      call dery(dvwzdy,vwz)
      call deryy(d2wzdy2,wz)

      call derx(aux, Txx)
      call dery(d2Txxdxdy, aux)

      call derx(aux, Tyy)
      call dery(d2Tyydxdy, aux)
      call deryy(d2Txydy2, Txy)
      call derxx(d2Txydx2, Txy)


      dvz = - duwzdx - dvwzdy + betann*(d2wzdx2 + d2wzdy2)/Rey+
     & ( d2Txxdxdy +  d2Txydy2 - d2Txydx2 - d2Tyydxdy ) -
     &  tfwz

!!!!!!#################################
!!!!!!drv_PSI
!!!!!!#################################
      call derx(duPsixxdx, uPsixx)
      call dery(dvPsixxdy, vPsixx)

      call derx(duPsixydx, uPsixy)
      call dery(dvPsixydy, vPsixy)

      call derx(duPsiyydx, uPsiyy)
      call dery(dvPsiyydy, vPsiyy)

      call derx(dudx, ux)
      call dery(dudy, ux)

      dvPsixx = - duPsixxdx - dvPsixxdy + 
     & 2.d0 * Omegalog * Psixy + 2.d0 * Bxx + Mxx / Wi -
     & tfpsixx

      dvPsixy = - duPsixydx - dvPsixydy + 
     & Omegalog * (Psiyy - Psixx) + 2.d0 * Bxy + Mxy / Wi -
     & tfpsixy

      dvPsiyy = - duPsiyydx - dvPsiyydy - 
     & 2.d0 * Omegalog * Psixy + 2.d0 * Byy + Myy / Wi -
     & tfpsiyy

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nlterms_logPsi(uwz, vwz, uPsixx, vPsixx, uPsixy, 
     &                       vPsixy, uPsiyy, vPsiyy)

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8  uPsixx(imax,jmax),   vPsixx(imax,jmax),
     &        uPsixy(imax,jmax),   vPsixy(imax,jmax),
     &        uPsiyy(imax,jmax),   vPsiyy(imax,jmax),
     &           uwz(imax,jmax),      vwz(imax,jmax) 

      uwz    = ux * wz
      vwz    = uy * wz
      uPsixx = ux * Psixx
      vPsixx = uy * Psixx
      uPsixy = ux * Psixy
      vPsixy = uy * Psixy
      uPsiyy = ux * Psiyy
      vPsiyy = uy * Psiyy

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of main subroutines                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
