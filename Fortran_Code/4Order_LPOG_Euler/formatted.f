ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutines transforms binary data in asc       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> Convert the binary output data.bin in formatted form spatial.dat
      program formated

      implicit none
      include 'par.nn'
      integer i, j, t 
      real*8 x, y,
     &      uxb(jmax), wzb(jmax), Txxb(jmax), Txyb(jmax), Tyyb(jmax),
     &       ux(imax,jmax),  uy(imax,jmax),  wz(imax,jmax),
     &      Txx(imax,jmax), Txy(imax,jmax), Tyy(imax,jmax)

      open(2,file='data.bin',form='unformatted')
      read(2) t
      read(2) ux, uy, wz, Txx, Txy, Tyy
      close (unit=2)

      ! writes data to spacial space to be open by tecplot
      open (3, file = 'spatial.dat',status = 'unknown')
      write(3,*) 'VARIABLES="x","y","velu","vely","vortz",
     &"Txx","Txy","Tyy"'
      write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

      do j = 1, jmax
        y       = dble(j-1) * dy
c       uxb(j)  = - y * y + 2.d0 * y
c       wzb(j)  = - 2.d0 * y + 2.d0
c       Txxb(j) = 8.d0 * (y-1.d0)*(y-1.d0) * Wi * (1.d0-beta) / Re
c       Txyb(j) = 2.d0 * (y-1.d0) * (beta-1.d0) / Re
c       Tyyb(j) = 0.d0
        do i = 1, imax
          x = dble(i-1) * dx
c         write(3,5) x, y, ux(i,j)-uxb(j), uy(i,j), 
c    &               wz(i,j)-wzb(j),
c    &               Txx(i,j) - Txxb(j), 
c    &               Txy(i,j) - Txyb(j), 
c    &               Tyy(i,j) - Tyyb(j)
          write(3,5) x, y, ux(i,j), uy(i,j), 
     &               wz(i,j), Txx(i,j), Txy(i,j), Tyy(i,j)
        end do
        write(3,5)
      end do

      close (unit=3)

    5 format(1x,2d14.6,6d17.9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   end subroutines transforms binary data in asc       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
