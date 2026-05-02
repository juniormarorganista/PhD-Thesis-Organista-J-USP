ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 subroutines writes the results        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_error_solutions(t)
 
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      character(200) nome1 , nome2, nome3
      real*8 r, t
      real*8 eu0, ev0, ewz0, epsi0, eTxx0, eTxy0, eTyy0
      real*8 eu1, ev1, ewz1, epsi1, eTxx1, eTxy1, eTyy1
      real*8 eu2, ev2, ewz2, epsi2, eTxx2, eTxy2, eTyy2
       
       r = t
       write(nome1,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a)')
     & 'Euler_err_inf_imax_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',ty_mms,'.dat'

       write(nome2,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a)')
     & 'Euler_err_1st_imax_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',ty_mms,'.dat'

       write(nome3,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a)')
     & 'Euler_err_2nd_imax_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',ty_mms,'.dat'
      
      eu0   = maxval(abs(ux  - uxe ))
      ev0   = maxval(abs(uy  - uye ))
      ewz0  = maxval(abs(wz  - wze ))
      epsi0 = maxval(abs(psi - psie))
      eTxx0 = maxval(abs(Txx - Txxe))
      eTxy0 = maxval(abs(Txy - Txye))
      eTyy0 = maxval(abs(Tyy - Tyye))
      
      call normerror_q(eu1,ev1,ewz1,epsi1,eTxx1,eTxy1,eTyy1,1)
      call normerror_q(eu2,ev2,ewz2,epsi2,eTxx2,eTxy2,eTyy2,2)

      open (1, file = nome1)
      open (2, file = nome2)
      open (3, file = nome3)
            write(1,*) eu0, ev0, ewz0, epsi0, eTxx0, eTxy0, eTyy0
            write(2,*) eu1, ev1, ewz1, epsi1, eTxx1, eTxy1, eTyy1
            write(3,*) eu2, ev2, ewz2, epsi2, eTxx2, eTxy2, eTyy2
      close (unit=1)
      close (unit=2)
      close (unit=3)
       
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_solutions_txt(ind)
      
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
       
      integer i, j, ind
      real*8 x, y
      
      character(200) nome1
      character(200) nome2
 
       write(nome1,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a,i0.6,a)')
     & 'sol_num_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',
     & ty_mms,'_ind_',ind,'.txt'
       
       write(nome2,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a,i0.6,a)')
     & 'sol_exa_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',
     & ty_mms,'_ind_',ind,'.txt'
      
      open(unit=1, file=nome1)
      open(unit=2, file=nome2)
      do i=1,imax
         x = dble(i-1)*dx
         do j=1,jmax
            y = dble(j-1)*dy
            write(1,*) i, j, x, y, ux(i,j) , uy(i,j) , wz(i,j) , 
     & psi(i,j), Txx(i,j), Txy(i,j), Tyy(i,j)
            write(2,*) i, j, x, y, uxe(i,j), uye(i,j), wze(i,j), 
     & psie(i,j), Txxe(i,j), Txye(i,j), Tyye(i,j)
         end do
      end do
      close(unit=1)
      close(unit=2)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_solutions_bin(ind)
      
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
       
      integer ind
      
      character(200) nome1
      character(200) nome2
 
       write(nome1,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a,i0.6,a)')
     & 'sol_num_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',
     & ty_mms,'_ind_',ind,'.bin'
       
       write(nome2,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,a,i0.2,a,i0.6,a)')
     & 'sol_exa_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,'_MMS_',
     & ty_mms,'_ind_',ind,'.bin'
      
      open(unit=1, file=nome1, form='unformatted')
      open(unit=2, file=nome2, form='unformatted')

      write(1) ux , uy , wz , psi , Txx , Txy , Tyy
      write(2) uxe, uye, wze, psie, Txxe, Txye, Tyye

      close(unit=1)
      close(unit=2)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_solutions_vtk(ind)
      
      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
       
      integer i, j, dim_z, ind, cont
      real*8 x(imax*jmax), y(imax*jmax)
      
      character(200) nome1, nome2, nome3, nome4, nome5
      
      dim_z = 1
      cont  = 1
      do j = 1, jmax
        do i = 1, imax
          x(cont) = x0 + dble(i-1)*dx
          y(cont) = y0 + dble(j-1)*dy
          cont    = cont + 1
        end do
      end do

       !!! Name Load
       write(nome1,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,
     & a,i0.2,a,i0.4,a)')
     & 'Euler_sol_num_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,
     & '_MMS_',ty_mms,'_step_',ind,'.bin'
       
       write(nome2,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,
     & a,i0.2,a,i0.4,a)')
     & 'Euler_sol_exa_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,
     & '_MMS_',ty_mms,'_step_',ind,'.bin'


       !!! Name vtk
       write(nome3,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,
     & a,i0.2,a,i0.4,a)')
     & 'Euler_sol_num_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,
     & '_MMS_',ty_mms,'_step_',ind,'.vtk'
       
       write(nome4,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,
     & a,i0.2,a,i0.4,a)')
     & 'Euler_sol_exa_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,
     & '_MMS_',ty_mms,'_step_',ind,'.vtk'

       write(nome5,'(a,i0.4,a,f0.8,a,f0.8,a,
     & f0.8,a,f0.8,a,f0.8,a,f0.8,
     & a,f0.8,a,f0.8,a,i0.2,
     & a,i0.2,a,i0.4,a)')
     & 'Euler_sol_comp_',imax,'_Re_',Rey,'_Betann_',betann,'_Wi_',
     & wi,'_epsilon_',epsylon,'_xi_',xi,'_alphaG_',alphaG,
     & '_Dt_',dt,'_at_',at,'_tipsim_',ty_sim,
     & '_MMS_',ty_mms,'_step_',ind,'.vtk'

      !!! Write file vtk format
      open(unit=3, file=nome3, form='formatted')
      write(3,9)'# vtk DataFile Version 2.0'
      write(3,9)'vtk output'
      write(3,9)'ASCII'
      write(3,9)'DATASET STRUCTURED_GRID'
      write(3,7)'DIMENSIONS',imax, jmax, dim_z
      write(3,8)'POINTS',imax*jmax,'double'
      write(3,9)''
     
      do i = 1, imax*jmax
        write(3,5) x(i), y(i), 0.0
      end do
     
      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS velocidade-u float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) ux(i,j)
        end do 
      end do
     
      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS velocidade-v float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) uy(i,j)
        end do 
      end do
     
      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS vorticidade-z float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) wz(i,j)
        end do 
      end do
     
      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS streamwise float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) psi(i,j)
        end do 
      end do
     
      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS Txx float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) Txx(i,j)
        end do 
      end do

      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS Txy float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) Txy(i,j)
        end do 
      end do

      write(3,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(3,10)'POINT_DATA', imax*jmax
      write(3,*)''
      write(3,9)'SCALARS Tyy float 1'
      write(3,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(3,6) Tyy(i,j)
        end do 
      end do
      close(unit=3)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(unit=4, file=nome4, form='formatted')
      write(4,9)'# vtk DataFile Version 2.0'
      write(4,9)'vtk output'
      write(4,9)'ASCII'
      write(4,9)'DATASET STRUCTURED_GRID'
      write(4,7)'DIMENSIONS',imax, jmax, dim_z
      write(4,8)'POINTS',imax*jmax,'double'
      write(4,9)''
     
      do i = 1, imax*jmax
        write(4,5) x(i), y(i), 0.0
      end do
     
      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS velocidade-ue float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) uxe(i,j)
        end do 
      end do
     
      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS velocidade-ve float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) uye(i,j)
        end do 
      end do
     
      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS vorticidade-ze float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) wze(i,j)
        end do 
      end do
     
      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS streamwisee float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) psie(i,j)
        end do 
      end do
     
      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS Txxe float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) Txxe(i,j)
        end do 
      end do

      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS Txye float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) Txye(i,j)
        end do 
      end do

      write(4,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(4,10)'POINT_DATA', imax*jmax
      write(4,*)''
      write(4,9)'SCALARS Tyye float 1'
      write(4,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(4,6) Tyye(i,j)
        end do 
      end do
      close(unit=4)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(unit=5, file=nome5, form='formatted')
      write(5,9)'# vtk DataFile Version 2.0'
      write(5,9)'vtk output'
      write(5,9)'ASCII'
      write(5,9)'DATASET STRUCTURED_GRID'
      write(5,7)'DIMENSIONS',imax, jmax, dim_z
      write(5,8)'POINTS',imax*jmax,'double'
      write(5,9)''
     
      do i = 1, imax*jmax
        write(5,5) x(i), y(i), 0.0
      end do
     
      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS velocidade-uc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(ux(i,j) - uxe(i,j))
        end do 
      end do
     
      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS velocidade-vc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(uy(i,j) - uye(i,j))
        end do 
      end do
     
      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS vorticidade-zc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(wz(i,j) - wze(i,j))
        end do 
      end do
     
      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS streamwisec float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(psi(i,j) - psie(i,j))
        end do 
      end do
     
      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS Txxc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(Txx(i,j) - Txxe(i,j))
        end do 
      end do

      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS Txyc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(Txy(i,j) - Txye(i,j))
        end do 
      end do

      write(5,10)'CELL_DATA',(imax-1)*(jmax-1)
      write(5,10)'POINT_DATA', imax*jmax
      write(5,*)''
      write(5,9)'SCALARS Tyyc float 1'
      write(5,9)'LOOKUP_TABLE default'
      do j = 1, jmax
        do i = 1, imax
          write(5,6) dabs(Tyy(i,j) - Tyye(i,j))
        end do 
      end do

      close(unit=5)
      
    5  format(F10.8, F11.8, F11.8)
    6  format(F13.8)
    7  format(A, I4, I4, I2)
    8  format(A, I6, A7)
    9  format(A)
   10  format(A, I6)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine write_script_plots(t)
c
c     implicit none
c     include 'par.nn'
c     character*80 nome0, nome1 , nome2, nome3, nome4, nome5 
c     integer i, j
c     real*8 x, y, t

c      write(nome0,'(a,a,a,a,a,a)')
c    & 'script','_','plots','_','new','.gnu'

c     write(nome1,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'sol_exa_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome2,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'sol_num_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome3,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'err_Max_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome4,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'err_Int_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome5,'(i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,".pdf'"

c     open(1, file = nome0)
c     write(1,*) 'reset'
c     write(1,*) '###'
c     write(1,*) 'FIL_ES = system("ls',
c    & ' -1 ',nome1,'")'
c     write(1,*) 'FIL_NS = system("ls',
c    & ' -1 ',nome2,'")'
c     write(1,*) 'FIL_EI = system("ls',
c    & ' -1 ',nome3,'")'
c     write(1,*) 'FIL_EM = system("ls',
c    & ' -1 ',nome4,'")'
c     write(1,*) 'LABEL  = system(" ")'
c     write(1,*) '### 3D plot ###'
c     write(1,*) '### Ux  ###'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_ux_",nome5
c     write(1,*) '#set multiplot layout 1,2'
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:4 with lines',
c    &' title "Ux Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:4 with lines',
c    &' title "Ux Exact" lc rgb "blue" lw 2'
c     write(1,*) '#unset multiplot'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Uy  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_uy_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:5 with lines',
c    &' title "Uy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:5 with lines',
c    &' title "Uy Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Wz  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_Wz_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:6 with lines',
c    &' title "Wz Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:6 with lines',
c    &' title "Wz Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Psi  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_Psi_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:7 with lines',
c    &' title "Psi Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:7 with lines',
c    &' title "Psi Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Txx  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_Txx_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:8 with lines',
c    &' title "Txx Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:8 with lines',
c    &' title "Txx Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Txy  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_Txy_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x" '
c     write(1,*) 'set ylabel "y" '
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:9 with lines',
c    &' title "Txy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:9 with lines',
c    &' title "Txy Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### Tyy  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_Tyy_",nome5
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    
c     write(1,*) 'splot FIL_NS using 2:3:10 with lines',
c    &' title "Tyy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:10 with lines',
c    &' title "Tyy Exact" lc rgb "blue" lw 2'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*) '### All  ###'
c     write(1,*) 'reset'
c     write(1,*) 'set terminal pdfcairo enhanced color'
c     write(1,*) 'set output ' ,"'3D_All_",nome5
c     write(1,*) '#set multiplot layout 3,5'
c     write(1,*) '#unset multiplot'
c     write(1,*) '#set logscale xy 2.718281' 
c     write(1,*) 'set dgrid3d ',imax,',',jmax
c     write(1,*) 'set hidden3d'
c     write(1,*) 'set pm3d'
c     write(1,*) 'set   autoscale'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*) '#set lmargin 18'
c     write(1,*) '#set bmargin 14'
c     write(1,*) '#set key 50.0,5.0'
c     write(1,*) '#set key font'    , " 'arial,18' "
c     write(1,*) '#set xlabel font' , " 'arial,18' "
c     write(1,*) '#set ylabel font' , " 'arial,18' "
c     write(1,*) '#set xtic font'   , " 'arial,16' "
c     write(1,*) '#set mxtics 2'
c     write(1,*) '#set ytic font'   , " 'arial,16' "
c     write(1,*) '#set mytics 2'
c     write(1,*) 'set xlabel "x"'
c     write(1,*) 'set ylabel "y"'
c     write(1,*) '#set key below right above' 
c     write(1,*) '#set key box inside right top'
c     write(1,*) 'set key inside right top'    

c     write(1,*) 'splot FIL_NS using 2:3:4 with lines',
c    &' title "Ux Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:4 with lines',
c    &' title "Ux Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:5 with lines',
c    &' title "Uy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:5 with lines',
c    &' title "Uy Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:6 with lines',
c    &' title "Wz Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:6 with lines',
c    &' title "Wz Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:7 with lines',
c    &' title "Psi Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:7 with lines',
c    &' title "Psi Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:8 with lines',
c    &' title "Txx Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:8 with lines',
c    &' title "Txx Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:9 with lines',
c    &' title "Txy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:9 with lines',
c    &' title "Txy Exact" lc rgb "blue" lw 2'

c     write(1,*) 'splot FIL_NS using 2:3:10 with lines',
c    &' title "Tyy Numerical" lc rgb "red" lw 2'
c     write(1,*) 'splot FIL_ES using 2:3:10 with lines',
c    &' title "Tyy Exact" lc rgb "blue" lw 2'
c     write(1,*) '#unset multiplot'
c     write(1,*) 'unset output'
c     write(1,*) '#pause -1 "enter para ver outro plot',
c    &' ou Ctrl+C para fechar o gnuplot" '

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_ux_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Ux - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:4',
c    & 'with image title "Ux Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Ux - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:4',
c    &'with image title "Ux Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_uy_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Uy - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:5',
c    & 'with image title "Uy Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Uy - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:5',
c    &'with image title "Uy Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_Wz_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Wz - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:6',
c    & 'with image title "Wz Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Wz - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:6',
c    &'with image title "Wz Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_Psi_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Psi - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:7',
c    & 'with image title "Psi Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Psi - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:7',
c    &'with image title "Psi Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_Txx_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Txx - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:8',
c    & 'with image title "Txx Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Txx - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:8',
c    &'with image title "Txx Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_Txy_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]'
c     write(1,*)'set cblabel "Txy - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:9',
c    & 'with image title "Txy Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'#set cbrange[-1:1]' 
c     write(1,*)'set cblabel "Ux - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:9',
c    &'with image title "Txy Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_maps_Tyy_",nome5
c     write(1,*)'set multiplot layout 1,2'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    &"2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'set cbrange[*:*]'
c     write(1,*)'set cblabel "Tyy - Numerical" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_NS using 2:3:10',
c    & ' with image title "Tyy Numerical"'
c     write(1,*)'set colorbox vertical'
c     write(1,*)'set xlabel "x"'
c     write(1,*)'set ylabel "y" norotate offset -1,0'
c     write(1,*)'# set palette defined (0.0 "#000fff",',
c    & "2 '#90ff70', 10 '#ee0000')"
c     write(1,*)'# set palette gray negative'
c     write(1,*)'# set palette gray'
c     write(1,*)'set palette defined ( 0 "green", 1 "blue",',
c    &'2 "orange", 3 "red")' 
c     write(1,*)'set cbrange[*:*]'
c     write(1,*)'set cblabel "Tyy - Exact" rotate by -90'
c     write(1,*) 'set xrange[',x0,':',x0+(imax-1)*dx,']'
c     write(1,*) 'set yrange[',y0,':',y0+(jmax-1)*dy,']'
c     write(1,*)'plot FIL_ES using 2:3:10',
c    &' with image title "Tyy Exact"'
c     write(1,*)'unset multiplot'
c     write(1,*)'unset output'

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_err_",nome5
c     write(1,*)'set logscale y 2.718281' 
c     write(1,*)'#set   autoscale'
c     write(1,*)'set xrange [*:*]'
c     write(1,*)'set yrange [*:*]'
c     write(1,*)'##set lmargin 18'
c     write(1,*)'##set bmargin 14'
c     write(1,*)'##set key 50.0,5.0'
c     write(1,*)'#### Fonte'
c     write(1,*)'#set key font'   ,"'arial,18'"
c     write(1,*)'#set xlabel font',"'arial,18'"
c     write(1,*)'#set ylabel font',"'arial,18'"
c     write(1,*)'#### Eixos coordenados'
c     write(1,*)'#set xtic font'," 'arial,16' "
c     write(1,*)'#set mxtics 2'     
c     write(1,*)'#set ytic font'," 'arial,16'"  
c     write(1,*)'#set mytics 2'     
c     write(1,*)'#### Legenda dos eixos'
c     write(1,*)'set xlabel "t" '
c     write(1,*)'set ylabel "Erro" '
c     write(1,*)'#### Plotagem dos graficos'
c     write(1,*)'set key below right above' 
c     write(1,*)'##set key box inside right top'
c     write(1,*)'##set title font', "'arial,22'"
c     write(1,*)'#set key inside right top' 

c     write(1,*)'plot FIL_EI using 1:2',
c    & ' title "Ux - Integral error" ',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:2',
c    & ' title "Ux - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:3',
c    & ' title "Uy - Integral error"',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:3',
c    & ' title "Uy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:4',
c    & ' title "Wz - Integral error"',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:4',
c    & ' title "Wz - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:5',
c    & ' title "Psi - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:5',
c    & ' title "Psi - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:6',
c    & ' title "Txx - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:6',
c    & ' title "Txx - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:7',
c    & ' title "Txy - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:7',
c    & ' title "Txy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:8', 
c    & ' title "Tyy - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:8',
c    & ' title "Tyy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'unset output'
c     write(1,*)'#pause -1 "enter para ver outro plot ou Ctrl+C,
c    & para fechar o gnuplot"'
c     write(1,*)'#reread'
c     close (unit=1)
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine write_script_plots_error(t)
c
c     implicit none
c     include 'par.nn'
c     character*80 nome0, nome1 , nome2, nome3, nome4, nome5 
c     integer i, j
c     real*8 x, y, t

c      write(nome0,'(a,a,a,a,a,a)')
c    & 'script','_','plots','_','new_error','.gnu'

c     write(nome1,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'sol_exa_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome2,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'sol_num_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome3,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'err_Max_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome4,'(a,i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & 'err_Int_',imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,'.dat'

c     write(nome5,'(i0.4,a,f0.2,a,f0.2,a,f0.2,a,f0.8,
c    & a,f0.8,a,i0.2,a,i0.2,a)')
c    & imax,'_',re,'_',beta,'_',
c    & wi,'_',alpha_G,'_',dt,'_',type_sim,
c    & '_',type_mms,".pdf'"

c     open(1, file = nome0)
c     write(1,*) 'reset'
c     write(1,*) '###'
c     write(1,*) 'FIL_ES = system("ls',
c    & ' -1 ',nome1,'")'
c     write(1,*) 'FIL_NS = system("ls',
c    & ' -1 ',nome2,'")'
c     write(1,*) 'FIL_EI = system("ls',
c    & ' -1 ',nome3,'")'
c     write(1,*) 'FIL_EM = system("ls',
c    & ' -1 ',nome4,'")'

c     write(1,*)'reset'
c     write(1,*)'set terminal pdfcairo enhanced color'
c     write(1,*)'set output',"'2D_err_",nome5
c     write(1,*)'set logscale y 2.718281' 
c     write(1,*)'#set   autoscale'
c     write(1,*)'set xrange [*:*]'
c     write(1,*)'set yrange [*:*]'
c     write(1,*)'##set lmargin 18'
c     write(1,*)'##set bmargin 14'
c     write(1,*)'##set key 50.0,5.0'
c     write(1,*)'#### Fonte'
c     write(1,*)'#set key font'   ,"'arial,18'"
c     write(1,*)'#set xlabel font',"'arial,18'"
c     write(1,*)'#set ylabel font',"'arial,18'"
c     write(1,*)'#### Eixos coordenados'
c     write(1,*)'#set xtic font'," 'arial,16' "
c     write(1,*)'#set mxtics 2'     
c     write(1,*)'#set ytic font'," 'arial,16'"  
c     write(1,*)'#set mytics 2'     
c     write(1,*)'#### Legenda dos eixos'
c     write(1,*)'set xlabel "t" '
c     write(1,*)'set ylabel "Erro" '
c     write(1,*)'#### Plotagem dos graficos'
c     write(1,*)'set key below right above' 
c     write(1,*)'##set key box inside right top'
c     write(1,*)'##set title font', "'arial,22'"
c     write(1,*)'#set key inside right top' 

c     write(1,*)'plot FIL_EI using 1:2',
c    & ' title "Ux - Integral error" ',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:2',
c    & ' title "Ux - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:3',
c    & ' title "Uy - Integral error"',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:3',
c    & ' title "Uy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:4',
c    & ' title "Wz - Integral error"',
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:4',
c    & ' title "Wz - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:5',
c    & ' title "Psi - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,',
c    & ' FIL_EM using 1:5',
c    & ' title "Psi - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:6',
c    & ' title "Txx - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:6',
c    & ' title "Txx - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:7',
c    & ' title "Txy - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:7',
c    & ' title "Txy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'plot FIL_EI using 1:8', 
c    & ' title "Tyy - Integral error"', 
c    & ' lc rgb "red" lw 1.5 lt 3 dt 5,', 
c    & ' FIL_EM using 1:8',
c    & ' title "Tyy - Maximal error"',
c    & ' lc rgb "blue" lw 1.5 lt 2 dt 5'

c     write(1,*)'unset output'
c     write(1,*)'#pause -1 "enter para ver outro plot ou Ctrl+C,
c    & para fechar o gnuplot"'
c     write(1,*)'#reread'
c     close (unit=1)
c     end
