c     this code compares energies in 3D vs original 1D profile
      program en_comparison

      implicit none

      real*8 Msun, Rsun, G
      parameter(Msun=1.9891d33,Rsun=6.9599d10,G=6.67390d-8) 
      
      integer nmax_sph
      parameter (nmax_sph=500000)

      integer nmax_vars, nvars
      parameter (nmax_vars=50)

      integer :: num_args, ix
      character(len=32), dimension(:), allocatable :: args
      
      integer i,j, k, l, modnum, maxn, maxl, it, im, ir
      integer maxnr, maxt
      integer num1, num2, numf
      integer max_var, max_mesa, max_m

      real*8  eint_sph(nmax_sph), epot_sph(nmax_sph), m_sph(nmax_sph)
      real*8  eint_mesa(10000),epot_mesa(10000),m_mesa(10000)
      real*8  ei_sp,ep_sp, dm, dummy7
      real*8  r_mesa(10000),dummy1,dummy2, dummy3,dummy4,dummy5,dummy6
      real*8  eint_accum, epot_accum, m_error
      real*8  eint_error(10000), epot_error(10000), dm_err
      
      real mcore, mcomp
      real ce(nmax_vars, nmax_sph) 
      real maccum, mbot,mtop, rbot, rtop
      integer cestat(nmax_sph), ces
      real*8 mass_min, mass_max, age_min, age_max
      real*8 rad_max, rad_min, dr
      real *8 z, fmin, fmax
      real wconv, w
      real xco, yco, zco, xcm, ycm, zcm
c      real*8 age, mcur, var(25) 
       real*8 rcur, vcur, drcur
       integer nir, ir1, ir2,irc

      character*60 file_sph, file_mesa,file_out
      character*3  cebin
      character*4  lif, prnt
      character*5 prefix

      real*8 epot, eint, menv, orb_sep
      real*8 e_orb_start, e_orb
      real q, rl, q3

      logical it_exists

c     attantion! check with sortit.f how many variables are not in the output
      nvars=26
      dm_err=0.01               ! dm fow which accumulated energy error is found
      
c default values      
      file_sph='sorted_0000.dat'
      file_mesa='mesa.starsmash'

c     parsing input
      num_args = command_argument_count()
      
      write(*,*) "number of arguments", num_args
      if(num_args.ge.0) then
         
         allocate(args(num_args))  
         
         do ix = 1, num_args
            call get_command_argument(ix,args(ix))
         end do               
         
         do ix = 1, num_args
            if(args(ix) == '-h') then
               write(*,*) "Energy comparison between SPH model and initial MESA model."
               write(*,*) "SPH model has to be sorted first by ./sorted "
               write(*,*) "MESA model used for SPH has to be parsed first by port_mesa_to_starsmash.py"
               write(*,*) "usage (in [] are the default values):"
               write(*,*) "-nsph  XXXX [sorted_0000.dat]input file after sorting "
               write(*,*) "-nmesa XXXX [mesa.starsmash] parsed 1D mesa profile"
               goto 42
            end if
            
            if(args(ix) == '-nsph') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for nsph"
                  goto 42
               end if
               read(args(ix+1),*) file_sph
               write(*,*) "use SPH file ", file_sph
            end if
            if(args(ix) == '-nmesa') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for nmesa"
                  goto 42
               end if
               read(args(ix+1),*) file_mesa
               write(*,*) "use SPH file ", file_mesa
            end if
             
         end do
      end if
      
      inquire(file=file_sph,exist=it_exists)
      if (it_exists) then 
         OPEN(1,file=file_sph,STATUS='old')
         write(*,*) 'opened ', file_sph
      else
         write(*,*) 'cannot open ', file_sph
         write(*,*) 'Run as "./enc -h"'
         goto 42
      endif
      
      inquire(file=file_mesa,exist=it_exists)
      if (it_exists) then 
         OPEN(2,file=file_mesa,STATUS='old')
         write(*,*) 'opened ', file_mesa
      else
         write(*,*) 'cannot open ', file_mesa
         write(*,*) 'Run as "./enc -h"'
         goto 42
      endif
      
      
c     input has that many lines (maximum)
      maxl=nmax_sph
      
      do i=1,maxl
         do j=1,nvars
            ce(j,i)=0
         end do
         cestat(i)=0
      end do


c     read this profile file
c     xco, --- compact obkect coordinate
c     xcm center of mass of the binay coordinate
      read(1,*) maxl, mcore, mcomp, xco, yco, zco, 
     &     xcm, ycm, zcm
      
      orb_sep = sqrt((xco)**2. + (yco)**2.+
     &     (zco)**2)
      write(*,*) "orbital separation ", orb_sep, mcomp
         
c     
      do i=1,maxl
         read(1, *) (ce(j,i), j=1,27), cestat(i)
!     current order:
!     mass, rad, dm_particle, P, rho, entropy, eu, ek, ep, ep_sphere, a.m., a.m.c.m.,
!     rad to c.m. of the binary, rad to z via c.o.of the binary, rad to rgb z,
!     w to binary, w to rg core status, dev, kappa, t, x,y,z
      end do
      close(1)

      eint_accum=0.
      epot_accum=0.
      do i=maxl,1,-1
         j=maxl-i+1
         m_sph(j)=ce(1,i)
         eint_accum=eint_accum+(ce(7,i)/1.e13) *ce(3,i)*(Msun/1.e33)
         eint_sph(j)=eint_accum
         epot_accum=epot_accum+(ce(10,i)*100.)*ce(3,i)*(Msun/1.e33)
         epot_sph(j)=epot_accum         
      end do

      
c     read mesa file
      eint_accum=0.
      epot_accum=0.
      max_mesa=20000
      do i=1, max_mesa
         read(2,*,end=99) m_mesa(i), r_mesa(i), dummy1, dummy2, dummy3,
     &        dummy4, dummy5, dummy6, ei_sp, dummy7
         if(i.gt.1) then
            dm=m_mesa(i-1)-m_mesa(i)
            ep_sp=-(G*Msun/1.e33)*(Msun/1.e13/Rsun)*m_mesa(i)/(10.**r_mesa(i))
            epot_accum=epot_accum+dm*ep_sp
            eint_accum=eint_accum+(ei_sp/1.e13)*dm*(Msun/1.e33)
            eint_mesa(i)=eint_accum
            epot_mesa(i)=epot_accum
         end if
      end do
      close(2)
99    max_mesa=i-1
      write(*,*) 'total lines in mesa star profile ', max_mesa

      
      
         
      open(2,file="energy_mesa.dat",status='unknown')
      do i=1,max_mesa       
         write(2,*) m_mesa(i), eint_mesa(i), epot_mesa(i)
      end do
      close(2)
      
      open(2,file="energy_sph.dat",status='unknown')
      do i=1,maxl
         write(2,*) m_sph(i), eint_sph(i), epot_sph(i)
      end do
      close(2)


      max_m=(m_sph(1)-m_sph(maxl))/dm_err-1
      
      open(2,file="energy_comparison.dat",status='unknown')
      j=1
      k=1
      do i=1,max_m
         m_error=m_sph(1)-i*dm_err
         do while (m_mesa(k).gt.m_error)
            k=k+1
         end do
         do while (m_sph(j).gt.m_mesa(k))
            j=j+1
         end do
         write(2,100) m_error,  m_sph(j), m_mesa(k),
     &        eint_sph(j), eint_mesa(k),
     &        (eint_sph(j)-eint_mesa(k))/eint_mesa(k),
     &        epot_sph(j), epot_mesa(k),
     &        (epot_sph(j)-epot_mesa(k))/epot_mesa(k)
      enddo
      close(2)
 100     format(9e13.5)
         
 42   continue
      
      stop
      end

      
