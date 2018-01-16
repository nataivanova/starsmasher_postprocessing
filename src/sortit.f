      program sortit

      implicit none

      integer nmax_sph
      parameter (nmax_sph=800000)
      integer max_im, max_in
      parameter (max_im=5000, max_in=1000)

      integer :: num_args, ix
      character(len=12), dimension(:), allocatable :: args
      integer :: ni,nf,ns
      
      integer i,j, k, l, modnum, maxn, maxl, it, im, n, m
      integer num1, num2, numf
      integer sort, sortcount 
      integer*4 sortcounttot

      real*8 dm, dt, dmloc
      real ce(25, nmax_sph) 
      real ces(25, nmax_sph) 
      real ceone(25), maccum, mbot,mtop
      integer cestat(nmax_sph)
      integer cestats(nmax_sph)
      integer list(nmax_sph)
      integer sr(max_im,max_in), listsr(max_im),innmax(max_im)
      real*8 mass_min, mass_max, age_min, age_max
      real wconv, w, rmin, rcur

      real epot1, rv
      integer icheck

      real xco, yco, zco, rzcb, rzcc, rcb
      real xcm, ycm, zcm, wcb, wcc

      logical it_exists,it_existsb
      
c      real*8 age, mcur, var(25) 

      character*60 filein, fileout
      character*3  cebin
      character*4  lif, prnt

      integer inmax, immax, inc, imc, counter, nmin
      real  tstart, tfinish
      
      call cpu_time(tstart)

c parsing input
c     defaults
      ni=0
      nf=1000000
      ns=1
      
c     parsing input
      num_args = command_argument_count()

c      write(*,*) "number of arguments", num_args
      if(num_args.ge.0) then

         allocate(args(num_args))  
      
         do ix = 1, num_args
            call get_command_argument(ix,args(ix))
         end do      
         

         do ix = 1, num_args
            if(args(ix) == '-h') then
               write(*,*) "Sortes out an sph object in an 1D profile"
               write(*,*) "usage (in [] are the default values):"
               write(*,*) "-ni XXXX [0]       start processing from outXXXX.sph "
               write(*,*) "-nf XXXX [1000000] end processing at outXXXX.sph "
               write(*,*) "-ns X [1] skip each X file"
               write(*,*) "          it skips all files it can't find, but its faster if ns is given"
               goto 42
            end if
            if(args(ix) == '-ni') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for ni"
                  goto 42
               end if
               read(args(ix+1),*) ni
               write(*,*) "start from file ", ni
               if(ni.lt.0) then
                  write(*,*) "ni: only can work with values >0"
                  goto 42
               end if
            end if
            if(args(ix) == '-nf') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for nf"
                  goto 42
               end if
               read(args(ix+1),*) nf
               write(*,*) "end at file ", nf
               if(nf.lt.0.or.nf.lt.ni) then
                  write(*,*) "nf: only can work with nf > ni"
                  goto 42
               end if
            end if
            if(args(ix) == '-ns') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for ns"
                  goto 42
               end if
               read(args(ix+1),*) ns
               write(*,*) "skip ", ns, "files"
               if(ns.lt.0) then
                  write(*,*) "ns: only can work with ns > 0"
                  goto 42
               end if
            end if
         end do
      end if
      
c     input has that many lines (maximum)
      maxl=nmax_sph


      do numf=ni,nf,ns

         write(LIF,'(I4)') numf
         read(LIF,'(4A)') prnt
         do i=1,4 
            if(prnt(i:i).EQ.' ') prnt(i:i)='0'
         end do
         filein='prof_'//prnt//'.dat'
         inquire(file=filein,exist=it_exists)
         if (it_exists) then 
            OPEN(1,file=filein,STATUS='old', err=1002)
            write(*,*) 'opened ', filein
         else
            write(*,*)'Skipping the file',filein
            cycle
         endif


c     read this profile file

         do i=1,nmax_sph
            do j=1,25
               ce(j,i)=0
            end do
            cestat(i)=0
            list(i)=i
         end do
         maxl=nmax_sph

         xco=0.
         yco=0.
         zco=0.
         read(1, *, end=99)
C     sanity check: this has to have same number of inputs as output from classifaction.f
         do i=1,maxl
            read(1, *, end=99) (ce(j,i), j=1,19), cebin
            cestat(i)=1
            if(cebin(1:1).EQ.'c') cestat(i)=2
            ce(25,i)=sqrt( (ce(1,i)-ce(1,1))**2 
     &           + (ce(2,i)-ce(2,1))**2
     &           + (ce(3,i)-ce(3,1))**2)
            if(ce(4,i).ge.0.049.and.i.gt.100.and.ce(5,i).le.1e-20) then
               maxl=i
               xco=ce(1,i)
               yco=ce(2,i)
               zco=ce(3,i)
               exit
            end if
         end do
         close(1)

 99      maxl=i-1
         write(*,*) 'total lines', maxl
         write(*,*) "core particle", ce(4,1)

         inmax=100      
         immax=(1000*ceiling(maxl/1000.))/inmax
         if(immax.gt.5000) then
            inmax=300
            immax=maxl/inmax+1
         end if
         write(*,*) "sorting using ", inmax, immax
         write(*,*) "reading from ", ni, " to ", nf , " skipping ", ns 

         do m=1, immax
            do n=1,inmax
               sr(m,n)=-1
            end do
            listsr(m)=1
         end do

         counter=1
         do m=1, immax
            do n=1,inmax
               counter =counter +1
               sr(m,n) = counter
               innmax(m)=n  
               if(counter.eq.(maxl-1)) exit
            end do
            if(counter.eq.(maxl-1)) exit
         end do
         inc=n
         imc=m
         if(counter.lt.maxl-1) then
            write(*,*) 'not succient array, recompile ', 
     &           inc, imc, counter, maxl
            stop
         end if

         write(*,*) "mapped to ", inc, imc, inmax,immax, counter

         do m=1, imc
            sortcount=0
            sortcounttot=0
            sort=0
            do 
               sort=1
               do n=1, innmax(m)-1
                  if(m.eq.imc.and.n.gt.inc) exit
                  if(ce(25,sr(m,n)).gt.ce(25,sr(m,n+1))) then
                     sort=sr(m,n+1)
                     sr(m,n+1)=sr(m,n)
                     sr(m,n)=sort
                     sort=0
                     sortcounttot=sortcounttot+1
                  end if
               end do
               if(sort.eq.1) then
                  exit
               end if
               sortcount=sortcount+1
            end do
c            write(*,*) "made ", sortcount, " ", sortcounttot, 
c     &           "steps for m=", m, inmax

c     debug checkup
c     goto 7001
            
            do i=1, inmax-1
               if(m.eq.imc.and.j.gt.inc) exit
               do j=i+1,inmax-2
                  if(m.eq.imc.and.j.gt.inc-1) exit                  
                  if(sr(m,i).eq.sr(m,j)) then
                     write(*,*) "sorting error for m=", m, list(i), list(j)
                     stop
                  end if
               end do
            end do
 7001       continue

         end do

         counter=1
         do 
            nmin=-1
            if(imc.gt.1) then
               rmin=1.e10
               do m=1, imc
                  if(sr(m,listsr(m)).ge.1) then
                     rcur = ce(25,sr(m,listsr(m)))
                     if(rcur.lt.rmin) then
                        rmin=rcur
                        nmin=m
                     end if
                  end if
               end do
            else
               nmin=1               
            end if
            if(nmin.eq.-1) then
c     this error comes out for some reasons if arrays are not sufficiantly mapped. To fix.
               write(*,*) "error in finding minimum element", 
     &              counter, maxl, m, rmin, rcur
               stop
            end if
            
            counter =counter +1
c            write(*,*) counter, nmin, listsr(nmin), 
c     &           sr(nmin,listsr(nmin)),
c     &           ce(13,sr(nmin,listsr(nmin))) 
            list(counter)=sr(nmin,listsr(nmin))
            listsr(nmin)=listsr(nmin)+1
            if(counter.eq.(maxl-1)) exit
         end do                     


c     this is rather self-consistency check, only do for small files
         if(maxl.le.100000) then 
            do i=2, maxl-2
               do j=i+1,maxl-2
                  if(list(j).eq.list(i)) then
                     write(*,*) "sorting error", list(i), list(j)
                     stop
                  end if
               end do
            end do
         end if


         do i=2, maxl
            if(ces(25,i-1).gt.ces(25,i)) then
               write(*,*) "sorting error", i-1, ces(25,i-1),ces(25,i)
            end if
         end do

         do i=2, maxl-1
            cestats(i)=cestat(list(i))
            do j=1,25
               ces(j,i)=ce(j,list(i))
            end do 
         end do

         write(*,*) "sorted ", filein, sortcount, sortcounttot, ce(25,2)
         write(*,*) "core particle", ce(4,1), ces(4,1)

         fileout='sorted_'//prnt//'.dat'
         OPEN(2,FILE=fileout,STATUS='UNKNOWN', err=1003)
         maccum=ce(4,1)
         xcm=ce(4,maxl)*xco/(ce(4,1) + ce(4,maxl))
         ycm=ce(4,maxl)*yco/(ce(4,1) + ce(4,maxl))
         zcm=ce(4,maxl)*zco/(ce(4,1) + ce(4,maxl))
         write(2,*) maxl-2, ce(4,1), ce(4,maxl), xco, yco, zco, 
     &              xcm, ycm, zcm
         wconv=1.989*6.96*6.96*10.
         ces(25,maxl)= sqrt( (ce(1,maxl))**2 
     &        + (ce(2,maxl))**2
     &        + (ce(3,maxl))**2)
         do i=2,maxl-1
            maccum=maccum+ces(4,i)
!************
            epot1=0 ! this is "inside a sphere" potential energy  
            do icheck=2,maxl
               if(i.ne.icheck) then
                  rv=sqrt( (ces(1,i)-ces(1,icheck))**2 
     &                 + (ces(2,i)-ces(2,icheck))**2
     &                 + (ces(3,i)-ces(3,icheck))**2)
                  if(ces(25,icheck).le.ces(25,i)) then
                     epot1=epot1-1.907*ces(4,icheck)/rv
                  else
                     exit
                  end if
               end if
            end do
            rv=sqrt( (ce(1,1)-ces(1,i))**2 
     &           + (ce(2,1)-ces(2,i))**2
     &           + (ce(3,1)-ces(3,i))**2)
            epot1=epot1-1.907*ce(4,1)/rv  

c     this is if there is a second core
            if(ces(25,maxl).le.ces(25,i)) then
               rv=sqrt( (ces(1,i)-ce(1,maxl))**2 
     &              + (ces(2,i)-ce(2,maxl))**2
     &              + (ces(3,i)-ce(3,maxl))**2)
               epot1=epot1-1.907*ce(4,maxl)/rv  
            end if 
          
c     distance to the z- axes (later should be replace by rotation axes) going through the center of mass of the binary
            rcb=sqrt( (ces(1,i)-xcm)**2. 
     &           + (ces(2,i)-ycm)**2. + (ces(3,i)-zcm)**2.)
            rzcb=sqrt( (ces(1,i)-xcm)**2. 
     &           + (ces(2,i)-ycm)**2.)
            rzcc=sqrt( (ces(1,i))**2. 
     &           + (ces(2,i))**2.)
            wcb=0.
            wcc=0.
            if(rzcb.gt.0) 
     &           wcb=ces(11,i)/wconv/rzcb/rzcb/ces(4,i)
            if(rzcc.gt.0) 
     &           wcc=ces(12,i)/wconv/rzcc/rzcc/ces(4,i)

            write(2,100) maccum, ces(25,i), ces(4,i), ! mass, rad, dm_particle
     &           ces(5,i),ces(6,i),ces(9,i), ! P, rho, entropy
     &           ces(7,i),       ! specific internal energy
     &           ces(8,i),       ! kinetic energy
     &           ces(10,i),      ! total potential energy
     &           epot1,                ! potential energy to the "inner" sphere
     &           ces(11,i), ces(12,i), ! z-a.m. to c.m., z-a.m. to r.g.,
     &           rcb, rzcb, rzcc,      ! distance to  c.m. of bin, z-rad_to cm of bin, z-rad to rg core,
     &           wcb, wcc, 
     &           ces(13,i), ces(14,i), ces(15,i), ! H+, He+, He++
     &           ces(16,i), ces(17,i), ces(18,i), ! div, kappa, T
     &           ces(1,i),ces(2,i),ces(3,i), ces(19,i),! x,y,z, h
     &           cestats(i)  ! w to c.o. of bin, w to the core, status     
         end do
         close(2)
 100     format(2e16.8, 25e14.6, i3)
         
      end do

      goto 10000

 1002 write(*,*) " could not find file", filein

      goto 10000

 1003 write(*,*) " could not open file", fileout

 42   continue

10000 call cpu_time(tfinish)

      write(*,*) "Spent ",tfinish-tstart , "seconds"
      

      STOP 
      END

      
