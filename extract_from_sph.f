      program main
      
      include 'common_sph_var.h'
      
      integer :: i,ni,nf,ns, flag_read
      integer :: case_eos, case_run
      integer :: num_args, ix
      character(len=12), dimension(:), allocatable :: args
      
      common/case_id_all/ case_eos, case_run


c     defaults
      case_run=0    ! 1 is relaxation run, etc, this is same as nrelax
      case_eos=0
      ni=0
      nf=1000000
      ns=1
           
c     parsing input
      num_args = command_argument_count()

      write(*,*) "number of arguments", num_args
      if(num_args.ge.0) then
         allocate(args(num_args))  
      
         do ix = 1, num_args
            call get_command_argument(ix,args(ix))
         end do      
         
            
         do ix = 1, num_args
            
c            write (*,*) "processing  argument ", ix, args(ix), num_args
            if(args(ix) == '-h') then
               write(*,*) "Extract a set of variables from binary outXXXX.sph files into ascii prof_XXXX.dat files"
               write(*,*) "usage (in [] are the default values):"
               write(*,*) "-eos [0],1] for tabulated EOS or not"
               write(*,*) "-run [0],1,2,3 same as nrelax"
               write(*,*) "-ni XXXX [0]       start processing from outXXXX.sph "
               write(*,*) "-nf XXXX [1000000] end processing at outXXXX.sph "
               write(*,*) "-ns X [1] skip each X file"
               goto 42
            end if
            if(args(ix) == '-eos') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for EOS", ix+1, num_args
                  goto 42
               end if
               read(args(ix+1),*) case_eos
               write(*,*) "using for EOS case ", case_eos
               if(case_eos.ne.0.and.case_eos.ne.1) then
                  write(*,*) "only can work with eos values 0 or 1"
                  goto 42
               end if
            end if
            if(args(ix) == '-run') then
               if((ix+1).gt.num_args) then
                  write(*,*) "need an argument for run"
                  goto 42
               end if
               read(args(ix+1),*) case_run
               write(*,*) "using for run type ", case_run
               if(case_run.ne.0.and.case_run.ne.1.and.case_run.ne.2.and.case_run.ne.3) then
                  write(*,*) "only can work with run type values 0, 1, 2 , 3"
                  goto 42
               end if
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
      
      if(case_eos.eq.0) then
         write(*,*) "Processing using tabulated EOS"
      end if

      if(case_run.eq.1) then
         write(*,*) "Processing a single star"
      else
         write(*,*) "Processing a star in a binary"
      end if      

      write(*,100)'Units(t [day],E [erg],v[km/s]):',tunit,eunit,vunit,lunit
 100  format(A,ES15.6,ES15.6,ES15.6,ES15.6)

      if (case_eos.eq.0) then
         call readineostable
         call opacfile
      endif

      open(10,file='orbit.dat',status='unknown')
 
      do i=ni,nf,ns
         call read_file(i,flag_read)
         if(flag_read.eq.0) then
            write(*,*) "process ", ni
            call classification
            if(case_run.eq.0) then
               call ejecta
            end if
         end if
      enddo

      close(10)

 42   continue ! happy end

      end
