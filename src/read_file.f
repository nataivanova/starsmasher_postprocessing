      subroutine read_file(ni,flag_read)
c this subroutine reads the sph output files
      include 'common_sph_var.h'
      integer ni,i, ncheck, flag_read
      real*8 dummy
      logical it_exists,it_existsb
      character*13 filename

      flag_read=0
      write(filename,100)ni
 100  format('out',I4.4,'.sph')
 101  format('out',I5.5,'.sph')
      inquire(file=filename,exist=it_exists)
      if(it_exists) then 
          open(10,file=filename,form='unformatted')
          print*,'Reading:',filename
      else
          write(filename,101)ni
          inquire(file=filename,exist=it_existsb)
          if (it_existsb) then 
               open(10,file=filename,form='unformatted')
               print*,'Reading:',filename
          else
               write(*,*)'Skipping the file',filename
               flag_read=1
               goto 42
          endif
      endif
c read the first line of out*.sph which contains the metadata
      read(10)ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,t,
     &          nav,alpha,beta,tskipahead,ngr,nrelax,trelax,dt,omega2
c if this is the first file, dump its output to the stdout
      if (nout.eq.0) then
          write(*,*) "File's header:"
          write(*,102)ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,t,
     &          nav,alpha,beta,tskipahead,ngr,nrelax,trelax,dt,omega2
 102  format(2I6,5ES15.6,2I6ES15.6,I6,3ES15.6,2I6,3ES15.6)
      endif
c     nata's brutal fix for file's numbering, no idea why Jose's way is getting broken
      nout = ni
      do i=1,ntot
         read(10) x(i),y(i),z(i),m(i),h(i),rho(i),vx(i),vy(i),vz(i),
     &            vxdot(i),vydot(i),vzdot(i),u(i),udot(i),grpot(i),
     &            dummy,cc(i),divv(i)
         mu(i)=dummy/mp
      enddo
      read(10) ncheck

      if (ncheck.ne.ntot) then 
          write(*,*)'The file is corrupt'
          stop
      endif

      write(*,*)filename,' is read! This file is controlled by nrelax=', nrelax
      close(10)

 42   continue ! the classic asnwer collects the exits
      
      end
