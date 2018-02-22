      subroutine opacfile
      !include 'common_sph_var.h'
      implicit none
      integer i,j
      logical fileexists
      integer nummaxR,nummaxT
      parameter(nummaxR=19,nummaxT=139)
      real*8 logR(nummaxR),logT(nummaxT)
      real*8 logopac(nummaxR,nummaxT)
      common/opac_low/logT,logR,logopac
      character*5 init
      character*40 datafile
 
      datafile = 'new_table_opal_kap.tron'
      inquire(file=trim(datafile),exist=fileexists)
      if (.not.fileexists) then
         write(*,*)' *** error: '//trim(datafile)//': file not found ***'
         call exit(1)
      endif

      open(12,file=datafile,status='old')

      read(12,*)init,(logr(j),j=1,nummaxR)
      do i=1,nummaxT
        read(12,*)logT(i),(logopac(j,i),j=1,nummaxR)
        !write(*,*)logT(i)
      enddo
      close(12)
      end
