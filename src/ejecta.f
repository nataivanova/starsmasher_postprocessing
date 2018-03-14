      subroutine ejecta
C this subroutine calculates the evolution quantities of a binary star
c this includes ejecta, circumbinary, and binary, as well as the separation
c note that it needs id(i) calculated by classification.f, the better the routine, the more accurate
c this subroutine
      include 'common_sph_var.h'
      integer i
C Bound material
      real*8 lxb,lyb,lzb,mb
      real*8 ekb,epb,eib
C Ejecta
      real*8 ep,munb,ekunb,epunb,eiunb
      real*8 lxunb,lyunb,lzunb
C Total angular momentum and energy
      real*8 lx,ly,lz,etot
      real*8 lxi,lyi,lzi
      real*8 epot, eint, ekin
C Velocities
      real*8 vtan,vrad,r
C Circumbinary disk
      real*8 mc,lxc,lyc,lzc,ekc,epc,eic
C Number of particles
      integer cmb,cmunb,cmc

      cmb=0
      cmunb=0
      cmc=0
      munb=0.d0
      ekunb=0.d0
      epunb=0.d0
      eiunb=0.d0
      ekb=0.d0
      epb=0.d0
      eib=0.d0
      etot = 0.d0
      lx=0.d0
      ly=0.d0
      lz=0.d0
      lxunb=0.d0
      lyunb=0.d0
      lzunb=0.d0
      lxb=0.d0
      lyb=0.d0
      lzb=0.d0
      mb=0.d0
      lxc=0.d0
      lyc=0.d0
      lzc=0.d0
      mc=0.d0
      ekc=0.d0
      epc=0.d0
      eic=0.d0

      do i=1,ntot
c energy for each particle
         epot = m(i) * grpot(i)
         eint = m(i) * u(i)
         ekin = 0.5d0 * m(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2)
c angular momentum for each particle
         lxi = m(i) * (y(i)*vz(i)-z(i)*vy(i))
         lyi = m(i) * (z(i)*vx(i)-x(i)*vz(i))
         lzi = m(i) * (x(i)*vy(i)-y(i)*vx(i))
c use the id(i) to determine where the particles are
         if (id(i).eq.4) then ! this will be ejecta
            munb = munb+m(i)
            ekunb = ekunb + ekin
            eiunb = eiunb + eint
            epunb = epunb + epot
            lxunb = lxunb + lxi
            lyunb = lyunb + lyi
            lzunb = lzunb + lzi
            cmunb = cmunb + 1
         elseif(id(i).le.2) then ! this will be binary star
            mb = mb + m(i)
            ekb = ekb + ekin
            eib = eib + eint
            epb = epb + epot
            lxb = lxb + lxi
            lyb = lyb + lyi
            lzb = lzb + lzi
            cmb = cmb + 1
         else !then circumbinary
            mc = mc + m(i)
            ekc = ekc + ekin
            eic = eic + eint
            epc = epc + epot
            lxc = lxc + lxi
            lyc = lyc + lyi
            lzc = lzc + lzi
            cmc = cmc + 1
         endif
C Get the total angular momentum and energy
         etot = etot + ekin + 0.5d0 * epot + eint
         lx = lx + lxi
         ly = ly + lyi
         lz = lz + lzi
      enddo
      write(50,*)t*tunit,etot*eunit,sqrt(lx**2+ly**2+lz**2)
      write(51,*)t*tunit,munb,ekunb*eunit,eiunb*eunit,epunb*eunit,cmunb
      write(52,*)t*tunit,mb,ekb*eunit,eib*eunit,epb*eunit,cmb
      write(55,*)t*tunit,mc,ekc*eunit,eic*eunit,epc*eunit,cmc
      write(53,*)t*tunit,lx,ly,lz,lxunb,lyunb,lzunb,lxb,lyb,lzb,lxc,lyc,lzc

      end
