      subroutine orbital_parameter
      include 'common_sph_var.h'
      integer i,nstar,idmax1,idmax2
      real*8 P,a,e2,rhomax1,rhomax2,m1,m2
      common/sep_max/a,m1,m2
      real*8 mtot,lz,lx,ly,l
      real*8 Eorb
      real*8 cmmbin,cmxbin,cmybin,cmzbin,cmvxbin,cmvybin,cmvzbin
      real*8 x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2
      common/l_mom/lx,ly,lz
      common/cofm_bin/cmmbin,cmxbin,cmybin,cmzbin,cmvxbin,cmvybin,cmvzbin
      common/bin_par/x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2
      mtot=m1+m2
C Period assuming Keplerian orbit
      P = sqrt(4.d0*pi**2*a**3/mtot)*tunit
      lx = m1*((y1-cmybin)*(vz1-cmvzbin)-(z1-cmzbin)*(vy1-cmvybin))
     &   + m2*((y2-cmybin)*(vz2-cmvzbin)-(z2-cmzbin)*(vy2-cmvybin))
      ly = m1*((z1-cmzbin)*(vx1-cmvxbin)-(x1-cmxbin)*(vz1-cmvzbin))
     &   + m2*((z2-cmzbin)*(vx2-cmvxbin)-(x2-cmxbin)*(vz2-cmvzbin))
      lz = m1*((x1-cmxbin)*(vy1-cmvybin)-(y1-cmybin)*(vx1-cmvxbin))
     &   + m2*((x2-cmxbin)*(vy2-cmvybin)-(y2-cmybin)*(vx2-cmvxbin))
      l = sqrt(lx**2+ly**2+lz**2)
      Eorb = 0.5d0*m1*(vx1**2+vy1**2+vz1**2)
     &      +0.5d0*m2*(vx2**2+vy2**2+vz2**2)
     &      -m1*m2/a      
      e2 = 1.d0 - l**2*mtot/(m1**2*m2**2*a)
      write(*,*)'Parameter ',m1,m2 
      write(*,*)'angular momentum (CM)',lx*lunit,ly*lunit,lz*lunit
      write(*,*)t*tunit,a,P,e2,Eorb,lx,ly,lz," checkecc1"
      e2 = 1.d0+2.0d0*Eorb*l**2*mtot/(m1**3*m2**3)
      write(*,*)t*tunit,a,P,e2,Eorb,lx,ly,lz," checkecc2"
      if (e2.lt.0.d0) then 
         write(*,*)'Eccectricity error',e2
         e2 = 0.d0
      endif 
      write(*,*)t*tunit,sqrt(e2),' ecc' 
      write(54,*)t*tunit,a,P,sqrt(e2)
      end
