      subroutine classification
      include 'common_sph_var.h'
      
      integer i,iter
      real*8 etest,sep,phi1,xL1
      real*8 xL1p,yL1p,theta,zL1p
      real*8 rcmi,r1,r2
      integer idmax1,idmax2 
      real*8 rhomax1,rhomax2
      real*8 r1i,r2i,m1,m2,c1,c2
      real*8 vxi,vyi,omeg
      real*8 epot1i,epot2i,cent 
      real*8 cmmbin,cmxbin,cmybin,cmzbin
      real*8 cmvxbin,cmvybin,cmvzbin
      real*8 x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 Ixx,Iyy,Izz,Ixy,Iyz,Ixz
      real*8 lx12,ly12,lz12,omegax,omegay,omegaz
      real*8 Pcgs,ucgs,TK,rhocgs,scgs,ekincgs
      real*8 vx1c,vy1c,vz1c,x1c,y1c,z1c
      real*8 vxric,vyric,vzric,xric,yric,zric
      real*8 vxricm,vyricm,vzricm,xricm,yricm,zricm
      real*8 xe,xh1,xhe1,xhe2,XXX,YYY
      real*8 vrad,rad,eradcgs
      real*8 lxi,lyi,lzi,lxic,lyic,lzic
      real*8 kappa,opac_rho_t
      external opac_rho_t
      real*8 useeostable
      character*15 filename
      character*5 sph_type, starLabel
      logical core_1,core_2
      integer case_eos, case_run
      common/case_id_all/case_eos, case_run
      common/sep_max/sep,m1,m2
      common/rho_max/idmax1,idmax2,rhomax1,rhomax2
      common/l_mom/lx12,ly12,lz12
      common/cofm_bin/cmmbin,cmxbin,cmybin,cmzbin,cmvxbin,cmvybin,cmvzbin 
      common/bin_par/x1,x2,y1,y2,z1,z2,vx1,vx2,vy1,vy2,vz1,vz2
      common/abundances/XXX,YYY

c     creating the output file
      if(nout.lt.10000) then 
         write(filename,200)nout
      else
         write(filename,201)nout
      endif
 200  format('prof_',I4.4,'.dat')
 201  format('prof_',I5.5,'.dat')
      open(10,file=trim(filename)) 
      
c     omega is not zero only for a corotating framce
      if(nrelax.ge.2) then 
         omeg = sqrt(omega2)
      else
         omeg = 0.d0
      endif
      
      rhomax1=-1d30
      rhomax2=-1d30 
      cmmunb = 0.d0
      cmxunb = 0.d0
      cmyunb = 0.d0
      cmzunb = 0.d0 
      cmvxunb = 0.d0
      cmvyunb = 0.d0
      cmvzunb = 0.d0
      cmmb = 0.d0
      cmxb = 0.d0
      cmyb = 0.d0
      cmzb = 0.d0
      cmvxb = 0.d0
      cmvyb = 0.d0
      cmvzb = 0.d0
      core_1 = .false.
      core_2 = .false.      

      x1 = 0.d0
      x2 = 0.d0
      y1 = 0.d0
      y2 = 0.0d0
      z1 = 0.d0
      z2 = 0.0d0
      vx1 = 0.d0
      vx2 = 0.d0
      vy1 = 0.d0
      vy2 = 0.d0
      vz1 = 0.d0
      vz2 =0.d0
 

      
c finds center of mass location and velocities, ids particles
      
      do i=1,ntot
         if(nrelax.gt.2) then 
            vxi = vx(i) - omeg * y(i)
            vyi = vy(i) + omeg * x(i)
         else
            vxi = vx(i)
            vyi = vy(i)
         endif
         
         etest = 0.5d0*m(i)*(vxi**2+vyi**2+vz(i)**2)+m(i)*grpot(i)+m(i)*u(i)
        
         if(u(i).eq.0.d0) then
             write(*,*)'energy for ',i,' is',etest 
          endif
          
          if (etest.gt.0.d0) then
c     unbound particles, their center of mass and velocity
             id(i) = 4
             cmmunb = cmmunb + m(i)
             cmxunb = cmxunb + x(i)*m(i)
             cmyunb = cmyunb + y(i)*m(i)
             cmzunb = cmzunb + z(i)*m(i)
             cmvxunb = cmvxunb + vxi*m(i)
             cmvyunb = cmvyunb + vyi*m(i)
             cmvzunb = cmvzunb + vz(i)*m(i)
         else
             if (cc(i).eq.cc(1)) then 
                id(i) = 1
                if (u(i).eq.0.d0) then
                   write(*,*)'Core point in ',1, m(i)
                   idmax1 = i
                   core_1 = .true.
                else
                   if(rhomax1.lt.rho(i)) then
                     idmax1 = i
                     rhomax1 = rho(i)
                   endif
                endif
             else
                id(i) = 2
                if (u(i).eq.0.d0) then
                   write(*,*)'Core point in ',2, m(i)
                   idmax2 = i
                   core_2 = .true.
                else
                   if(rhomax2.lt.rho(i)) then
                      idmax2 = i
                      rhomax2 = rho(i)
                   endif
                endif
             endif
             cmmb = cmmb + m(i)
             cmxb = cmxb + x(i)*m(i)
             cmyb = cmyb + y(i)*m(i)
             cmzb = cmzb + z(i)*m(i)
             cmvxb = cmvxb + vxi*m(i)
             cmvyb = cmvyb + vyi*m(i)
             cmvzb = cmvzb + vz(i)*m(i)        
         endif
      enddo
      
      if (cmmunb.gt.0.d0) then 
         cmxunb = cmxunb / cmmunb
         cmyunb = cmyunb / cmmunb
         cmzunb = cmzunb / cmmunb
         cmvxunb = cmvxunb / cmmunb
         cmvyunb = cmvyunb / cmmunb
         cmvzunb = cmvzunb / cmmunb
      endif
      if (cmmb.gt.0.d0) then
         cmxb = cmxb / cmmb
         cmyb = cmyb / cmmb
         cmzb = cmzb / cmmb
         cmvxb = cmvxb / cmmb
         cmvyb = cmvyb / cmmb
         cmvzb = cmvzb / cmmb
      endif

      write(*,*) t*tunit,cmmunb,cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb, 'unbmat'
      write(*,*) t*tunit,cmmb,cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,
     &     "boundmat"


      
      sep=0
      if(case_run.ne.1.and.core_1.and.core_2) then
         if(core_1) idmax1 = 1
         if(core_2) idmax2 = ntot
         sep = sqrt((x(idmax1)-x(idmax2))**2+(y(idmax1)-y(idmax2))**2
     &        +(z(idmax1)-z(idmax2))**2)
c         write(70,*)nout,t*tunit,sep
         c1 = m(idmax1)
         c2 = m(idmax2)
         
         
         do i=1,ntot
            rcmi = sqrt((x(i)-cmxb)**2+(y(i)-cmyb)**2+(z(i)-cmzb)**2)

            if (rcmi.lt.sep .and. id(i).le.2) then
               if((i.ne.idmax1).and.(i.ne.idmax2)) then 
                  r1i = sqrt((x(i)-x(idmax1))**2+(y(i)-y(idmax1))**2+(z(i)-z(idmax1))**2)
                  r2i = sqrt((x(i)-x(idmax2))**2+(y(i)-y(idmax2))**2+(z(i)-z(idmax2))**2)

                  epot1i = -m(idmax1)*m(i)/r1i
                  epot2i = -m(idmax2)*m(i)/r2i
                  if(epot1i.lt.epot2i) then 
                     id(i) = 1
                  else
                     id(i) = 2
                  endif
               endif
            else
               if (id(i).ne.4 .and. (i.ne.idmax1 .and. i.ne.idmax2)) then  
                  id(i) = 3
               endif
            endif
         enddo
      end if
      
      m1 = 0.d0
      m2 = 0.d0
      cmmbin = 0.d0
      cmxbin = 0.d0
      cmybin = 0.d0
      cmzbin = 0.d0
      cmvxbin = 0.d0
      cmvybin = 0.d0
      cmvzbin = 0.d0
      Ixx = 0.d0
      Iyy = 0.d0
      Izz = 0.d0
      Ixy = 0.d0
      Iyz = 0.d0
      Ixz = 0.d0
      lx12 = 0.d0
      ly12 = 0.d0
      lz12 = 0.d0

c     skip for a single star relaxation
      if(case_run.ne.1) then
      do i=1,ntot
         if(nrelax.gt.2) then
           vxi = vx(i) - omeg * y(i)
           vyi = vy(i) + omeg * x(i)
         else
           vxi = vx(i)
           vyi = vy(i)
         endif
         if (id(i).eq.1) then 
            m1 = m1 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vxi*m(i)
            cmvybin = cmvybin + vyi*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
         elseif (id(i).eq.2) then  
            m2 = m2 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vxi*m(i)
            cmvybin = cmvybin + vyi*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
         endif
         if (id(i).le.2) then 
            Ixx = Ixx + m(i)*(y(i)**2+z(i)**2)
            Iyy = Iyy + m(i)*(x(i)**2+z(i)**2)
            Izz = Izz + m(i)*(x(i)**2+y(i)**2)
            Ixy = Ixy - m(i)*x(i)*y(i)
            Iyz = Iyz - m(i)*z(i)*y(i)
            Ixz = Ixz - m(i)*x(i)*z(i)
            lx12 = lx12 + m(i)*(y(i)*vz(i)-z(i)*vyi)
            ly12 = ly12 + m(i)*(z(i)*vxi-x(i)*vz(i))
            lz12 = lz12 + m(i)*(x(i)*vyi-y(i)*vxi)
         endif
      enddo
      if (cmmbin.ne.0.d0) then 
         cmxbin = cmxbin / cmmbin
         cmybin = cmybin / cmmbin
         cmzbin = cmzbin / cmmbin
         cmvxbin = cmvxbin / cmmbin
         cmvybin = cmvybin / cmmbin
         cmvzbin = cmvzbin / cmmbin
      endif


      if(core_1) then 
         x1 = x(idmax1)
         y1 = y(idmax1)
         z1 = z(idmax1)
      end if
      if(core_2) then
         y2 = y(idmax2)
         x2 = x(idmax2)
         z2 = z(idmax2)
      end if
      if(core_1.and.core_2) then
         if(nrelax.gt.2) then   
            vx1 = vx(idmax1)-omeg*y(idmax1)
            vx2 = vx(idmax2)-omeg*y(idmax2)    
            vy1 = vy(idmax1)+omeg*x(idmax1)
            vy2 = vy(idmax2)+omeg*x(idmax2)
         else
            vx1 = vx(idmax1)
            vx2 = vx(idmax2)
            vy1 = vy(idmax1)
            vy2 = vy(idmax2)
         endif
         vz1 = vz(idmax1)
         vz2 = vz(idmax2)
         
         write(*,*)'I tensor',Ixx,Iyy,Izz,Ixy,Ixz,Iyz
         write(*,*)'L vec',lx12,ly12,lz12
         call omega_vec(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx12,ly12,lz12,omegax,omegay,omegaz)
         write(*,*)'Omega (x,y,z)',omegax/tunit,omegay/tunit,omegaz/tunit
         write(*,*)t*tunit,m1,m2,cmmbin,cmxbin,cmybin,cmzbin,
     &        cmvxbin,cmvybin,cmvzbin
      end if
      end if

      if(case_run.eq.0) then 
         call phi(phi1,xL1,m1,m2,sep)
         write(*,*)'m1,m2,sep',m1,m2,sep 
         write(*,*)'Potential,xL1: ',phi1,xL1
         if ((y2.gt.0.d0).and.(x2.gt.0.d0)) then 
            theta = atan(y2/x2)
         elseif ((y2.gt.0.d0).and.(x2.lt.0.d0)) then
            theta = pi + atan(y2/x2)
         elseif ((y2.lt.0.d0).and.(x2.lt.0.d0)) then
            theta = pi + atan(y2/x2)
         elseif ((y2.lt.0.d0).and.(x2.gt.0.d0)) then
            theta = 2.d0*pi + atan(y2/x2)
         endif
         xL1p = xL1*cos(theta) + cmxbin
         yL1p = xL1*sin(theta) + cmybin
         zL1p = cmzbin
         write(*,*)'xL1 prime (x,y)',xL1p,yL1p,theta
         write(*,*)'x1 (x,y,z)(vx,vy,vz)',x1,y1,z1,vx1,vy1,vz1
         write(*,*)'x2 (x,y,z)(vx,vy,vz)',x2,y2,z2,vx2,vy2,vz2
         r1 = sqrt((x1-xL1p)**2+(y1-yL1p)**2+(z1-zL1p)**2)
         r2 = sqrt((x2-xL1p)**2+(y2-yL1p)**2+(z2-zL1p)**2)
         write(*,*)'Spheres 1,2',r1,r2
      end if

         
CC Recalculate I, L, omega, and centre of mass for the particles in 
C the sphere r1 and r2 which gives a better approximation of the
C particles around m1 and m2
      m1 = 0.d0
      m2 = 0.d0
      cmmbin = 0.d0
      cmxbin = 0.d0
      cmybin = 0.d0
      cmzbin = 0.d0
      cmvxbin = 0.d0
      cmvybin = 0.d0
      cmvzbin = 0.d0
      Ixx = 0.d0
      Iyy = 0.d0
      Izz = 0.d0
      Ixy = 0.d0
      Iyz = 0.d0
      Ixz = 0.d0
      lx12 = 0.d0
      ly12 = 0.d0
      lz12 = 0.d0
      do i=1,ntot
         if(nrelax.gt.2) then
           vxi = vx(i) - omeg * y(i)
           vyi = vy(i) + omeg * x(i)
         else
           vxi = vx(i)
           vyi = vy(i)
         endif
       if(i.ne.idmax1 .and. i.ne.idmax2) then 
         r1i = sqrt((x(i)-x1)**2+(y(i)-y1)**2+(z(i)-z1)**2)
         r2i = sqrt((x(i)-x2)**2+(y(i)-y2)**2+(z(i)-z2)**2)     
         etest = 0.5d0*m(i)*(vxi**2+vyi**2+vz(i)**2)+m(i)*grpot(i)+m(i)*u(i)
         if (r1i.lt.r1 .and. etest.lt.0.d0) then 
             id(i)=1
         elseif(r2i.lt.r2 .and. etest.lt.0.d0) then 
             id(i)=2
         elseif(etest .lt. 0.d0) then 
             id(i)=3
         else
             id(i)=4
         endif
       endif
       if (id(i).eq.1) then
            m1 = m1 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vxi*m(i)
            cmvybin = cmvybin + vyi*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
         elseif (id(i).eq.2) then
            m2 = m2 + m(i)
            cmmbin = cmmbin + m(i)
            cmxbin = cmxbin + x(i)*m(i)
            cmybin = cmybin + y(i)*m(i)
            cmzbin = cmzbin + z(i)*m(i)
            cmvxbin = cmvxbin + vxi*m(i)
            cmvybin = cmvybin + vyi*m(i)
            cmvzbin = cmvzbin + vz(i)*m(i)
       endif
      enddo
      
      if (cmmbin.ne.0.d0) then
          cmxbin = cmxbin / cmmbin
          cmybin = cmybin / cmmbin
          cmzbin = cmzbin / cmmbin
          cmvxbin = cmvxbin / cmmbin
          cmvybin = cmvybin / cmmbin
          cmvzbin = cmvzbin / cmmbin
      endif

      write(*,*)t*tunit,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,'bef'

      x1c = x(1)
      y1c = y(1)
      z1c = z(1)
      vx1c = vx(1) - omeg * y1c 
      vy1c = vy(1) + omeg * x1c
      vz1c = vz(1)
      
      write(20,*) nout,t*tunit,sep,x1c,y1c,z1c,vx1c*vunit,vy1c*vunit,vz1c*vunit

      write(*,*)"fraction of mass",XXX,YYY, nrelax
      
      do i=1,ntot
         if(nrelax.ge.2) then
           vxi = vx(i) - omeg * y(i)
           vyi = vy(i) + omeg * x(i)
         else
           vxi = vx(i)
           vyi = vy(i)
         endif
         rhocgs = rho(i)*Msun/Rsun**3.d0
         ucgs = u(i)*G*Msun/Rsun
         vxric =  vxi - vx1c
         vyric =  vyi - vy1c
         vzric = vz(i)- vz1c
         xric = x(i) - x1c
         yric = y(i) - y1c
         zric = z(i) - z1c
         vxricm =  vxi - cmvxbin
         vyricm =  vyi - cmvybin
         vzricm = vz(i)- cmvzbin
         xricm = x(i) - cmxbin
         yricm = y(i) - cmybin
         zricm = z(i) - cmzbin
         rad = sqrt(xric**2+yric**2+zric**2)
         vrad = (xric*vxric+yric*vyric+zric+vzric)/rad
         ekincgs = 0.5*(vxi**2+vyi**2+vz(i)**2)*eunitm
         eradcgs = m(i)*vrad**2*eunit
         lxi = m(i)*(yric*vzric-zric*vyric)
         lyi = m(i)*(zric*vxric-xric*vzric)
         lzi = m(i)*(xric*vyric-yric*vxric)
         lxic = m(i)*(yricm*vzricm-zricm*vyricm)
         lyic = m(i)*(zricm*vxricm-xricm*vzricm)
         lzic = m(i)*(xricm*vyricm-yricm*vxricm)
         if(case_eos.eq.0) then 
            if(u(i).gt.0.d0) then
               Pcgs = useeostable(ucgs,rhocgs,3)
               TK = useeostable(ucgs,rhocgs,1)
               scgs = useeostable(ucgs,rhocgs,4)              
               call deg_of_ion(rhocgs,TK,XXX,YYY,xe,xh1,xhe1,xhe2,iter)
               kappa = opac_rho_t(rhocgs,TK)
            else
               Pcgs = 0.d0
               TK = 0.d0
               scgs = 0.d0
               xh1=1.d0
               xhe1=0.d0
               xhe2=1.d0
               xe=(XXX*xh1+YYY/4.0d0*(xhe1+2.d0*xhe2))/(XXX+YYY/4.0d0)
               iter=0
               kappa=0.d0
            endif 
            if(kappa.ne.0.d0) kappa=10**(kappa)
         else
            if(u(i).gt.0.d0) then
               call getTemperature(qconst*rhocgs/(mu(i)*mp),
     $              -ucgs*rhocgs/arad,TK)
               scgs = 3.0*log(TK/rhocgs**(2.0/3.0))/(2.0*(mu(i)*mp)/mH)
     $              + 4.0*arad*TK**3.0/(3.0*rhocgs)/boltz/NA
               Pcgs = rhocgs*boltz*TK/(mu(i)*mp) + arad*TK**4/3.d0
            else
               TK = 0.0d0
               scgs = 0.d0
               Pcgs = 0.d0
            endif
         endif


         if(id(i).le.2) sph_type=' bin '
         if(id(i).eq.3) sph_type=' cir '
         if(id(i).eq.4) sph_type=' eje '
         

         if( (1000*(i/1000)).eq.i) write(*,*) "debug", i, id(i), sph_type, m(i)
         
         if (id(i).le.2) then
            
            Ixx = Ixx + m(i)*((y(i)-cmybin)**2+(z(i)-cmzbin)**2)
            Iyy = Iyy + m(i)*((x(i)-cmxbin)**2+(z(i)-cmzbin)**2)
            Izz = Izz + m(i)*((x(i)-cmxbin)**2+(y(i)-cmybin)**2)
            Ixy = Ixy - m(i)*(x(i)-cmzbin)*(y(i)-cmybin)
            Iyz = Iyz - m(i)*(z(i)-cmzbin)*(y(i)-cmybin)
            Ixz = Ixz - m(i)*(x(i)-cmxbin)*(z(i)-cmzbin)
            lx12 = lx12 + m(i)*((y(i)-cmybin)*(vz(i)-cmvzbin)-(z(i)-cmzbin)*(vyi-cmvybin))
            ly12 = ly12 + m(i)*((z(i)-cmzbin)*(vxi-cmvxbin)-(x(i)-cmxbin)*(vz(i)-cmvzbin))
            lz12 = lz12 + m(i)*((x(i)-cmxbin)*(vyi-cmvybin)-(y(i)-cmybin)*(vxi-cmvxbin))

         end if

         if (cc(i).eq.cc(1)) starLabel = " 0 "
         if (cc(i).ne.cc(1)) starLabel = " 1 "
            
         if(case_eos.eq.0) then
c     outputs: 19 doubles and a character
c     if you change here, do chage the reading in sorted.f!
            write(10,100)   xric,yric,zric,m(i),
     &           Pcgs,rhocgs,ucgs,ekincgs*1d15,scgs,
     &           grpot(i)*eunitm*1d15,lzic*lunit,lzi*lunit,
     &           xh1,xhe1,xhe2,
     &           divv(i)/tunits,kappa,TK, h(i), sph_type, starLabel
         else
            write(10,110)xric,yric,zric,m(i),Pcgs,rhocgs,ucgs,ekincgs*1d15,scgs,
     &           grpot(i)*eunitm*1d15,lzic*lunit,lzi*lunit,
     &           divv(i)/tunits,TK,  sph_type, starLabel
         endif
 100  format(19e16.8,A,A)
 110  format(14e16.8,A,A)
 

         if(id(i).eq.1) then 
            x1 = x1 + x(i)*m(i)
            y1 = y1 + y(i)*m(i)
            z1 = z1 + z(i)*m(i)
            vx1 = vx1 + vxi*m(i)
            vy1 = vy1 + vyi*m(i)
            vz1 = vz1 + vz(i)*m(i)
         endif
         if(id(i).eq.2) then
            x2 = x2 + x(i)*m(i)
            y2 = y2 + y(i)*m(i)
            z2 = z2 + z(i)*m(i)
            vx2 = vx2 + vxi*m(i)
            vy2 = vy2 + vyi*m(i)
            vz2 = vz2 + vz(i)*m(i)
         endif
      enddo


      write(*,*) "profile swapped out"
      close(10)
      
      if(m1.gt.0) then
         x1=x1/m1;y1=y1/m1;z1=z1/m1;vx1=vx1/m1;vy1=vy1/m1;vz1=vz1/m1
      end if
      if(m2.gt.0) then
         x2=x2/m2;y2=y2/m2;z2=z2/m2;vx2=vx2/m2;vy2=vy2/m2;vz2=vz2/m2
      end if
         
      write(*,*)t*tunit,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,'aft' 

      if(case_run.eq.0.and.core_1.and.core_2) then 
         write(*,*)'I tensor',Ixx,Iyy,Izz,Ixy,Ixz,Iyz
         write(*,*)'L vec',lx12,ly12,lz12
         call omega_vec(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx12,ly12,lz12,omegax,omegay,omegaz)
         write(*,*)'Omega (x,y,z)',omegax/tunit,omegay/tunit,omegaz/tunit
         write(*,*)t*tunit,m1,m2,cmmbin,cmxbin,cmybin,cmzbin,
     &        cmvxbin*vunit,cmvybin*vunit,cmvzbin*vunit,' bin_evol'
         write(*,*)t*tunit,2.d0*pi/(sqrt(omegax**2+omegay**2+omegaz**2))*tunit,
     &        'Period'
         write(*,*)t*tunit,m1,x1-cmxbin,y1-cmybin,z1-cmzbin,
     &        (vx1-cmvxbin)*vunit,(vy1-cmvybin)*vunit,(vz1-cmvzbin)*vunit,
     &        ' star_1'
         write(*,*)t*tunit,m2,x2-cmxbin,y2-cmybin,z2-cmzbin,
     &        (vx2-cmvxbin)*vunit,(vy2-cmvybin)*vunit,(vz2-cmvzbin)*vunit,
     &        ' star_2'
      end if

      
 42   continue
      end
