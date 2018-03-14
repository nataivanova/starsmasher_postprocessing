c we want to have all the variables explicitly declared
      implicit none
c Define the maximum number of particles
      integer nmax
      parameter(nmax=1000000)
c Definition of StarSmasher variables, these are to match up with the variables names
      integer ntot,nnopt,nout,nit,nav,ngr,nrelax,cc(nmax)
      real*8 hmin,hmax,sep0,tf,dtout,t,alpha,beta,
     &       tskipahead,trelax,dt,omega2
      real*8 x(nmax),y(nmax),z(nmax),m(nmax),h(nmax),rho(nmax),vx(nmax),vy(nmax),vz(nmax),
     &       vxdot(nmax),vydot(nmax),vzdot(nmax),u(nmax),udot(nmax),grpot(nmax),
     &       mu(nmax),divv(nmax) 
      integer id(nmax)
      real*8 cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,cmmb
      real*8 cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb,cmmunb

c Common variables
      common/nconst/ntot,nnopt,nout,nit,nav,ngr,nrelax
      common/init_real/hmin,hmax,sep0,tf,dtout,t,alpha,beta,
     &       tskipahead,trelax,dt,omega2      
      common/dynamics/m,x,y,z,vx,vy,vz,vxdot,vydot,vzdot,h,grpot,divv
      common/thermodynamics/rho,u,udot,mu,cc
      common/id_tracks/id
      common/bound_ma/cmxb,cmyb,cmzb,cmvxb,cmvyb,cmvzb,cmmb
      common/unb_ma/cmxunb,cmyunb,cmzunb,cmvxunb,cmvyunb,cmvzunb,cmmunb

c Constants (note that all constant are in cgs units)
      real*8 pi,mp,G,Msun,Rsun,eunit,tunit,eunitm
      real*8 tunits,lunits,lunit,rhounit,vunit
      real*8 kelvin,gram,sec,cm,erg,boltz,crad,planck,crad2,sigma,arad,
     &       qconst,punit
      real*8 mH,NA
c defining fundamental constants
      parameter(pi=3.1415926535897932384626d0)
      parameter(mp=1.67262158d-24,G=6.67390d-8)
      parameter(Msun=1.9891d33,Rsun=6.9599d10)
      parameter(mH=1.67262158d-24)
      parameter(NA=6.0221417d23)

c definining  derived constants
      parameter(eunitm = G*Msun/Rsun*1.0d-15)
      parameter(eunit = G*Msun**2/Rsun*1.0d-46)
      parameter(tunits = sqrt(Rsun**3/G/Msun))
      parameter(tunit = tunits/60.d0/60.d0/24.d0) !in days
      parameter(lunits = Rsun**2/tunits)
      parameter(lunit = Msun*Rsun**2/tunits*1.0d-52)
      parameter(rhounit = Msun/Rsun**3)
      parameter(vunit = sqrt(G*Msun/Rsun)*1e-5)

      parameter(gram=1.d0,sec=1.d0,cm=1.d0,kelvin=1.d0)
      parameter(erg = gram * cm**2/sec**2)
      parameter(boltz = 1.380658d-16*erg/kelvin)
      parameter(crad = 2.997924580d+10*cm/sec,
     &          planck = 6.6260755d-27*gram*cm**2/sec)
      parameter(crad2 = crad**2,
     &          sigma = pi**2*boltz*(boltz*2.0d0*pi/planck)**3/60.0d0/crad2,
     &          arad  = 4.0d0*sigma/crad,
     &          qconst = 1.5d0*boltz/arad)
