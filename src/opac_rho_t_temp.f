      real*8 function opac_rho_t(rho,t)
      implicit none
      real*8 rho,T,T6,R,logRi,logTi
      integer nummaxR,nummaxT,i
      integer id_T,id_R
      real*8 R1,R2,Q1,Q2,Q3,Q4,x1,x2,y1,y2
      real*8 x,y
      parameter(nummaxR=19,nummaxT=139)
      real*8 logR(nummaxR),logT(nummaxT)
      real*8 logopac(nummaxR,nummaxT)
      common/opac_low/logT,logR,logopac
      T6 = T*1.d-6
      r = rho/(T6**3.d0)
      logRi = log10(r)
      logTi = log10(T)
     
      x = logTi
      y = logRi      

      if (logRi.lt.-8.d0) then 
          id_r = 1
      else
          
          do i=1,nummaxR-1
             if ((logR(i).le.logRi).and.(logRi.le.logR(i+1))) then
                id_r = i
                exit
             else
                id_r = -1
             endif
          enddo
      endif
      do i=1,nummaxT-1
         if ((logT(i).le.logTi).and.(logTi.le.logT(i+1))) then
            id_T = i
            exit
         else
            id_T = -1
         endif
      enddo
      if ((id_r.ne.-1).and.(id_t.ne.-1)) then 
          x1 = logT(id_T)
          x2 = logT(id_T+1)
          y1 = logR(id_r)
          y2 = logR(id_r+1)
          Q1 = logopac(id_r,id_T)
          Q2 = logopac(id_r,id_T+1)
          Q3 = logopac(id_r+1,id_T)
          Q4 = logopac(id_r+1,id_T+1)
          R1 = Q2*(x-x1)/(x2-x1)+Q1*(x2-x)/(x2-x1)
          R2 = Q4*(x-x1)/(x2-x1)+Q3*(x2-x)/(x2-x1)
          opac_rho_t = R2*(y-y1)/(y2-y1)+R1*(y2-y)/(y2-y1)
          !write(*,*)x1,x2,y1,y2,Q1,Q2,Q3,Q4,R1,R2 
      !opac_rho_t=1.d0
      else
          opac_rho_t=0.d0
      endif
  
      !write(69,*)id_r,id_t,rho,T,logRi,logTi,opac_rho_t
 
      return
      end
