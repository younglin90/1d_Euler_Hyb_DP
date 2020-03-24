
module modFaceValues

    use modValues
    implicit none
    
    contains
    
    subroutine getFaceValues(ReconMethod,FluxMethod,&
        calfgP,calfgN,calfUx,calfCM1,calfCM2,calfRho,calfDrDp)
        implicit none
        character(len=*),optional,intent(in) :: &
            calfgP,calfgN,calfUx,calfCM1,calfCM2,calfRho,calfDrDp
        
        character(len=*) :: ReconMethod, FluxMethod
        
        
        if( trim(ReconMethod) == '1stUpwind' ) then
            CALL FirUpwind
        elseif( trim(ReconMethod) == '2ndUpwind' ) then
            CALL SecUpwind
        elseif( trim(ReconMethod) == '2ndCentral' ) then
            CALL SecCentral
        elseif( trim(ReconMethod) == '3thMUSCL' ) then
            CALL ThrMUSCL
        else
            print*,'  ReconMethod not defined'
            stop
        endif
        
        if( trim(FluxMethod) == 'AUSM' ) then
            CALL calAUSM(calfgP,calfgN,calfUx,calfCM1,calfCM2,calfRho,calfDrDp)
        endif
        
        
    endsubroutine

    !================================
    subroutine FirUpwind
        implicit none
        
        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(i,1) = cRho(lCell) 
            wL(i,2) = cU(lCell) 
            wL(i,3) = cP(lCell)
            
            wR(i,1) = cRho(rCell) 
            wR(i,2) = cU(rCell) 
            wR(i,3) = cP(rCell)
        enddo
        
    endsubroutine
    

    !================================
    subroutine SecUpwind
        implicit none
        
        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(i,1) = cRho(lCell) + 0.5d0*(cRho(rCell)-cRho(lCell))
            wL(i,2) = cU(lCell) + 0.5d0*(cU(rCell)-cU(lCell))
            wL(i,3) = cP(lCell) + 0.5d0*(cP(rCell)-cP(lCell))
            
            wR(i,1) = cRho(rCell) - 0.5d0*(cRho(rCell)-cRho(lCell))
            wR(i,2) = cU(rCell) - 0.5d0*(cU(rCell)-cU(lCell))
            wR(i,3) = cP(rCell) - 0.5d0*(cP(rCell)-cP(lCell))
        enddo
        
    endsubroutine
    

    !================================
    subroutine SecCentral
        implicit none
        
        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(i,1) = 0.5d0*(cRho(lCell)+cRho(rCell))
            wL(i,2) = 0.5d0*(cU(lCell)+cU(rCell))
            wL(i,3) = 0.5d0*(cP(lCell)+cP(rCell))
            
            wR(i,1) = 0.5d0*(cRho(lCell)+cRho(rCell))
            wR(i,2) = 0.5d0*(cU(lCell)+cU(rCell))
            wR(i,3) = 0.5d0*(cP(lCell)+cP(rCell))
        enddo
        
    endsubroutine
    
    !================================
    subroutine ThrMUSCL
        implicit none
        integer :: lmCell,rpCell
        
        do i=nFaceStr+1,nFaceEnd-1
            lmCell = i-1; lCell = i; rCell = i+1; rpCell = i+2
            CALL MUSCL_3th(cRho(lmCell),cRho(lCell),cRho(rCell),cRho(rpCell),&
                wL(i,1),wR(i,1))
            CALL MUSCL_3th(cU(lmCell),cU(lCell),cU(rCell),cU(rpCell),&
                wL(i,2),wR(i,2))
            CALL MUSCL_3th(cP(lmCell),cP(lCell),cP(rCell),cP(rpCell),&
                wL(i,3),wR(i,3))
        enddo
        
        i=nFaceStr
        lmCell = i; lCell = i; rCell = i+1; rpCell = i+2
        CALL MUSCL_3th(cRho(lmCell),cRho(lCell),cRho(rCell),cRho(rpCell),&
            wL(i,1),wR(i,1))
        CALL MUSCL_3th(cU(lmCell),cU(lCell),cU(rCell),cU(rpCell),&
            wL(i,2),wR(i,2))
        CALL MUSCL_3th(cP(lmCell),cP(lCell),cP(rCell),cP(rpCell),&
            wL(i,3),wR(i,3))
        
        i=nFaceEnd
        lmCell = i-1; lCell = i; rCell = i+1; rpCell = i+1
        CALL MUSCL_3th(cRho(lmCell),cRho(lCell),cRho(rCell),cRho(rpCell),&
            wL(i,1),wR(i,1))
        CALL MUSCL_3th(cU(lmCell),cU(lCell),cU(rCell),cU(rpCell),&
            wL(i,2),wR(i,2))
        CALL MUSCL_3th(cP(lmCell),cP(lCell),cP(rCell),cP(rpCell),&
            wL(i,3),wR(i,3))
        
    endsubroutine
    
    
    !================================
    subroutine calAUSM(calfgP,calfgN,calfUx,calfCM1,calfCM2,calfRho,calfDrDp)
        implicit none
        character(len=*),optional,intent(in) :: &
            calfgP,calfgN,calfUx,calfCM1,calfCM2,calfRho,calfDrDp
        
        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            CALL calfaceValAUSM(wL(i,:),wR(i,:),eosGamma, fMdot(i),fU(i),fP(i))
            fMdot(i) = fMdot(i)*area(i)
            if(present(calfgP)) fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            if(present(calfgN)) fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
            if(present(calfUx)) fUx(i) = fgP(i)*cU(lCell) + fgN(i)*cU(rCell)
            if(present(calfCM1)) fcoeffMom1(i) = fgP(i)*coeffMom1(lCell) + fgN(i)*coeffMom1(rCell)
            if(present(calfCM2)) fcoeffMom2(i) = fgP(i)*coeffMom2(lCell) + fgN(i)*coeffMom2(rCell)
            if(present(calfRho)) fRho(i) = fgP(i)*cRho(lCell) + fgN(i)*cRho(rCell)
            if(present(calfDrDp)) fDrDp(i) = fgP(i)*cDrDp(lCell) + fgN(i)*cDrDp(rCell)
            
        enddo
        
        
    endsubroutine
    


    subroutine calfaceValAUSM(wL,wR,gamma, mLR,uLR,pLR)

        implicit none
    
        real(8), dimension(3), intent(in) :: wL, wR
        real(8), intent(in) :: gamma
        real(8), intent(out):: mLR,uLR,pLR
        real(8) :: rho_L,rho_R,u_L,u_R,p_L,p_R,c_L,c_R, &
        H_L,H_R,chat,M_L,M_R,MLP,preP,MRM,preM,mdot
    
        ! AUSM+ scheme
        ! wL : left value of cell interface (face)
        ! wR : right
        rho_L = wL(1); rho_R = wR(1)
        u_L = wL(2); u_R = wR(2)
        p_L = wL(3); p_R = wR(3)
    
        c_L = sqrt(gamma*p_L/rho_L); c_R = sqrt(gamma*p_R/rho_R)
        H_L = gamma/(gamma-1.d0)*p_L/rho_L+0.5d0*u_L**2.d0
        H_R = gamma/(gamma-1.d0)*p_R/rho_R+0.5d0*u_R**2.d0
        !----------------------------------------------------
    
        chat = 0.5d0*(c_L+c_R)
        M_L = u_L/chat; M_R = u_R/chat
    
        !------ numerical mass flux & pressure flux in AUSM+ scheme ------
        if( abs(M_L) > 1.0 ) then
            MLP = 0.5*(M_L+abs(M_L))
            preP = 0.5*(1.0+sign(1.0,M_L))
        else
            MLP = 0.25*(M_L+1.0)**2.0
            preP = 0.25*(M_L+1.0)**2.0*(2.0-M_L)
        endif
    
        if( abs(M_R) > 1.0 ) then
            MRM = 0.5*(M_R-abs(M_R))
            preM = 0.5*(1.0-sign(1.0,M_R))
        else
            MRM = -0.25*(M_R-1.0)**2.0
            preM = 0.25*(M_R-1.0)**2.0*(2.0+M_R)
        endif
    
    

        uLR = (MLP + MRM)*chat
    
        if( uLR>=0.d0 ) then
            mLR = rho_L*uLR
        else
            mLR = rho_R*uLR
        endif
    
        pLR = preP*p_L+preM*p_R



    endsubroutine

    
    endmodule
    
    
    
    
    
        subroutine MUSCL_3th(qm,qi,qp,qp2,ql,qr)
      
            implicit none

            real(8)  qm,qi,qp,qp2,ql,qr
            real(8)  kpa,phi4,gbq,dbq,kam,kap,dq1,dq2,dq3,betat
      
            kpa = 1.0d0/3.0d0
            phi4  = 1.0d0*.250d0
  
            !betat=(3.0d0-kpa)/(1.0d0-kpa)
            betat=1.0d0
  
            kap   = 1.0d0+kpa
            kam   = 1.0d0-kpa      

            dq1 = ( qi - qm)
            dq2 = ( qp - qi)
            dq3 = (qp2 - qp)

            gbq = minmod(betat*dq1,dq2)
            dbq = minmod(betat*dq2,dq1)
            ql  = qi+phi4*(kam*gbq+kap*dbq)

            gbq = minmod(dq2,dq3)
            dbq = minmod(dq3,dq2)
            qr  = qp-phi4*(kap*gbq+kam*dbq)
          
        contains
            function minmod(xx,yy,zz) result(mu)
                implicit none
                real(8), intent(in) :: xx,yy  !Input
                real(8),optional, intent(in) :: zz
                real(8)             :: mu        !output
                real(8) :: aa,bb,cc,dd,ee

                if(present(zz)) then 
                    aa = dsign(1.0d0,xx)
                    bb = dabs(xx)
                    cc = aa*yy
                    dd = aa*zz
                    ee = dmin1(bb,cc,dd)
                    mu = aa*dmax1(0.d0,ee)
                else
                    aa = dsign(1.0d0,xx)
                    bb = dabs(xx)
                    cc = aa*yy
                    dd = dmin1(bb,cc)
                    mu = aa*dmax1(0.d0,dd)      
                endif
            endfunction

        end subroutine
    