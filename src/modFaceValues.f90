
module modFaceValues

    use modValues
    use modEOS
    implicit none
    
    contains
    
    subroutine getFaceValues(&
        FaceMethodU,FaceMethodP,&
        calfUx,calfCM1,calfCM2,calfRho,calfDrDp)
        implicit none
        character(len=*),optional,intent(in) :: &
            FaceMethodU,FaceMethodP,&
            calfUx,calfCM1,calfCM2,calfRho,calfDrDp
        integer :: lmCell,rpCell
        real(8) :: leftcRho,leftcC,leftcHt,rightcRho,rightcC,rightcHt,eosY(nspe)
        
        
        do i=nFaceStr,nFaceEnd
            lmCell = max(nFaceStr,i-1); lCell = i
            rCell = i+1; rpCell = min(nFaceEnd,i+2)
            
            CALL getRecon(chReconP,cP(lmCell),cP(lCell),cP(rCell),cP(rpCell),&
                        wL(i,1),wR(i,1))
            CALL getRecon(chReconU,cU(lmCell),cU(lCell),cU(rCell),cU(rpCell),&
                        wL(i,2),wR(i,2))
            CALL getRecon(chReconT,cT(lmCell),cT(lCell),cT(rCell),cT(rpCell),&
                        wL(i,3),wR(i,3))
            do nq=1,nspe-1
                CALL getRecon(chReconY,cY(lmCell,nq),cY(lCell,nq),&
                    cY(rCell,nq),cY(rpCell,nq),&
                    wL(i,3+nq),wR(i,3+nq))
            enddo
            
            forall(nq=1:nspe-1) wL(i,3+nq) = dmax1( 0.d0 , dmin1( 1.d0, wL(i,3+nq) ) )
            wL(i,3+nspe) = 1.d0 - sum( wL(i,4:2+nspe) )
            wL(i,3+nspe) = dmax1( 0.d0 , dmin1( 1.d0, wL(i,3+nspe) ) )
            forall(nq=1:nspe-1) wR(i,3+nq) = dmax1( 0.d0 , dmin1( 1.d0, wR(i,3+nq) ) )
            wR(i,3+nspe) = 1.d0 - sum( wR(i,4:2+nspe) )
            wR(i,3+nspe) = dmax1( 0.d0 , dmin1( 1.d0, wR(i,3+nspe) ) )
            
            !> get Rho, C, Ht
            CALL getRCH(wL(i,1),wL(i,2),wL(i,3),wL(i,4:3+nspe), leftcRho,leftcC,leftcHt)
            CALL getRCH(wR(i,1),wR(i,2),wR(i,3),wR(i,4:3+nspe), rightcRho,rightcC,rightcHt)
            !
            !> face U, Mdot
            if( trim(FaceMethodU) == 'AUSM' ) then
                CALL calfaceVelAUSM(&
                    wL(i,1),wL(i,2),wL(i,3),leftcRho,leftcC,&
                    wR(i,1),wR(i,2),wR(i,3),rightcRho,rightcC, fU(i))
            elseif( trim(FaceMethodU) == 'SLAU' ) then
                CALL calfaceVelSLAU(&
                    wL(i,1),wL(i,2),wL(i,3),leftcRho,leftcC,&
                    wR(i,1),wR(i,2),wR(i,3),rightcRho,rightcC, fU(i))
            else
                print*,'  ReconMethod not defined'
                stop
            endif
            !> upwind coefficient
            fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
            !> Mdot, upwind
            fMdot(i) = leftcRho*fU(i)*fgP(i)*area(i) &
                     + rightcRho*fU(i)*fgN(i)*area(i)
            
            !> face P
            if( present(FaceMethodP) ) then
                if( trim(FaceMethodP) == 'AUSM' ) then
                    CALL calFacePressAUSM(&
                        wL(i,1),wL(i,2),wL(i,3),leftcRho,leftcC,&
                        wR(i,1),wR(i,2),wR(i,3),rightcRho,rightcC, fP(i))
                elseif( trim(FaceMethodP) == 'SLAU2' ) then
                    CALL calFacePressSLAU(&
                        wL(i,1),wL(i,2),wL(i,3),leftcRho,leftcC,&
                        wR(i,1),wR(i,2),wR(i,3),rightcRho,rightcC, fP(i))
                elseif( trim(FaceMethodP) == 'SecOrder' ) then
                    CALL calFacePressSecOrder(&
                        wL(i,1),wL(i,2),wL(i,3),leftcRho,leftcC,&
                        wR(i,1),wR(i,2),wR(i,3),rightcRho,rightcC, fP(i))
                else
                    print*,'  ReconMethod not defined'
                    stop
                endif
            endif
            
            
            if(present(calfUx)) &
                fUx(i) = fgP(i)*cU(lCell) + fgN(i)*cU(rCell)
            if(present(calfCM1)) &
                fcoeffMom1(i) = fgP(i)*coeffMom1(lCell) + fgN(i)*coeffMom1(rCell)
            if(present(calfCM2)) &
                fcoeffMom2(i) = fgP(i)*coeffMom2(lCell) + fgN(i)*coeffMom2(rCell)
            if(present(calfRho)) &
                fRho(i) = fgP(i)*cRho(lCell) + fgN(i)*cRho(rCell)
            if(present(calfDrDp)) &
                fDrDp(i) = fgP(i)*cDrDp(lCell) + fgN(i)*cDrDp(rCell)
            
        enddo
    endsubroutine

    !================================
    subroutine getRecon(ReconMethod,cm2Phi,cm1Phi,cp1Phi,cp2Phi, leftPhi,rightPhi)
        implicit none
        character(len=*) :: ReconMethod
        real(8) :: cm2Phi,cm1Phi,cp1Phi,cp2Phi
        real(8) :: leftPhi,rightPhi
        
        if( trim(ReconMethod) == '1stOrder' ) then
            leftPhi = cm1Phi
            rightPhi = cp1Phi
        elseif( trim(ReconMethod) == '3thMUSCL' ) then
            CALL MUSCL_3th(cm2Phi,cm1Phi,cp1Phi,cp2Phi,leftPhi,rightPhi)
        else
            print*,'  ReconMethod not defined'
            stop
        endif
    
    endsubroutine
        
        
    !================================
    subroutine calFacePressSecOrder(fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR, pLR)

        implicit none
    
        real(8), intent(in) :: fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR
        real(8), intent(out):: pLR
        real(8) :: revAmatMom,cPrevAmatMom,revLAmatMom,revRAmatMom
        integer :: lmCell,rpCell

        lmCell = max(nFaceStr,i-1); lCell = i
        rCell = i+1; rpCell = min(nFaceEnd,i+2)

        
        revLAmatMom = 1.d0/( cDRDP(lCell)/timestep )
        revRAmatMom = 1.d0/( cDRDP(rCell)/timestep )
        revAmatMom = revLAmatMom + revRAmatMom
        cPrevAmatMom = fPL*revLAmatMom + fPR*revRAmatMom
        
        pLR = cPrevAmatMom/revAmatMom
        
    endsubroutine
       
    


    !================================
    subroutine calFaceVelAUSM(fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR, uLR)

        implicit none
    
        real(8) :: fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR
        real(8) :: uLR,M_L,M_R,MLP,MRM,chat
    
        ! AUSM+ scheme
        ! wL : left value of cell interface (face)
        ! wR : right
        !----------------------------------------------------
    
        chat = 0.5d0*(fCL+fCR)
        M_L = fUL/chat; M_R = fUR/chat
    
        !------ numerical mass flux & pressure flux in AUSM+ scheme ------
        if( abs(M_L) > 1.0 ) then
            MLP = 0.5*(M_L+abs(M_L))
        else
            MLP = 0.25*(M_L+1.0)**2.0
        endif
        if( abs(M_R) > 1.0 ) then
            MRM = 0.5*(M_R-abs(M_R))
        else
            MRM = -0.25*(M_R-1.0)**2.0
        endif
    
        uLR = (MLP + MRM)*chat

    endsubroutine
    

    !================================
    subroutine calFaceVelSLAU(fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR, uLR)

        implicit none
    
        real(8) :: fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR
        real(8) :: uLR
        real(8) :: rho_L,rho_R,u_L,u_R,p_L,p_R,c_L,c_R, &
         chat,M_L,M_R,MLP,preP,MRM,preM,RhouLR,Mcy,phi_c,g_c,Mbar,&
         D_L,D_R,D_rho,D_P,KLR
    
        ! AUSM+ scheme
        ! wL : left value of cell interface (face)
        ! wR : right
        !----------------------------------------------------
    
        chat = 0.5d0*(fCL+fCR)
        M_L = fUL/chat; M_R = fUR/chat
    
	    ! !> SLAU
        KLR = dsqrt(0.5d0*(fUL*fUL+fUR*fUR))
	    Mcy = dmin1(1.d0,KLR/chat)
	    phi_c = (1.d0-Mcy)**2.d0 
	    g_c = 1.d0 + dmax1( dmin1(M_L,0.d0), -1.d0 )*dmin1( dmax1(M_R,0.d0), 1.d0 ) 
        Mbar = ( fRhoL*dabs(M_L)+fRhoR*dabs(M_R) ) / ( fRhoL + fRhoR )

        D_L = M_L+(1.d0-g_c)*dabs(M_L) 
        D_R = M_R-(1.d0-g_c)*dabs(M_R)
        D_rho = Mbar*g_c
        D_P = 1.d0/chat**2.d0*phi_c

        RhouLR = 0.5d0*chat*(D_L*fRhoL+D_R*fRhoR-D_rho*(fRhoR-fRhoL)-D_P*(fPR-fPL))
	    if( RhouLR >= 0.d0 ) then
            uLR = RhouLR/fRhoL
        else
            uLR = RhouLR/fRhoR
        endif

    endsubroutine
    
    !================================
    subroutine calFacePressAUSM(fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR, pLR)

        implicit none
    
        real(8), intent(in) :: fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR
        real(8), intent(out):: pLR
        real(8) :: rho_L,rho_R,u_L,u_R,p_L,p_R,c_L,c_R, &
        H_L,H_R,chat,M_L,M_R,MLP,preP,MRM,preM,mdot
    
        ! AUSM+ scheme
        ! wL : left value of cell interface (face)
        ! wR : right
        !----------------------------------------------------
    
        chat = 0.5d0*(fCL+fCR)
        M_L = fUL/chat; M_R = fUR/chat
    
        !------ numerical mass flux & pressure flux in AUSM+ scheme ------
        if( abs(M_L) > 1.0 ) then
            preP = 0.5*(1.0+sign(1.0,M_L))
        else
            preP = 0.25*(M_L+1.0)**2.0*(2.0-M_L)
        endif
        if( abs(M_R) > 1.0 ) then
            preM = 0.5*(1.0-sign(1.0,M_R))
        else
            preM = 0.25*(M_R-1.0)**2.0*(2.0+M_R)
        endif
        pLR = preP*fPL+preM*fPR

    endsubroutine

    
    !================================
    subroutine calFacePressSLAU(fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR, pLR)

        implicit none
    
        real(8), intent(in) :: fPL,fUL,fTL,fRhoL,fCL,fPR,fUR,fTR,fRhoR,fCR
        real(8), intent(out):: pLR
        real(8) :: rho_L,rho_R,u_L,u_R,p_L,p_R,c_L,c_R, &
        H_L,H_R,chat,M_L,M_R,MLP,preP,MRM,preM,mdot,KLR,rhohat
    
        ! AUSM+ scheme
        ! wL : left value of cell interface (face)
        ! wR : right
        !----------------------------------------------------
    
        chat = 0.5d0*(fCL+fCR)
        M_L = fUL/chat; M_R = fUR/chat
    
        !------ numerical mass flux & pressure flux in AUSM+ scheme ------
        if( abs(M_L) > 1.0 ) then
            preP = 0.5*(1.0+sign(1.0,M_L))
        else
            preP = 0.25*(M_L+1.0)**2.0*(2.0-M_L)
        endif
        if( abs(M_R) > 1.0 ) then
            preM = 0.5*(1.0-sign(1.0,M_R))
        else
            preM = 0.25*(M_R-1.0)**2.0*(2.0+M_R)
        endif
        KLR = dsqrt(0.5d0*(fUL*fUL+fUR*fUR))
        rhohat = 0.5d0*(fRhoL+fRhoR)
	    PLR = 0.5d0*(fPL+fPR) &
              + 0.5d0*(fPL-fPR)*(preP-preM)   &
		      + KLR*rhohat*chat*(preP+preM-1.d0)  

    endsubroutine

    endmodule
    
    
    
    

    !================================
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
    