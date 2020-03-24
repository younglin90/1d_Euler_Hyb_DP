
module modContinuity

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getContinuity(solMethod)
        implicit none
        
        character(len=*) :: solMethod

        if( trim(solMethod) == 'explicit' ) then
            !CALL getRho_explicit
            stop
        elseif( trim(solMethod) == 'implicit' ) then
            CALL getRho_implicit
        endif

    endsubroutine
    

    !========================================================
    
    subroutine getRho_implicit
        use modFaceValues
        implicit none

        !do i=nFaceStr,nFaceEnd
        !    lCell = i; rCell = i+1
        !    wL(i,:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
        !    wR(i,:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
        !    CALL calfaceValAUSM(wL(i,:),wR(i,:),eosGamma, fMdot(i),fU(i),fP(i))
        !    fMdot(i) = fMdot(i)*area(i)
        !enddo
        CALL getFaceValues(ReconMethod=chReconR,&
            FluxMethod=chFluxR)
        
        !> Update to Density
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
                
            savecDummy = cRho(i)
            cRho(i) = oldcRho(i) + timestep/cVol(i)*( fMdot(lFace) - fMdot(rFace) )
            cDr(i) =  cRho(i) - savecDummy
        enddo
        
        
        !!> Variables initialization AUSM+
        !do i=nFaceStr,nFaceEnd
        !    lCell = i; rCell = i+1
        !    wL(i,:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
        !    wR(i,:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
        !    CALL calfaceValAUSM(wL(i,:),wR(i,:),eosGamma, fMdot(i),fU(i),fP(i))
        !    fMdot(i) = fMdot(i)*area(i)
        !    fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
        !    fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
        !enddo
        CALL getFaceValues(ReconMethod=chReconR,&
            FluxMethod=chFluxR,calfgP='on',calfgN='on')
        
        !> Update to Temperature
        AmatCon = 0.d0
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i 
            
            AmatCon(i,i) = 1.d0/timestep &
                - fgN(lFace)*fU(lFace)*area(lFace)/cVol(i) &
                + fgP(rFace)*fU(rFace)*area(rFace)/cVol(i)
            if( i/=nCellStr ) then
                AmatCon(i,i-1) = &
                    - fgP(lFace)*fU(lFace)*area(lFace)/cVol(i) 
            endif
            if( i/=nCellEnd ) then
                AmatCon(i,i+1) = &
                    + fgN(rFace)*fU(rFace)*area(rFace)/cVol(i)
            endif
            
            BmatCon(i) = &
                -( cRho(i) - oldCons(i,1) )/timestep &
                +1.d0/cVol(i)*( fMdot(lFace) - fMdot(rFace) )
        enddo
        !> inverse of A matrix
        CALL matInverse(AmatCon,nCellEnd)
        
        !> dR update
        cDR(nCellStr:nCellEnd) = &
            underRelaxFactR &
            *matmul(AmatCon(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                    BmatCon(nCellStr:nCellEnd))
        
        !> Rho update
        cRho(nCellStr:nCellEnd) = cRho(nCellStr:nCellEnd) + cDR(nCellStr:nCellEnd)
        
        

    endsubroutine

    !!========================================================
    !subroutine getRho_explicit
    !    implicit none
    !
    !    do i=nFaceStr,nFaceEnd
    !        lCell = i; rCell = i+1
    !        wL(i,:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
    !        wR(i,:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
    !        CALL calfaceValAUSM(wL(i,:),wR(i,:),eosGamma, fMdot(i),fU(i),fP(i))
    !        fMdot(i) = fMdot(i)*area(i)
    !    enddo
    !    !CALL getFaceValues(ReconMethod='FirUpwind',&
    !    !    FluxMethod='AUSM')
    !    
    !    !> Update to Density
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !            
    !        savecDummy = cRho(i)
    !        cRho(i) = oldcRho(i) + timestep/cVol(i)*( fMdot(lFace) - fMdot(rFace) )
    !        cDr(i) =  cRho(i) - savecDummy
    !    enddo
    !
    !endsubroutine

    
endmodule
    