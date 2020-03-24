
module modEnergy

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getEnergy(solMethod)
        implicit none
        
        character(len=*) :: solMethod

        if( trim(solMethod) == 'explicit' ) then
            CALL getT_explicit
        elseif( trim(solMethod) == 'implicit' ) then
            CALL getT_implicit
        endif

    endsubroutine
    

    !========================================================
    subroutine getT_implicit
        use modFaceValues
        implicit none
        
        real(8) :: save_dummy(ncell)
        real(8) :: fL_cdp,fR_cdp,fL_drdp,fR_drdp,Rgas,lfs,rfs
        real(8) :: lfdt,rfdt,dhdT,dedT,lfHt,rfHt,calHt
        
        CALL getFaceValues(&
            FaceMethodU=chFluxU &
            )
        
        !> Update to Temperature
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i 
            
            dhdT = eosCp
            !dedT = eosCp/eosGamma
            
            AmatEnr(i,i) = cRho(i)*dhdT/timestep &
                - fgN(lFace)*dhdT*fMdot(lFace)/cVol(i) &
                + fgP(rFace)*dhdT*fMdot(rFace)/cVol(i)
            if( i/=nCellStr ) then
                AmatEnr(i,i-1) = &
                    - fgP(lFace)*dhdT*fMdot(lFace)/cVol(i) 
            endif
            if( i/=nCellEnd ) then
                AmatEnr(i,i+1) = &
                    + fgN(rFace)*dhdT*fMdot(rFace)/cVol(i)
            endif
            
            lfHt = fgP(lFace)*oldcHt(i-1)  + fgN(lFace)*oldcHt(i)
            rfHt = fgP(rFace)*oldcHt(i)  + fgN(rFace)*oldcHt(i+1)
            
            !> @todo : calc. Et 
            calHt = eosCp*cT(i) + 0.5d0*(cU(i)**2.d0)
            BmatEnr(i) = &
                -( cRho(i)*calHt - oldcRho(i)*oldcHt(i) )/timestep &
                +( cP(i) - oldcP(i) )/timestep &
                +1.d0/cVol(i)*( fMdot(lFace)*lfHt - fMdot(rFace)*rfHt )
        enddo
        !> inverse of A matrix
        CALL matInverse(AmatEnr,nCellEnd)
        
        !> dT update
        cDT(nCellStr:nCellEnd) = &
            underRelaxFactT &
            *matmul(AmatEnr(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                    BmatEnr(nCellStr:nCellEnd))
        
        !> T update
        cT(nCellStr:nCellEnd) = cT(nCellStr:nCellEnd) + cDT(nCellStr:nCellEnd)
        
        
    endsubroutine
    

    !========================================================
    subroutine getT_explicit
        use modFaceValues
        implicit none
        
        real(8) :: save_dummy(ncell)
        real(8) :: fL_cdp,fR_cdp,fL_drdp,fR_drdp,Rgas,lfs,rfs
        real(8) :: lfdt,rfdt,dhdT,dedT,lfHt,rfHt,calHt
        
        CALL getFaceValues(&
            FaceMethodU=chFluxU &
            )
        
        !> Update to Density
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
                
            dhdT = eosCp
            
            lfHt = fgP(lFace)*(cT(i-1)*eosCp+0.5d0*cU(i-1)**2.d0) &
                 + fgN(lFace)*(cT(i)*eosCp+0.5d0*cU(i)**2.d0)
            rfHt = fgP(rFace)*(cT(i)*eosCp+0.5d0*cU(i)**2.d0) &
                 + fgN(rFace)*(cT(i+1)*eosCp+0.5d0*cU(i+1)**2.d0)
            
            cDT(i) = &
                -(cRho(i)*cHt(i)-oldcRho(i)*oldcHt(i))/timestep&
                +(cP(i)-oldcP(i))/timestep &
                +1.d0/cVol(i)*( fMdot(lFace)*lfHt - fMdot(rFace)*rfHt )
            cDT(i) = underRelaxFactT*cDT(i)/cRho(i)/dhdT*timestep
            cT(i) =  cT(i) + cDT(i)
        enddo
        
    endsubroutine
    
    

    !!========================================================
    !subroutine getTemperature_explicit
    !    implicit none
    !    
    !    real(8) :: save_dummy(ncell)
    !    real(8) :: fL_cdp,fR_cdp,fL_drdp,fR_drdp,Rgas,lfs,rfs
    !    real(8) :: lPFlux,rPFlux
    !    
    !    !> Variables initialization AUSM+
    !    do i=nFaceStr,nFaceEnd
    !        lCell = i; rCell = i+1
    !        wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
    !        wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
    !        CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
    !        fMdot(i) = fMdot(i)*area(i)
    !        fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
    !        fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
    !        fHt(i) = fgP(i)*cHt(lCell)  + fgN(i)*cHt(rCell)
    !    enddo
    !    !> Update to Temperature
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !        
    !        savecDummy = cT(i)
    !        cT(i) = oldCons(i,3) + timestep/cVol(i) &
    !                *( fMdot(lFace)*fHt(lFace) - fMdot(rFace)*fHt(rFace) )
    !        cT(i) = ( cT(i)/cRho(i) - 0.5d0*(cU(i)**2.d0) )/(eosCp/eosGamma)
    !        cDT(i) = cT(i) - savecDummy
    !    enddo
    !    
    !    
    !endsubroutine
    !
    !
    !
    !!========================================================
    !subroutine getTemperature_implicit ! density constancy
    !    implicit none
    !    
    !    real(8) :: save_dummy(ncell)
    !    real(8) :: fL_cdp,fR_cdp,fL_drdp,fR_drdp,Rgas,lfs,rfs
    !    real(8) :: lfHt,rfHt,dhdT
    !    
    !    !> Variables initialization AUSM+
    !    do i=nFaceStr,nFaceEnd
    !        lCell = i; rCell = i+1
    !        wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
    !        wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
    !        CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
    !        fMdot(i) = fMdot(i)*area(i)
    !        fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
    !        fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
    !    enddo
    !    !> Update to Temperature
    !    AmatEnr = 0.d0
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !        
    !        lfHt = fgP(lFace)*oldcHt(i-1)  + fgN(lFace)*oldcHt(i)
    !        rfHt = fgP(rFace)*oldcHt(i)  + fgN(rFace)*oldcHt(i+1)
    !        
    !        dhdT = eosCp
    !        AmatEnr(i,i) = cRho(i)*dhdT/timestep &
    !            - fgN(lFace)*dhdT*fMdot(lFace)/cVol(i) &
    !            + fgP(rFace)*dhdT*fMdot(rFace)/cVol(i)
    !        if( i/=nCellStr ) then
    !            AmatEnr(i,i-1) = &
    !                - fgP(lFace)*dhdT*fMdot(lFace)/cVol(i) 
    !        endif
    !        if( i/=nCellEnd ) then
    !            AmatEnr(i,i+1) = &
    !                + fgN(rFace)*dhdT*fMdot(rFace)/cVol(i)
    !        endif
    !        
    !        BmatEnr(i) = &
    !            +( cP(i)-oldcP(i) + (oldcRho(i)-cRho(i))*oldcHt(i) )/timestep &
    !            +1.d0/cVol(i)*( fMdot(lFace)*lfHt - fMdot(rFace)*rfHt )
    !    enddo
    !    
    !    save_dummy(nCellStr:nCellEnd) = cT(nCellStr:nCellEnd)
    !    
    !    cT(nCellStr:nCellEnd) = cT(nCellStr:nCellEnd) &
    !                    +  matmul(AmatEnr(nCellStr:nCellEnd,nCellStr:nCellEnd),&
    !                              BmatEnr(nCellStr:nCellEnd))
    !    cDT(nCellStr:nCellEnd) = &
    !        cT(nCellStr:nCellEnd) - save_dummy(nCellStr:nCellEnd)
    !    
    !    
    !endsubroutine
    
endmodule
    