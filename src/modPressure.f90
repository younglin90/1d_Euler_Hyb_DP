
module modPressure

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getPressure(solMethod)
        implicit none
        
        character(len=*) :: solMethod

        if( trim(solMethod) == 'explicit' ) then
         !   CALL getP_explicit
            stop
        elseif( trim(solMethod) == 'implicit' ) then
            CALL getP_implicit
        endif

    endsubroutine
    

    !========================================================
    subroutine getP_implicit
        use modFaceValues
        implicit none
        
        real(8) :: save_dummy(ncell)
        real(8) :: fL_cdp,fR_cdp,fL_drdp,fR_drdp,Rgas,lfs,rfs
        real(8) :: lPFlux,rPFlux,rRhoflux,lRhoflux
        
        CALL getFaceValues(&
            FaceMethodU=chFluxU &
            )
        
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
                
            dPN = dx
            
            !Rgas = eosCp*(1.d0-1.d0/eosGamma)
            
            !> @todo : 1/R/T => drdp
            lFs = fU(lFace)*( fgP(lFace)*cDRDP(i-1) &
                             +fgN(lFace)*cDRDP(i) )
            rFs = fU(rFace)*( fgP(rFace)*cDRDP(i) &
                             +fgN(rFace)*cDRDP(i+1) )
            
            AmatPres(i,i) = cDRDP(i)/timestep
            AmatPres(i,i) = AmatPres(i,i) &
                + 1.d0/cVol(i)*( rFs*(fgP(rFace))*area(rFace) &
                               - lFs*(fgN(lFace))*area(lFace) )&
                - 1.d0/cVol(i)*( (-1.d0)*area(rFace) &
                               - (+1.d0)*area(lFace) ) &
                              /dPN*timestep
            if( i /= nCellStr ) then
                AmatPres(i,i-1) =  &
                    + 1.d0/cVol(i)*( -lFs*(fgP(lFace))*area(lFace) )&
                    - 1.d0/cVol(i)*( -(-1.d0)*area(lFace) ) &
                                  /dPN*timestep
            endif
            if( i /= nCellEnd ) then
                AmatPres(i,i+1) =  &
                    + 1.d0/cVol(i)*( +rFs*(fgN(rFace))*area(rFace) )&
                    - 1.d0/cVol(i)*( +(+1.d0)*area(rFace) ) &
                                  /dPN*timestep
            endif
            
            BmatPres(i) = &
                - ( cRho(i) - oldCons(i,1) )/timestep &
                +1.d0/cVol(i)*( fMdot(lFace) - fMdot(rFace) )

        enddo
            
        !> Update Pressure
        CALL matInverse(AmatPres,nCellEnd)
        
        cDP(nCellStr:nCellEnd) = &
            underRelaxFactP &
            *matmul(AmatPres(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                    BmatPres(nCellStr:nCellEnd))
        cP(nCellStr:nCellEnd) = cP(nCellStr:nCellEnd) + cDP(nCellStr:nCellEnd)
        
    endsubroutine
    

    !!========================================================
    !subroutine getP_Xisto
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
    !        fUx(i) = fgP(i)*cU(lCell)  + fgN(i)*cU(rCell)
    !        fcoeffMom1(i) = fgP(i)*coeffMom1(lCell) + fgN(i)*coeffMom1(rCell)
    !        fcoeffMom2(i) = fgP(i)*coeffMom2(lCell) + fgN(i)*coeffMom2(rCell)
    !        fRho(i) = fgP(i)*cRho(lCell) + fgN(i)*cRho(rCell)
    !        fDrDp(i) = fgP(i)*cDrDp(lCell) + fgN(i)*cDrDp(rCell)
    !    enddo
    !    
    !    AmatPres = 0.d0
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !            
    !        dPN = dx
    !        
    !        Rgas = eosCp*(1.d0-1.d0/eosGamma)
    !        
    !        !> @todo : 1/R/T => drdp
    !        lFs = fU(lFace)*( fgP(lFace)*(1.d0/Rgas/oldcT(i-1)) &
    !                         +fgN(lFace)*(1.d0/Rgas/oldcT(i)) )
    !        rFs = fU(rFace)*( fgP(rFace)*(1.d0/Rgas/oldcT(i)) &
    !                         +fgN(rFace)*(1.d0/Rgas/oldcT(i+1)) )
    !        
    !        AmatPres(i,i) = (1.d0/Rgas/oldcT(i))/timestep
    !        AmatPres(i,i) = AmatPres(i,i) &
    !            + 1.d0/cVol(i)*( rFs*(fgP(rFace))*area(rFace) &
    !                           - lFs*(fgN(lFace))*area(lFace) )&
    !            - 1.d0/cVol(i)*( (-1.d0)*area(rFace) &
    !                           - (+1.d0)*area(lFace) ) &
    !                          /dPN*timestep
    !        if( i /= nCellStr ) then
    !        AmatPres(i,i-1) =  &
    !            + 1.d0/cVol(i)*( -lFs*(fgP(lFace))*area(lFace) )&
    !            - 1.d0/cVol(i)*( -(-1.d0)*area(lFace) ) &
    !                          /dPN*timestep
    !        endif
    !        if( i /= nCellEnd ) then
    !        AmatPres(i,i+1) =  &
    !            + 1.d0/cVol(i)*( +rFs*(fgN(rFace))*area(rFace) )&
    !            - 1.d0/cVol(i)*( +(+1.d0)*area(rFace) ) &
    !                          /dPN*timestep
    !        endif
    !        
    !        rPFlux = fgP(rFace)*oldcP(i) + fgN(rFace)*oldcP(i+1)
    !        lPFlux = fgP(lFace)*oldcP(i-1) + fgN(lFace)*oldcP(i)
    !        !BmatPres(i) = oldcRho(i)/timestep
    !        BmatPres(i) = 1.d0/cVol(i)*( lFs*lPFlux*area(lFace) &
    !                                    -rFs*rPFlux*area(rFace) ) &
    !                     -1.d0/cVol(i)*( (oldcP(i)-oldcP(i-1))*area(lFace) &
    !                                    -(oldcP(i+1)-oldcP(i))*area(rFace) )/dPN*timestep
    !
    !    enddo
    !        
    !    !> Update Pressure
    !    CALL matInverse(AmatPres,nCellEnd)
    !    save_dummy(nCellStr:nCellEnd) = cP(nCellStr:nCellEnd)
    !    !cP(nCellStr:nCellEnd) = &
    !    !    matmul(AmatPres(nCellStr:nCellEnd,nCellStr:nCellEnd),&
    !    !           BmatPres(nCellStr:nCellEnd))
    !    !cDp(nCellStr:nCellEnd) = 0.5 &
    !    !                        *cDp(nCellStr:nCellEnd)
    !    cP(nCellStr:nCellEnd) = oldcP(nCellStr:nCellEnd) &
    !                           +matmul(AmatPres(nCellStr:nCellEnd,nCellStr:nCellEnd),&
    !                                   BmatPres(nCellStr:nCellEnd))
    !    cDp(nCellStr:nCellEnd) = cP(nCellStr:nCellEnd) &
    !                           - save_dummy(nCellStr:nCellEnd)
    !    
    !endsubroutine
    !
    !
    !!========================================================
    !subroutine getP_Xiao
    !    implicit none
    !    
    !    
    !    do i=nFaceStr,nFaceEnd
    !        lCell = i; rCell = i+1
    !        wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
    !        wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
    !        CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
    !        fMdot(i) = fMdot(i)*area(i)
    !        fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
    !        fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
    !        fUx(i) = fgP(i)*cU(lCell)  + fgN(i)*cU(rCell)
    !        fcoeffMom1(i) = fgP(i)*coeffMom1(lCell) + fgN(i)*coeffMom1(rCell)
    !        fcoeffMom2(i) = fgP(i)*coeffMom2(lCell) + fgN(i)*coeffMom2(rCell)
    !        fRho(i) = fgP(i)*cRho(lCell) + fgN(i)*cRho(rCell)
    !        fDrDp(i) = fgP(i)*cDrDp(lCell) + fgN(i)*cDrDp(rCell)
    !    enddo
    !    
    !    AmatPres(:,:) = 0.d0
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !            
    !        dPN = dx
    !
    !        segAlp1 = area(lFace)*timestep &
    !                    /(0.5d0*(cRho(i-1)+cRho(i))) &
    !                    /dPN
    !        segAlp2 = area(rFace)*timestep &
    !                /(0.5d0*(cRho(i+1)+cRho(i))) &
    !                    /dPN
    !        AmatPres(i,i)   = AmatPres(i,i) -segAlp1
    !        AmatPres(i,i)   = AmatPres(i,i) -segAlp2
    !        if( i /= nCellStr ) then
    !            AmatPres(i,i-1) = segAlp1
    !        endif
    !        if( i /= nCellEnd ) then
    !            AmatPres(i,i+1) = segAlp2
    !        endif
    !        AmatPres(i,i) = AmatPres(i,i) - 1.d0/cRho(i)*cDrDp(i)/timestep
    !            
    !        BmatPres(i) = area(rFace) * ( fU(rFace) - timestep/fRho(rFace) &
    !                        *(cDp(i-1)-cDp(i))/dPN ) &
    !                    - area(lFace) * ( fU(lFace) - timestep/fRho(lFace) &
    !                        *(cDp(i+1)-cDp(i))/dPN )
    !    enddo
    !        
    !    !> Update Pressure
    !    CALL matInverse(AmatPres,nCellEnd)
    !    cDp(nCellStr:nCellEnd) = &
    !        matmul(AmatPres(nCellStr:nCellEnd,nCellStr:nCellEnd),&
    !                BmatPres(nCellStr:nCellEnd))
    !    !cDp(nCellStr:nCellEnd) = 0.5 &
    !    !                        *cDp(nCellStr:nCellEnd)
    !    cP(nCellStr:nCellEnd) = cP(nCellStr:nCellEnd) &
    !                            +cDp(nCellStr:nCellEnd)
    !    
    !endsubroutine
    !
    !
    !!========================================================
    !subroutine gePressure_YLL
    !    implicit none
    !
    !    do i=nFaceStr,nFaceEnd
    !        lCell = i; rCell = i+1
    !        wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
    !        wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
    !        CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
    !        fMdot(i) = fMdot(i)*area(i)
    !        fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
    !        fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
    !        fUx(i) = fgP(i)*cU(lCell)  + fgN(i)*cU(rCell)
    !        fcoeffMom1(i) = fgP(i)*coeffMom1(lCell) + fgN(i)*coeffMom1(rCell)
    !        fcoeffMom2(i) = fgP(i)*coeffMom2(lCell) + fgN(i)*coeffMom2(rCell)
    !        fRho(i) = fgP(i)*cRho(lCell) + fgN(i)*cRho(rCell)
    !        fDrDp(i) = fgP(i)*cDrDp(lCell) + fgN(i)*cDrDp(rCell)
    !    enddo
    !    AmatPres(:,:) = 0.d0
    !    do i=nCellStr,nCellEnd
    !        lFace = i-1; rFace = i
    !            
    !        dPN = dx
    !            
    !        
    !        segAlp1 = fDrDp(lFace)*fcoeffMom1(lFace)/cVol(i)*area(lFace)
    !        !segAlp2 = fcoeffMom2(lFace)/cVol(i)/timestep*area(lFace) 
    !        segAlp2 = timestep/cVol(i)*area(lFace) 
    !        
    !        
    !        AmatPres(i,i) = cDrDp(i)/timestep &
    !            - segAlp1*fgN(lFace) + segAlp1*fgP(rFace) &
    !            - segAlp2/dPN - segAlp2/dPN
    !        
    !        if( i/=nCellStr ) then
    !            AmatPres(i,i-1) = AmatPres(i,i-1) &
    !                - segAlp1*fgP(lFace) &
    !                + segAlp2/dPN
    !            
    !        endif
    !        if( i/=nCellEnd ) then
    !            AmatPres(i,i+1) = AmatPres(i,i+1) &
    !                + segAlp1*fgN(rFace) &
    !                + segAlp2/dPN
    !        endif
    !        
    !        
    !        BmatPres(i) = -segAlp2 &
    !            *( -( cDp(i) - cDp(i-1) )/dPN + ( cDp(i+1) - cDp(i) )/dPN ) &
    !            + fMdot(lFace) - fMdot(rFace)
    !            
    !    enddo
    !    
    !    !> Update Pressure
    !    CALL matInverse(AmatPres,nCellEnd)
    !    cDp(nCellStr:nCellEnd) = &
    !        matmul(AmatPres(nCellStr:nCellEnd,nCellStr:nCellEnd),&
    !                BmatPres(nCellStr:nCellEnd))
    !    !cDp(nCellStr:nCellEnd) = 0.001 &
    !    !                        *cDp(nCellStr:nCellEnd)
    !    cP(nCellStr:nCellEnd) = cP(nCellStr:nCellEnd) &
    !                            +cDp(nCellStr:nCellEnd)
    !            
    !    
    !
    !endsubroutine

    
endmodule
    