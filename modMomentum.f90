
module modMomentum

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getMomentum(solMethod)
        implicit none
        
        character(len=*) :: solMethod

        if( trim(solMethod) == 'impCoeff' ) then
            CALL getMomentumCoeff
        elseif( trim(solMethod) == 'dp_to_u' ) then
            CALL getVelocity_explicit
        endif
        

    endsubroutine
    
    subroutine getMomentumCoeff
        implicit none
        real(8) :: save_dummy(ncell),save_Bmat(ncell)
    
        CALL getCoeff
        
        save_dummy(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd)

        !> Update to Velocity
        coeffMom1(nCellStr:nCellEnd) = &
                        matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                               BmatMom1(nCellStr:nCellEnd))
        
        coeffMom2(nCellStr:nCellEnd) = &
                        matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                               BmatMom2(nCellStr:nCellEnd))
        
        cDU(nCellStr:nCellEnd) = ( coeffMom1(nCellStr:nCellEnd) &
                                 +coeffMom2(nCellStr:nCellEnd) ) &
                               - save_dummy(nCellStr:nCellEnd)
        
        cU(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd) &
            + underRelaxFactU*cDU(nCellStr:nCellEnd)
    
    endsubroutine
    
    
    subroutine getCoeff
        implicit none
    
            !> Variables initialization AUSM+
            do i=nFaceStr,nFaceEnd
                lCell = i; rCell = i+1
                wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
                wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
                CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
                fMdot(i) = fMdot(i)*area(i)
                fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
                fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
                fUx(i) = fgP(i)*cU(lCell) + fgN(i)*cU(rCell)
            enddo
            !> Calc. Matix of Linear System about Momentum Eq.
            !> A*x=B
            AmatMom(:,:) = 0.d0
            do i=nCellStr,nCellEnd
                lFace = i-1; rFace = i
                    !> | (i-1) |-> (i) |-> (i+1) |
                AmatMom(i,i) = &
                                + cRho(i)/timestep &
                                - fgN(lFace)*fMdot(lFace) / cVol(i) &
                                + fgP(rFace)*fMdot(rFace) / cVol(i)
                if( i/=nCellStr ) then
                    AmatMom(i,i-1) = &
                                - fgP(lFace)*fMdot(lFace) / cVol(i)
                endif
                if( i/=nCellEnd ) then
                    AmatMom(i,i+1) = &
                                + fgN(rFace)*fMdot(rFace) / cVol(i)
                endif
                
                BmatMom1(i)   = oldCons(i,2) / timestep
                BmatMom2(i)   = ( fP(lFace)*area(lFace) &
                                - fP(rFace)*area(rFace) ) / cVol(i)
            enddo
            !> inverse of AmatMom
            CALL matInverse(AmatMom,nCellEnd)
    
    endsubroutine
    

    !========================================================
    
    subroutine getVelocity_explicit
        implicit none
        
        real(8) :: save_dummy(ncell),save_Bmat(ncell)
        real(8) :: lfVelx,rfVelx,lfPx,rfPx

        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
            wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
            CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
            fMdot(i) = fMdot(i)*area(i)
            fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
        enddo
        !> Calc. Matix of Linear System about Momentum Eq.
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
            
            save_dummy(i) = cU(i)
            
            lfPx = fP(lFace)
            rfPx = fP(rFace)
            
            lfPx = fgP(lFace)*cDp(i-1)  + fgN(lFace)*cDp(i)
            rfPx = fgP(rFace)*cDp(i)  + fgN(rFace)*cDp(i+1)
            
            save_Bmat(i) = ( lfPx*area(lFace) - rfPx*area(rFace) )
            
            
        enddo
        
            cDU(nCellStr:nCellEnd) = &
                underRelaxFactU &
                *matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                        save_Bmat(nCellStr:nCellEnd))
        
            cU(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd) + cDU(nCellStr:nCellEnd)
    
    endsubroutine
    
    
    
    

    subroutine getVelocity3
        implicit none

        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
            wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
            CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
            fMdot(i) = fMdot(i)*area(i)
            fgP(i) = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            fgN(i) = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
            fUx(i) = fgP(i)*cU(lCell) + fgN(i)*cU(rCell)
        enddo
        !> Calc. Matix of Linear System about Momentum Eq.
        !> A*x=B
        AmatMom(:,:) = 0.d0
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
                !> | (i-1) |-> (i) |-> (i+1) |
            AmatMom(i,i) = AmatMom(i,i) &
                            + cRho(i)/timestep &
                            + fMdot(rFace)*fgP(rFace) / cVol(i) &
                            - fMdot(lFace)*fgN(lFace) / cVol(i)
            if( i/=nCellStr ) then
                AmatMom(i,i-1) = AmatMom(i,i-1) &
                                - fMdot(lFace)*fgP(lFace) / cVol(i)
            endif
            if( i/=nCellEnd ) then
                AmatMom(i,i+1) = AmatMom(i,i+1) &
                                + fMdot(rFace)*fgN(rFace) / cVol(i)
            endif
            BmatMom1(i)   = oldcRho(i)*oldcU(i)/timestep
            BmatMom2(i)   = ( fP(lFace)*area(lFace) &
                            - fP(rFace)*area(rFace) ) / cVol(i)
        enddo
        !> inverse of AmatMom
        CALL matInverse(AmatMom,nCellEnd)
        !> Update to Velocity
        coeffMom1(nCellStr:nCellEnd) = &
                        matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                                BmatMom1(nCellStr:nCellEnd))
        coeffMom2(nCellStr:nCellEnd) = &
                        matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                                BmatMom2(nCellStr:nCellEnd))
        cDu(nCellStr:nCellEnd) = ( coeffMom1(nCellStr:nCellEnd) &
                                    -coeffMom2(nCellStr:nCellEnd) ) &
                                    - cU(nCellStr:nCellEnd) 
        cU(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd) &
                                + cDu(nCellStr:nCellEnd)
    
    
    endsubroutine
    
    
    
    
    subroutine getVelocity2
        implicit none
        
        real(8) :: save_dummy(ncell)

        do i=nFaceStr,nFaceEnd
            lCell = i; rCell = i+1
            wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
            wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
            CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
        enddo
        !> Calc. Matix of Linear System about Momentum Eq.
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
            !BmatMom1(i)   = oldcRho(i)*oldcU(i)/timestep
            BmatMom2(i)   = ( fP(lFace)*area(lFace) &
                            - fP(rFace)*area(rFace) ) / cVol(i)
        enddo
        !> inverse of AmatMom
        !CALL matInverse(AmatMom,nCellEnd)
        !> Update to Velocity
        !coeffMom1(nCellStr:nCellEnd) = &
        !                matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
        !                        BmatMom1(nCellStr:nCellEnd))
        
        save_dummy(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd)
        
        coeffMom2(nCellStr:nCellEnd) = &
                        matmul(AmatMom(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                               BmatMom2(nCellStr:nCellEnd))
        
        cU(nCellStr:nCellEnd) = ( coeffMom1(nCellStr:nCellEnd) &
                                 +coeffMom2(nCellStr:nCellEnd) ) 
        
        cDu(nCellStr:nCellEnd) = cU(nCellStr:nCellEnd) &
                                -save_dummy(nCellStr:nCellEnd)
    
    
    endsubroutine
    
    
    
endmodule
    