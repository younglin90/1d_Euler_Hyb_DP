
module modContinuity

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getContinuity(solMethod)
        implicit none
        
        character(len=*) :: solMethod

        if( trim(solMethod) == 'explicit' ) then
            CALL getRho_explicit
            !stop
        elseif( trim(solMethod) == 'implicit' ) then
            CALL getRho_implicit
        endif

    endsubroutine
    

    !========================================================
    
    subroutine getRho_implicit
        use modFaceValues
        implicit none

        CALL getFaceValues(&
            FaceMethodU=chFluxU &
            )
        
        !> Update to Density
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

    !========================================================
    subroutine getRho_explicit
        use modFaceValues
        implicit none
    
        CALL getFaceValues(&
            FaceMethodU=chFluxU &
            )
        
        !> Update to Density
        do i=nCellStr,nCellEnd
            lFace = i-1; rFace = i
                
            cDR(i) = &
                - (cRho(i)-oldcRho(i)) &
                + timestep/cVol(i)*( fMdot(lFace) - fMdot(rFace) )
            cRho(i) = cRho(i) + underRelaxFactR*cDR(i)
        enddo
    
    endsubroutine

    
endmodule
    