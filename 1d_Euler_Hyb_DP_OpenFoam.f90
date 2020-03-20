!======================
! 1D euler eq.
!=====================
program Euler

    implicit none

    real(8) :: CFL, time_max, eosGamma, eosCp, len, dx, &
     time, charVel, wL(3), wR(3), vecSf, dPN, avgDp, underRelaxFactP, &
     timestep, timestepMax, resiP, fVol, gP, gF, resiU, &
     maxP,minP,maxU,minU
    real(8), dimension(:), allocatable :: cX, cRho, cU, cP, cT, &
     cVol, cDt, area, fMdot, fU, fP, fT, cAcRU, cAfRU, cBcRU, &
     cAcUU, cAfUU, cAfUP, cBcUU, cDrDp, cHt, &
     fDrDp, fRho, fFirPE, fAcUU, fSecPE, fDp_fS, &
     cAcPDp, cAfPR, cAfPDp, cBcPDp, cDp, &
     fHt, cAcHH, cAfHU, cBcHH, fThrPE, fDpDn, fDp, cDU, &
     fUx
    integer :: nCell, nCellStr, nCellEnd, nTCellStr, nTCellEnd, &
     nFace, nFaceStr, nFaceEnd, i, nstep, lFace, rFace, iterSIMPLE, &
     iterMaxSIMPLE, lCell, rCell, iterPISO, iterMaxPISO
    
    
    
    CFL = 0.1
    nCell = 21
    time_max = 0.012 !sec
    eosGamma = 1.4 ! air
    eosCp = 1006.d0
    len = 1.0 ! m, x-length
    iterMaxSIMPLE = 10
    iterMaxPISO = 100
    underRelaxFactP = 0.1d0
    timestepMax = 1.d-2
    
    nCellStr = 1
    nCellEnd = nCell
    
    nTCellStr = 0
    nTCellEnd = nCell + 1
    allocate( cX(nTCellStr:nTCellEnd), &
              cRho(nTCellStr:nTCellEnd), &
              cU(nTCellStr:nTCellEnd), &
              cP(nTCellStr:nTCellEnd), &
              cT(nTCellStr:nTCellEnd), &
              cDt(nTCellStr:nTCellEnd), &
              cVol(nTCellStr:nTCellEnd), &
              cAcRU(nTCellStr:nTCellEnd), &
              cAfRU(nTCellStr:nTCellEnd), &
              cBcRU(nTCellStr:nTCellEnd), &
              cAcUU(nTCellStr:nTCellEnd), &
              cAfUU(nTCellStr:nTCellEnd), &
              cAfUP(nTCellStr:nTCellEnd), &
              cBcUU(nTCellStr:nTCellEnd), &
              cDrDp(nTCellStr:nTCellEnd), &
              cHt(nTCellStr:nTCellEnd), &
              cAcPDp(nTCellStr:nTCellEnd), &
              cAfPR(nTCellStr:nTCellEnd), &
              cAfPDp(nTCellStr:nTCellEnd), &
              cBcPDp(nTCellStr:nTCellEnd), &
              cDp(nTCellStr:nTCellEnd), &
              cAcHH(nTCellStr:nTCellEnd), &
              cAfHU(nTCellStr:nTCellEnd), &
              cBcHH(nTCellStr:nTCellEnd), &
              cDu(nTCellStr:nTCellEnd) )
    cDP(:) = 0.d0
    
    nFace = nCell+1
    nFaceStr = 0
    nFaceEnd = nCell
    allocate( fP(nFaceStr:nFaceEnd), &
              fU(nFaceStr:nFaceEnd), &
              fT(nFaceStr:nFaceEnd), &
              area(nFaceStr:nFaceEnd), &
              fMdot(nFaceStr:nFaceEnd), &
              fDrDp(nFaceStr:nFaceEnd), &
              fRho(nFaceStr:nFaceEnd), &
              fFirPE(nFaceStr:nFaceEnd), &
              fAcUU(nFaceStr:nFaceEnd), &
              fSecPE(nFaceStr:nFaceEnd), &
              fDp_fS(nFaceStr:nFaceEnd), &
              fHt(nFaceStr:nFaceEnd), &
              fDpDn(nFaceStr:nFaceEnd), &
              fThrPE(nFaceStr:nFaceEnd), &
              fDp(nFaceStr:nFaceEnd), &
              fUx(nFaceStr:nFaceEnd) )
    
    
    !-------- one-dimensional Mesh
    dx = len/real(ncell,8)
    do i=nTCellStr,nTCellEnd
      ! x value of cell center point
      cX(i) = 0.5d0*dx + (real(i,8)-1.d0)*dx
      cVol(i) = dx*1.d0*1.d0
    enddo
    do i=nFaceStr,nFaceEnd
      area(i) = 1.d0
    enddo
    
    !-------- initial condition
    ! w=(rho,u,p), primitive variables
    do i=nTCellStr,nTCellEnd
      !> sod shock condition
      if( cX(i) < 0.8d0 ) then
          cRho(i) = 1.000d0
          cU(i)   = 0.000d0
          cP(i)   = 1.000d0
      else
          cRho(i) = 0.125d0
          cU(i)   = 0.000d0
          cP(i)   = 0.100d0
      endif
    enddo
    
    time = 0.d0
    nstep=0
 
    do while( time <= time_max )

        nstep=nstep+1

        !---- calculation time-step
        do i=nCellStr,nCellEnd
            lFace = i-1
            rFace = i
            charVel = dabs(cU(i))*dmax1(area(lFace),area(rFace))
            cDt(i) = dmin1( timestepMax, CFL*cVol(i)/charVel ) ! sec
        enddo 
        timestep = minval(cDt(nCellStr:nCellEnd))
        cDt(:) = timestep
        
        SIMPLE: do iterSIMPLE = 1,iterMaxSIMPLE
        
            !!> Variables initialization AUSM+
            !!> mass flow rate, velocity, pressure of cell interface
            !do i=nFaceStr,nFaceEnd
            !    lCell = i
            !    rCell = i+1
            !    wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
            !    wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
            !    CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
            !    fMdot(i) = fMdot(i)*area(i)
            !    gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            !    gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
            !    fUx(i) = gP*cU(lCell)  + gF*cU(rCell)
            !enddo
            !
            !!=======================================================
            !!> Density prediction (explicit)
            !!> Continuity equation
            !do i=nCellStr,nCellEnd
            !    lFace = i-1
            !    rFace = i
            !    
            !    cAcRU(i) = 1.d0 *cVol(i)/cDt(i)
            !    cAfRU(i) = fMdot(rFace) &
            !             - fMdot(lFace)
            !    cBcRU(i) = cRho(i) *cVol(i)/cDt(i)
            !    
            !    cRho(i) = ( cBcRU(i) - cAfRU(i)) / cAcRU(i)
            !enddo
            !!> B.C.
            !cRho(0) = cRho(1)
            !cRho(ncell+1) = cRho(ncell)
            !!=======================================================
            !
            !!=======================================================
            !!> Velocity prediction (implicit)
            !!> Momentum equation
            !do i=nCellStr,nCellEnd
            !    lFace = i-1
            !    rFace = i
            !    
            !    cAcUU(i) = cRho(i) *cVol(i)/cDt(i)
            !    cAfUU(i) = fMdot(rFace)*fUx(rFace) &
            !             - fMdot(lFace)*fUx(lFace)
            !    cAfUP(i) = fP(rFace)*area(rFace) &
            !             - fP(lFace)*area(lFace)
            !    cBcUU(i) = cRho(i)*cU(i) *cVol(i)/cDt(i)
            !    
            !    cU(i) = ( cBcUU(i) - cAfUU(i) - cAfUP(i) ) / cAcUU(i)
            !    
            !enddo
            !!> B.C.
            !cU(0) = cU(1)
            !cU(ncell+1) = cU(ncell)
            !cAcUU(0) = cAcUU(1)
            !cAcUU(ncell+1) = cAcUU(ncell)
            !!=======================================================
            
        
            !> PISO iteration start
            PISO: do iterPISO = 1,iterMaxPISO
                
                
            !> Variables initialization AUSM+
            !> mass flow rate, velocity, pressure of cell interface
            do i=nFaceStr,nFaceEnd
                lCell = i
                rCell = i+1
                wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
                wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
                CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
                fMdot(i) = fMdot(i)*area(i)
                gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
                gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
                fUx(i) = gP*cU(lCell)  + gF*cU(rCell)
            enddo
        
            !=======================================================
            !> Density prediction (explicit)
            !> Continuity equation
            do i=nCellStr,nCellEnd
                lFace = i-1
                rFace = i
                
                cAcRU(i) = 1.d0 *cVol(i)/cDt(i)
                cAfRU(i) = fMdot(rFace) &
                         - fMdot(lFace)
                cBcRU(i) = cRho(i) *cVol(i)/cDt(i)
                
                cRho(i) = ( cBcRU(i) - cAfRU(i)) / cAcRU(i)
            enddo
            !> B.C.
            cRho(0) = cRho(1)
            cRho(ncell+1) = cRho(ncell)
            !=======================================================
        
            !=======================================================
            !> Velocity prediction (implicit)
            !> Momentum equation
            do i=nCellStr,nCellEnd
                lFace = i-1
                rFace = i
                
                cAcUU(i) = cRho(i) *cVol(i)/cDt(i)
                cAfUU(i) = fMdot(rFace)*fUx(rFace) &
                         - fMdot(lFace)*fUx(lFace)
                cAfUP(i) = fP(rFace)*area(rFace) &
                         - fP(lFace)*area(lFace)
                cBcUU(i) = cRho(i)*cU(i) *cVol(i)/cDt(i)
                
                cU(i) = ( cBcUU(i) - cAfUU(i) - cAfUP(i) ) / cAcUU(i)
                
            enddo
            !> B.C.
            cU(0) = cU(1)
            cU(ncell+1) = cU(ncell)
            cAcUU(0) = cAcUU(1)
            cAcUU(ncell+1) = cAcUU(ncell)
            !=======================================================
                
                
                
                
                !!> Velocity prediction (implicit)
                !!> Momentum equation, neglecting the gradient of pressure
                !if( iterPISO /= 1 ) then
                !    !> Variables initialization AUSM+
                !    !> mass flow rate, velocity, pressure of cell interface
                !    do i=nFaceStr,nFaceEnd
                !        lCell = i
                !        rCell = i+1
                !        wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
                !        wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
                !        CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
                !    enddo
                !    
                !    do i=nCellStr,nCellEnd
                !        lFace = i-1
                !        rFace = i
                !
                !        cAcUU(i) = cRho(i) *cVol(i)/cDt(i)
                !        cAfUU(i) = fMdot(rFace)*fU(rFace)*area(rFace) &
                !                 - fMdot(lFace)*fU(lFace)*area(lFace)
                !        cBcUU(i) = cRho(i)*cU(i) *cVol(i)/cDt(i)
                !        
                !        cU(i) =  ( cBcUU(i) - cAfUU(i) ) / cAcUU(i)
                !    enddo
                !    !> B.C.
                !    cU(0) = cU(1)
                !    cU(ncell+1) = cU(ncell)
                !    cAcUU(0) = cAcUU(1)
                !    cAcUU(ncell+1) = cAcUU(ncell)
                !endif
                !
                
                
                !=======================================================
                !> Variables initialization AUSM+
                !> mass flow rate, velocity, pressure of cell interface
                do i=nFaceStr,nFaceEnd
                    lCell = i
                    rCell = i+1
                    wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
                    wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
                    CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
                    fMdot(i) = fMdot(i)*area(i)
                enddo
                !=======================================================
                
                
                !=======================================================
                !> Calc. Compressible Coefficient from EOS
                !> Calc. Total Enthalpy from EOS
                TCellEOS : do i=nTCellStr,nTCellEnd
                    cT(i) = eosGamma/(eosGamma-1.d0)/eosCp/cRho(i)*cP(i)
                    cDrDp(i) = eosGamma/(eosGamma-1.d0)/eosCp/cT(i)
                    cHt(i) = eosCp*cT(i) + 0.5d0*(cU(i)**2.d0)
                enddo TCellEOS
                !=======================================================
                
        
                !=======================================================
                !> Pressure Correction (implicit)
                !> Pressure equation : Continuity + Momentum equations
                FacePress : do i=nFaceStr,nFaceEnd
                    lCell = i
                    rCell = i+1
                    
                    gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
                    gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
                    fDrDp(i) = gP*cDrDp(lCell)  + gF*cDrDp(rCell)
                    fRho(i)  = gP*cRho(lCell)   + gF*cRho(rCell)
                    fAcUU(i) = gP*cAcUU(lCell)  + gF*cAcUU(rCell)
                    fDp(i)   = gP*cDp(lCell)    + gF*cDp(rCell)
                    fUx(i)   = gP*cU(lCell)  + gF*cU(rCell)
                    
                    fVol = 0.5d0*( cVol(lCell) + cVol(rCell) )
                    
                    fDpDn(i) = ( cDp(rCell) - cDp(lCell) ) / ( cX(rCell) - cX(lCell) )

                    
                    fFirPE(i) = fRho(i)*fUx(i)*area(i)
                    fSecPE(i) = fDrDp(i)*fDp(i)*fUx(i)*area(i) &
                              - fRho(i)*fVol*fDpDn(i)/fAcUU(i)*area(i)
                enddo FacePress
                CellPress: do i=nCellStr,nCellEnd
                    lFace = i-1
                    rFace = i
                    
                    !cAcPDp(i) = cDrDp(i) *cVol(i)/cDt(i)
                    !cAfPR(i)  = fFirPE(rFace)*area(rFace) - fFirPE(lFace)*area(lFace)
                    !cAfPDp(i) = fSecPE(rFace)*fDp_fS(rFace) - fSecPE(lFace)*fDp_fS(lFace)
                    !cBcPDp(i) = 0.d0
                    !
                    !cDp(i) = ( cBcPDp(i) - cAfPR(i) - cAfPDp(i) ) / cAcPDp(i)
                    !
                    !cDp(i) = underRelaxFactP*cDp(i)
                    !cP(i) = cP(i) + cDp(i)
                    
                    cAcPDp(i) = cDrDp(i) *cVol(i)/cDt(i)
                    cAfPR(i)  = fFirPE(rFace) - fFirPE(lFace)
                    cAfPDp(i) = fSecPE(rFace) - fSecPE(lFace)
                    cBcPDp(i) = 0.d0
                    
                    cDp(i) = ( cBcPDp(i) - cAfPR(i) - cAfPDp(i) ) / cAcPDp(i)
                
                    cDp(i) = underRelaxFactP*cDp(i)
                    cP(i) = cP(i) + cDp(i)
                    
                enddo CellPress
                !> B.C.
                cP(0) = cP(1)
                cP(ncell+1) = cP(ncell)
                !=======================================================
                
                
                !=======================================================
                !> Velocity Correction (explicit)
                !> Momentum equation, neglecting the pressure term
                do i=nFaceStr,nFaceEnd
                    lCell = i
                    rCell = i+1
                    gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
                    gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
                    fDp(i) = gP*cDp(lCell) + gF*cDp(rCell)
                enddo
                do i=nCellStr,nCellEnd
                    lFace = i-1
                    rFace = i
                    cDU(i) = -cVol(i)/cAcUU(i) &
                        * ( fDp(lFace)*area(lFace) &
                           -fDp(rFace)*area(rFace) )/cVol(i)
                    cU(i) = cU(i) + cDU(i)
                enddo
                !> B.C.
                cU(0) = cU(1)
                cU(ncell+1) = cU(ncell)
                !=======================================================
                
                !=======================================================
                !> Mass flux Correction (explicit)
                do i=nFaceStr,nFaceEnd
                    lCell = i
                    rCell = i+1
                    gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
                    gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
                    fDrDp(i) = gP*cDrDp(lCell)  + gF*cDrDp(rCell)
                    fRho(i)  = gP*cRho(lCell)   + gF*cRho(rCell)
                    fAcUU(i) = gP*cAcUU(lCell)  + gF*cAcUU(rCell)
                    fDp(i)   = gP*cDp(lCell)    + gF*cDp(rCell)
                    fUx(i)   = gP*cU(lCell)  + gF*cU(rCell)
                    
                    fSecPE(i) = fDrDp(i)*fDp(i)*fUx(i)*area(i) &
                              - fRho(i)*fVol*fDpDn(i)/fAcUU(i)*area(i)
                    
                    fMdot(i) = fMdot(i) - fSecPE(i)
                enddo
                !=======================================================
                

                !> Pressure Residual
                resiP = 0.d0
                maxP = -1.d8
                minP = 1.d8
                do i=nCellStr,nCellEnd
                    resiP = resiP + cDp(i)**2.d0
                    maxP = dmax1(maxP,cP(i))
                    minP = dmin1(minP,cP(i))
                enddo
                resiP = dsqrt(resiP) / (dmax1(maxP,0.d0)-dmin1(minP,0.d0)+1.d-15)
                
                !> Velocity Residual
                resiU = 0.d0
                maxU = -1.d8
                minU = 1.d8
                do i=nCellStr,nCellEnd
                    resiU = resiU + cDU(i)**2.d0
                    maxU = dmax1(maxU,cU(i))
                    minU = dmin1(minU,cU(i))
                enddo
                resiU = dsqrt(resiU) / (dmax1(maxU,0.d0)-dmin1(minU,0.d0)+1.d-15)
                
                
                write(*,100)iterSIMPLE,iterPISO,resiP,resiU
                if(iterPISO==99) stop
                

            enddo PISO !> End iteration PISO
            
            stop
            
            
            !> Temperature Correction (explicit)
            !> Energy equation
            do i=nFaceStr,nFaceEnd
                lCell = i
                rCell = i+1
                    
                fHt(i) = 0.5d0*( (1.d0+dsign(1.d0,fU(i)))*cHt(lCell) &
                                +(1.d0-dsign(1.d0,fU(i)))*cHt(rCell) )
            enddo
            do i=nCellStr,nCellEnd
                lFace = i-1
                rFace = i
                    
                !> Calc. T to Ht from the EOS
                
                cAcHH(i) = cRho(i) *cVol(i)/cDt(i)
                cAfHU(i) = fMdot(rFace)*fHt(rFace)*area(rFace) &
                            - fMdot(lFace)*fHt(lFace)*area(lFace)
                cBcHH(i) = cRho(i)*cHt(i) *cVol(i)/cDt(i) &
                            + cDp(i) *cVol(i)/cDt(i) 
                    
                cHt(i) = ( cBcHH(i) - cAfHU(i) ) / cAcHH(i)
                    
                !> Calc. Ht to T from the EOS
                cT(i) = ( cHt(i) - 0.5d0*(cU(i)**2.d0) ) / eosCp
            enddo
            !> B.C.
            cT(0) = cT(1)
            cT(ncell+1) = cT(ncell)
            

        enddo SIMPLE !> End iteration SIMPLE
        
        
        time = time + timestep
        
    !/////////////////////////////////
        print*,nstep,time
    !/////////////////////////////////
    
    enddo
    
    
    
    !/////////////////////////////////
    open(20,file='a.dat')!,status='old',position='append')
    do i=nCellStr,nCellEnd
        WRITE( 20, '(20(F20.9, 1X))' ) cX(i),cRho(i),cU(i),cP(i)
    enddo
    close(20)
    !/////////////////////////////////
    
100 format(i7,1x,i7,1x,6(1pe10.2))
 
end program Euler


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
