!======================
! 1D euler eq.
!=====================
program Euler

    implicit none

    real(8) :: CFL, time_max, eosGamma, eosCp, len, dx, &
     time, charVel, wL(3), wR(3), vecSf, dPN, avgDp, underRelaxFactP, &
     timestep, timestepMax, resiP, fVol, gP, gF, resiU, &
     maxP,minP,maxU,minU,segAlp1,segAlp2,gN,resiR,maxR,minR
    real(8), dimension(:), allocatable :: cX, cRho, cU, cP, cT, &
     cVol, cDt, area, fMdot, fU, fP, fT, cAcRU, cAfRU, cBcRU, &
     cAcUU, cAfUU, cAfUP, cBcUU, cDrDp, cHt, &
     fDrDp, fRho, fFirPE, fAcUU, fSecPE, fDp_fS, &
     cAcPDp, cAfPR, cAfPDp, cBcPDp, cDp, &
     fHt, cAcHH, cAfHU, cBcHH, fThrPE, fDpDn, fDp, cDU, &
     fUx, cEt, fEt, cC, cGECon, cGEMom, cGEEnr, &
     fGEConFlux, fGEMomFlux, fGEEnrFlux, segBmat,cDr
    real(8), dimension(:,:), allocatable :: segAmat
    integer :: nCell, nCellStr, nCellEnd, nTCellStr, nTCellEnd, &
     nFace, nFaceStr, nFaceEnd, i, nstep, lFace, rFace, iterSIMPLE, &
     iterMaxSIMPLE, lCell, rCell, iterPISO, iterMaxPISO
    
    
    
    CFL = 0.1
    nCell = 100
    time_max = 0.012 !sec
    eosGamma = 1.4 ! air
    eosCp = 1006.d0
    len = 1.0 ! m, x-length
    iterMaxSIMPLE = 10000
    underRelaxFactP = 0.1d0
    !timestepMax = 1.d-2
    timestep = 2.d-3
    
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
              cDu(nTCellStr:nTCellEnd), &
              cEt(nTCellStr:nTCellEnd), &
              cC(nTCellStr:nTCellEnd), &
              cGECon(nTCellStr:nTCellEnd), &
              cGEMom(nTCellStr:nTCellEnd), &
              cGEEnr(nTCellStr:nTCellEnd), &
              cDr(nTCellStr:nTCellEnd) )
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
              fUx(nFaceStr:nFaceEnd), &
              fEt(nFaceStr:nFaceEnd), &
              fGEConFlux(nFaceStr:nFaceEnd), &
              fGEMomFlux(nFaceStr:nFaceEnd), &
              fGEEnrFlux(nFaceStr:nFaceEnd) )
    
    allocate( segAmat(nCellEnd,nCellEnd), &
              segBmat(nCellEnd) )
    
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
    !--------- Ht, Et, T, c from EOS
    do i=nTCellStr,nTCellEnd
        cT(i) = cP(i)/cRho(i)/eosCp/(1.d0-1.d0/eosGamma)
        cHt(i) = eosCp*cT(i) + 0.5d0*(cU(i)**2.d0)
        cEt(i) = cHt(i) - cP(i)/cRho(i)
        cC(i) = dsqrt(eosGamma*cP(i)/cRho(i))
        cDrDp(i) = 1.d0/eosCp/(1.d0-1.d0/eosGamma)/cT(i)
    enddo
    
    time = 0.d0
    nstep=0
 
    do while( time <= time_max )
        
        nstep=nstep+1

        !!---- calculation time-step
        !do i=nCellStr,nCellEnd
        !    lFace = i-1
        !    rFace = i
        !    charVel = (dabs(cU(i))+cC(i))*dmax1(area(lFace),area(rFace))
        !    cDt(i) = cVol(i)/charVel ! sec
        !enddo 
        !timestep = minval(cDt(nCellStr:nCellEnd))
        !print*,timestep
        !stop
        
        cDt(:) = timestep
        
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        !------------ Density based solver ----
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
            gP = 0.5d0*(1.d0+dsign(1.d0,fU(i)))
            gF = 0.5d0*(1.d0-dsign(1.d0,fU(i)))
            fUx(i) = gP*cU(lCell)  + gF*cU(rCell)
            
            if( fMdot(i) >= 0.d0 ) then
                fGEConFlux(i) = fMdot(i)
                fGEMomFlux(i) = fMdot(i)*cU(lCell) !+ fP(i)*area(i)
                fGEEnrFlux(i) = fMdot(i)*cHt(lCell)
            else
                fGEConFlux(i) = fMdot(i)
                fGEMomFlux(i) = fMdot(i)*cU(rCell) !+ fP(i)*area(i)
                fGEEnrFlux(i) = fMdot(i)*cHt(rCell)
            endif
        enddo
        !=======================================================
        !> Update to Gorverning Equations
        do i=nCellStr,nCellEnd
            lFace = i-1
            rFace = i
                
            cGECon(i) = cRho(i) &
                      + cDt(i)/cVol(i)*( fGEConFlux(rFace) - fGEConFlux(lFace) )
            cGEMom(i) = cRho(i)*cU(i) &
                      + cDt(i)/cVol(i)*( fGEMomFlux(rFace) - fGEMomFlux(lFace) )
            !cGEEnr(i) = cRho(i)*cEt(i) &
            !          + cDt(i)/cVol(i)*( fGEEnrFlux(rFace) - fGEEnrFlux(lFace) )
            
        enddo
        !=======================================================
        !> Pressure, Density
        do i=nCellStr,nCellEnd
            cRho(i) = cGECon(i)
            cU(i) = cGEMom(i)/cRho(i)
        enddo
        
        
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        !----------------- Segregated Method ----------------
        Segregated: do iterSIMPLE = 1,iterMaxSIMPLE
        
            
            !=======================================================
            do i=nFaceStr,nFaceEnd
                lCell = i
                rCell = i+1
                
                wL(:) = (/cRho(lCell),cU(lCell),cP(lCell)/)
                wR(:) = (/cRho(rCell),cU(rCell),cP(rCell)/)
                CALL calfaceValAUSM(wL,wR,eosGamma, fMdot(i),fU(i),fP(i))
                fMdot(i) = fMdot(i)*area(i)
                
                if( fMdot(i) >= 0.d0 ) then
                    fRho(i) = cRho(lCell)
                    fDp(i) = cDp(lCell)
                    fDrDp(i) = cDrDp(lCell)
                else
                    fRho(i) = cRho(rCell)
                    fDp(i) = cDp(rCell)
                    fDrDp(i) = cDrDp(rCell)
                endif
            enddo
            !=======================================================
            segAmat(:,:) = 0.d0
            do i=nCellStr,nCellEnd
                lFace = i-1
                rFace = i
                
                dPN = dx
                
                segAlp1 = area(lFace)*cDt(i-1)/dPN/cVol(i)
                segAmat(i,i)   = segAmat(i,i) -segAlp1
                segAlp2 = area(lFace)*fDrDp(lFace)*fU(lFace)/cVol(i)
                gP = 0.5d0*(1.d0+dsign(1.d0,fU(lFace)))
                gN = 0.5d0*(1.d0-dsign(1.d0,fU(lFace)))
                segAmat(i,i) = segAmat(i,i) + segAlp2*gP*cDp(i)
                if( i-1 >= nCellStr ) then
                    segAmat(i,i-1) = segAmat(i,i-1) + segAlp1
                    segAmat(i,i-1) = segAmat(i,i-1) + segAlp2*gN*cDp(i-1)
                endif
                
                segAlp1 = area(rFace)*cDt(i+1)/dPN/cVol(i)
                segAmat(i,i)   = segAmat(i,i) -segAlp1
                segAlp2 = area(rFace)*fDrDp(rFace)*fU(rFace)/cVol(i)
                gP = 0.5d0*(1.d0+dsign(1.d0,fU(rFace)))
                gN = 0.5d0*(1.d0-dsign(1.d0,fU(rFace)))
                segAmat(i,i) = segAmat(i,i) + segAlp2*gP*cDp(i)
                if( i+1 <= nCellEnd ) then
                    segAmat(i,i+1) = segAmat(i,i+1) + segAlp1
                    segAmat(i,i+1) = segAmat(i,i+1) + segAlp2*gN*cDp(i+1)
                endif
                segAmat(i,i) = segAmat(i,i) + cDrDp(i)/cDt(i)
                
                segBmat(i) = area(lFace) * ( - cDt(i)*( (cDp(i-1)-cDp(i))/dPN ) ) &
                           - area(rFace) * ( - cDt(i)*( (cDp(i+1)-cDp(i))/dPN ) ) &
                           + fMdot(lFace) - fMdot(rFace)
                segBmat(i) = segBmat(i) &
                           + area(rFace)*cDt(i)*(cP(i+1)-cP(i))/dPN &
                           - area(lFace)*cDt(i)*(cP(i)-cP(i-1))/dPN
                segBmat(i) = segBmat(i) / cVol(i)
                
            enddo
            !=======================================================
            
            
            CALL matInverse(segAmat,nCellEnd)
            
            cDp(nCellStr:nCellEnd) = &
                matmul(segAmat(nCellStr:nCellEnd,nCellStr:nCellEnd),&
                       segBmat(nCellStr:nCellEnd))
            cDp(nCellStr:nCellEnd) = 0.1d0*cDp(nCellStr:nCellEnd)
            
            cDr(nCellStr:nCellEnd) = cDrDp(nCellStr:nCellEnd)*cDp(nCellStr:nCellEnd)
            
            do i=nCellStr,nCellEnd
                lFace = i-1
                rFace = i
                gP = 0.5d0!*(1.d0+dsign(1.d0,fU(rFace)))
                gN = 0.5d0!*(1.d0-dsign(1.d0,fU(rFace)))
                cDu(i) = -cDt(i)/cRho(i)/cVol(i) &
                        *area(rFace)*(gP*cDp(i)+gN*cDp(i+1)) 
                
                gP = 0.5d0!*(1.d0+dsign(1.d0,fU(lFace)))
                gN = 0.5d0!*(1.d0-dsign(1.d0,fU(lFace)))
                cDu(i) = cDu(i) +cDt(i)/cRho(i)/cVol(i) &
                        *area(lFace)*(gP*cDp(i)+gN*cDp(i-1))
                
            enddo
            
            do i=nCellStr,nCellEnd
                cP(i) = cP(i) + cDp(i)
                cU(i) = cU(i) + cDu(i)
                cRho(i) = cRho(i) + cDr(i)
            enddo
            
            !=======================================================
            
            resiP = 0.d0
            resiU = 0.d0
            resiR = 0.d0
            do i=nCellStr,nCellEnd
                resiP=resiP+cDp(i)**2.d0
                resiU=resiU+cDu(i)**2.d0
                resiR=resiR+cDr(i)**2.d0
            enddo
            resiP = dsqrt(resiP)
            resiU = dsqrt(resiU)
            resiR = dsqrt(resiR)
            
            !!> Pressure Residual
            !resiP = 0.d0
            !maxP = -1.d8
            !minP = 1.d8
            !do i=nCellStr,nCellEnd
            !    resiP = resiP + cDp(i)**2.d0
            !    maxP = dmax1(maxP,cP(i))
            !    minP = dmin1(minP,cP(i))
            !enddo
            !resiP = dsqrt(resiP) / (dmax1(maxP,0.d0)-dmin1(minP,0.d0)+1.d-15)
            !    
            !!> Velocity Residual
            !resiU = 0.d0
            !maxU = -1.d8
            !minU = 1.d8
            !do i=nCellStr,nCellEnd
            !    resiU = resiU + cDU(i)**2.d0
            !    maxU = dmax1(maxU,cU(i))
            !    minU = dmin1(minU,cU(i))
            !enddo
            !resiU = dsqrt(resiU) / (dmax1(maxU,0.d0)-dmin1(minU,0.d0)+1.d-15)
            !
            !!> Density Residual
            !resiR = 0.d0
            !maxR = -1.d8
            !minR = 1.d8
            !do i=nCellStr,nCellEnd
            !    resiR = resiR + cDr(i)**2.d0
            !    maxR = dmax1(maxU,cRho(i))
            !    minR = dmin1(minU,cRho(i))
            !enddo
            !resiR = dsqrt(resiR) / (dmax1(maxR,0.d0)-dmin1(minR,0.d0)+1.d-15)
            
            
            write(*,100) iterSIMPLE,resiP,resiU,resiR
            !stop
            
            !=======================================================

        enddo Segregated !> End iteration SIMPLE
        
        
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
    
100 format(i7,1x,6(1pe10.2))
 
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

    
    
    
!============================================================
!> Inverse matrix
!! Method: Based on Doolittle LU factorization for Ax=b
!! Original code by Alex G. December 2009
!! modified by Jongchan Kim 2018
!!-----------------------------------------------------------
!! input ...
!! a(n,n) - array of coefficients for matrix A
!! n      - dimension
!! output ...
!! a(n,n) - inverse matrix of A
!! comments ...
!! the original matrix a(n,n) will be destroyed 
!! during the calculation
!!===========================================================
subroutine matInverse(a,n)

    implicit none
    integer :: n
    real(8) :: a(n,n)
    real(8) :: L(n,n), U(n,n), b(n), d(n), x(n)
    real(8) :: coeff
    integer :: i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.d0
    U=0.d0
    b=0.d0

    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.d0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
        b(k)=1.d0
        d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do
        ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
        end do
        ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
            a(i,k) = x(i)
        end do
        b(k)=0.d0
    end do
end subroutine matInverse
   