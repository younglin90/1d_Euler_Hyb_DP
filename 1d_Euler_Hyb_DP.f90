!======================
! 1D euler eq.
!=====================
program Euler

    use modValues
    use modPressure
    use modContinuity
    use modEnergy
    use modMomentum
    use modResidual
    use modEOS
    use modBC
    implicit none

    
    !CFL = 0.1
    nCell = 100
    time_max = 0.2 !sec
    eosGamma = 1.4 ! air
    eosCp = 1004.d0
    len = 1.0 ! m, x-length
    iterMaxSIMPLE = 3
    iterMaxPISO = 5
    underRelaxFactP = 0.3d0
    underRelaxFactT = 0.3d0
    underRelaxFactR = 0.5d0
    underRelaxFactU = 0.5d0
    !timestepMax = 1.d-2
    timestep = 3.d-3
    
    CALL initValues()
    
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
      if( cX(i) < 0.5d0 ) then
          cRho(i) = 1.000d0
          cU(i)   = 0.000d0
          cP(i)   = 1.000d0
      else
          cRho(i) = 0.125d0
          cU(i)   = 0.000d0
          cP(i)   = 0.100d0
      endif
    enddo
    !--------- BC
    CALL getBC(leftBC='supBC',rightBC='supBC')
    
    !--------- Ht, Et, T, c from EOS
    CALL getEOS(eosT='on',eosC='on',eosEt='on',eosHt='on',eosDRDP='on')
    
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
        
        
        
        !=======================================================
        !> Save old Consvervative Variables
        !> Save old Primitive Variables
        do i=nTCellStr,nTCellEnd
            oldCons(i,1) = cRho(i) 
            oldCons(i,2) = cRho(i)*cU(i) 
            oldCons(i,3) = cRho(i)*cEt(i) 
            
            oldcRho(i) = cRho(i)
            oldcU(i)   = cU(i)
            oldcP(i)   = cP(i)
            oldcT(i)   = cT(i)
            oldcEt(i)  = cEt(i)
            oldcHt(i)  = cHt(i)
        enddo
        
        
        !======================== 1 ==========================
        !> Continuity Equation -> Density Prediction
        CALL getContinuity('implicit')
        

        SIMPLE: do iterSIMPLE = 1,iterMaxSIMPLE
        
            !======================== 2 ==========================
            !> Momentum Equation -> Velocity Prediction
            CALL getMomentum('impCoeff')
            

            PISO: do iterPISO = 1,iterMaxPISO
        
                !======================== 3 ==========================
                !> Energy Equation -> Temperature Prediction
                CALL getEnergy('implicit')

                !======================== 4 ==========================
                !> Equation of State -> update Thermodynamic variables
                CALL getEOS(eosRho='on',eosEt='on',eosHt='on')
                
                !======================== 5 ==========================
                !> Pressure Equation (Continuity + Momentum) 
                !> => update Pressure
                CALL getPressure('implicit')

                !======================== 7 ==========================
                !> Momentum Equation => update Velocity
                CALL getMomentum('dp_to_u')
                
                !======================== 6 ==========================
                !> Equation of State => update Density
                CALL getEOS(eosRho='on',eosEt='on',eosHt='on')

                !=======================
                !CALL getEnergy('implicit')
                !CALL getContinuity('implicit')
        
                !====================================================
                !> Monitor of Residuals
                CALL seeMonitor


            enddo PISO !> End iteration
            

            !====================================================
            CALL getEnergy('implicit')
            
            !====================================================
            CALL getContinuity('implicit')
    
    
        enddo SIMPLE !> End iteration
        
        
        !=======================================================
        !> Boundary Conditions
        CALL getBC(leftBC='supBC',rightBC='supBC')
        
        !=======================================================
        !> Equation of State => update thermodynamics
        CALL getEOS(eosEt='on',eosHt='on',eosDRDP='on')
                
                
        time = time + timestep
        
    
    enddo
    
    
    
    !/////////////////////////////////
    open(20,file='a.dat')!,status='old',position='append')
    do i=nCellStr,nCellEnd
        WRITE( 20, '(20(F20.9, 1X))' ) cX(i),cRho(i),cU(i),cP(i),cT(i)
    enddo
    close(20)
    !/////////////////////////////////
    
    print*,'   =================== program end ================='

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
   