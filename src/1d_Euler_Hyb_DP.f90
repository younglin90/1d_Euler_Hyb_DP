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
    use modInput
    !use modFaceValues
    implicit none
    LOGICAL:: there

    
    INQUIRE (FILE = 'residual.txt', EXIST = there)
    IF (there) then
        OPEN(20, FILE = 'residual.txt', STATUS = 'REPLACE')
    ELSE
        OPEN(20, FILE = 'residual.txt', STATUS = 'NEW')
    END IF
    close(20)
    
    CALL getInput()
    
    CALL getAllocate()
    
    CALL getInputSpec()
    
    
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
          cP(i)   = 1.000d0
          cU(i)   = 0.000d0
          cT(i)   = 0.003486055811614
      else
          cP(i)   = 0.100d0
          cU(i)   = 0.000d0
          cT(i)   = 0.002788844649291
      endif
    enddo
    !--------- BC
    CALL getBC(leftBC='supBC',rightBC='supBC')
    
    !--------- Ht, Et, T, c from EOS
    CALL getEOS(eosRho='on',eosHt='on',eosDRDP='on',eosDRDT='on',&
                eosDHDP='on',eosDHDT='on',eosC='on')
    time = 0.d0
    nstep = 0
    ntimestep=0
 
    do while( time <= time_max )
        
        ntimestep=ntimestep+1
        time = time + timestep
        
        !=======================================================
        !> Save old Consvervative Variables
        !> Save old Primitive Variables
        do i=nTCellStr,nTCellEnd
            oldCons(i,1) = cRho(i) 
            oldCons(i,2) = cRho(i)*cU(i) 
            oldCons(i,3) = cRho(i)*cHt(i)-cP(i)
            
            oldcRho(i) = cRho(i)
            oldcU(i)   = cU(i)
            oldcP(i)   = cP(i)
            oldcT(i)   = cT(i)
            oldcHt(i)  = cHt(i)
        enddo
        
        
        !======================== 1 ==========================
        !> Continuity Equation -> Density Prediction
        CALL getContinuity('explicit')
        

        SIMPLE: do iterSIMPLE = 1,iterMaxSIMPLE
        
            !======================== 2 ==========================
            !> Momentum Equation -> Velocity Prediction
            CALL getMomentum('impCoeff')
            

            PISO: do iterPISO = 1,iterMaxPISO
                
                nstep=nstep+1
        
                !======================== 3 ==========================
                !> Energy Equation -> Temperature Prediction
                CALL getEnergy('explicit')

                !======================== 4 ==========================
                !> Equation of State -> update Thermodynamic variables
                CALL getEOS(eosRho='on',eosHt='on',eosDRDP='on')
                
                !======================== 5 ==========================
                !> Pressure Equation (Continuity + Momentum) 
                !> => update Pressure
                CALL getPressure('implicit')

                !======================== 7 ==========================
                !> Momentum Equation => update Velocity
                CALL getMomentum('dP_to_U')
                
                !======================== 6 ==========================
                !> Equation of State => update Density
                CALL getEOS(eosRho='on',eosHt='on',eosDHDT='on')
        
                !====================================================
                !> Monitor of Residuals
                CALL seeMonitor


            enddo PISO !> End iteration

            !====================================================
            CALL getEnergy('explicit')
            
            !====================================================
            CALL getContinuity('explicit')
    
    
        enddo SIMPLE !> End iteration
        
        
        !=======================================================
        !> Boundary Conditions
        CALL getBC(leftBC='supBC',rightBC='supBC')
        
        !=======================================================
        !> Equation of State => update thermodynamics
        CALL getEOS(eosHt='on',eosDRDP='on')
        
    
    enddo
    
    
    
    !/////////////////////////////////
    open(20,file='./plot/result.dat')!,status='old',position='append')
    do i=nCellStr,nCellEnd
        WRITE( 20, '(20(F20.9, 1X))' ) cX(i),cRho(i),cU(i),cP(i),cT(i)
    enddo
    close(20)
    !/////////////////////////////////
    
    print*,'   =================== program end ================='

end program Euler

    
    
    