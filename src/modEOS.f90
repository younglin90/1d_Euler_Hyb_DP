
module modEOS

    use modValues
    implicit none
    integer :: nq
    
    contains
    
    !========= get EOS ============
    subroutine getEOS(eosC,eosRho,eosHt,eosDRDP,eosDRDT,eosDHDP,eosDHDT)
        implicit none
        character(len=*),optional,intent(in) :: eosC,eosRho,&
            eosHt,eosDRDP,eosDRDT,eosDHDP,eosDHDT
        real(8) :: eosRhon(nspe),eosCn(nspe),eosHtn(nspe),&
            eosDRDPn(nspe),eosDRDTn(nspe),eosDHDPn(nspe),eosDHDTn(nspe),&
            eosAlp(nspe),eosY(nspe)
        
        do i=nTCellStr,nTCellEnd
            
            
            do nq=1,nspe
                
                SELECTCASE(chEOS(nq))
                CASE('Ideal')
                        CALL getEOSIdealRho(cP(i),cT(i),eosRhon(nq))
                    if(present(eosHt))   &
                        CALL getEOSIdealHt(cU(i),cT(i),eosHtn(nq))
                    if(present(eosDRDP)) &
                        CALL getEOSIdealDRDP(cT(i),eosDRDPn(nq))
                    if(present(eosDRDT)) &
                        CALL getEOSIdealDRDT(cP(i),cT(i),eosDRDTn(nq))
                    if(present(eosDHDP)) &
                        CALL getEOSIdealDHDP(eosDHDPn(nq))
                    if(present(eosDHDT)) &
                        CALL getEOSIdealDHDT(eosDHDTn(nq))
                CASE('NASG')
                        CALL getEOSNASGRho(cP(i),cT(i),eosRhon(nq))
                    if(present(eosHt))   &
                        CALL getEOSNASGHt(cP(i),cU(i),cT(i),eosHtn(nq))
                    if(present(eosDRDP)) &
                        CALL getEOSNASGDRDP(cP(i),cT(i),eosDRDPn(nq))
                    if(present(eosDRDT)) &
                        CALL getEOSNASGDRDT(cP(i),cT(i),eosDRDTn(nq))
                    if(present(eosDHDP)) &
                        CALL getEOSNASGDHDP(eosDHDPn(nq))
                    if(present(eosDHDT)) &
                        CALL getEOSNASGDHDT(eosDHDTn(nq))
                ENDSELECT
                
            enddo
            
            
            eosY(1:nspe) = cY(i,1:nspe)

            cRho(i) = 1.d0 / sum(eosY(1:nspe)/eosRhon(1:nspe))
            
            forall(nq=1:nspe-1) &
                eosAlp(nq) = dmax1( 0.d0 , dmin1( 1.d0, cRho(i)*eosY(nq)/eosRhon(nq) ) )
            eosAlp(nspe) = 1.d0 - sum( eosAlp(1:nspe-1) )
            eosAlp(nspe) = dmax1( 0.d0 , dmin1( 1.d0, eosAlp(nspe) ) )
            
            if(present(eosHt))   &
                cHt(i) = sum(eosY(1:nspe)*eosHtn(1:nspe))
            if(present(eosDRDP)) &
                cDRDP(i) = sum(eosAlp(1:nspe)*eosDRDPn(1:nspe))
            if(present(eosDRDT)) &
                cDRDT(i) = sum(eosAlp(1:nspe)*eosDRDTn(1:nspe))
            if(present(eosDHDP)) &
                cDHDP(i) = sum(eosY(1:nspe)*eosDHDPn(1:nspe))
            if(present(eosDHDT)) &
                cDHDT(i) = sum(eosY(1:nspe)*eosDHDTn(1:nspe))
            if(present(eosC)) &
                cC(i) = dsqrt( 1.d0/(cDRDP(i) &
		          + 1.d0/cRho(i)*cDRDT(i)/cDHDT(i)*(1.d0-cRho(i)*cDHDP(i))) )
            
        enddo

    endsubroutine
    
    
    
    !========= get speed of sound ============
    subroutine getRCH(inpP,inpU,inpT,inpY, resultRho,resultC,resultHt)
        implicit none
        real(8) :: tempDRDP,tempDRDT,tempDHDP,tempDHDT,&
            inpY(nspe),inpP,inpU,inpT, resultRho,resultC,resultHt
        real(8) :: eosRhon(nspe),eosCn(nspe),eosHtn(nspe),&
            eosDRDPn(nspe),eosDRDTn(nspe),eosDHDPn(nspe),eosDHDTn(nspe),&
            eosAlp(nspe),eosY(nspe)
        
        do nq=1,nspe
                
            SELECTCASE(chEOS(nq))
            CASE('Ideal')
                CALL getEOSIdealRho(inpP,inpT,eosRhon(nq))
                CALL getEOSIdealHt(inpU,inpT,eosHtn(nq))
                CALL getEOSIdealDRDP(inpT,eosDRDPn(nq))
                CALL getEOSIdealDRDT(inpP,inpT,eosDRDTn(nq))
                CALL getEOSIdealDHDP(eosDHDPn(nq))
                CALL getEOSIdealDHDT(eosDHDTn(nq))
            CASE('NASG')
                CALL getEOSNASGRho(inpP,inpT,eosRhon(nq))
                CALL getEOSNASGHt(inpP,inpU,inpT,eosHtn(nq))
                CALL getEOSNASGDRDP(inpP,inpT,eosDRDPn(nq))
                CALL getEOSNASGDRDT(inpP,inpT,eosDRDTn(nq))
                CALL getEOSNASGDHDP(eosDHDPn(nq))
                CALL getEOSNASGDHDT(eosDHDTn(nq))
            ENDSELECT
                
        enddo
            
        forall(nq=1:nspe) eosY(nq) = inpY(nq)
            
        resultRho = 1.d0 / sum(eosY(1:nspe)/eosRhon(1:nspe))
            
        forall(nq=1:nspe-1) &
            eosAlp(nq) = dmax1( 0.d0 , dmin1( 1.d0, cRho(i)*eosY(nq)/eosRhon(nq) ) )
        eosAlp(nspe) = 1.d0 - sum( eosAlp(1:nspe-1) )
        eosAlp(nspe) = dmax1( 0.d0 , dmin1( 1.d0, eosAlp(nspe) ) )
            
        resultHt = sum(eosY(1:nspe)*eosHtn(1:nspe))
        tempDRDP = sum(eosAlp(1:nspe)*eosDRDPn(1:nspe))
        tempDRDT = sum(eosAlp(1:nspe)*eosDRDTn(1:nspe))
        tempDHDP = sum(eosY(1:nspe)*eosDHDPn(1:nspe))
        tempDHDT = sum(eosY(1:nspe)*eosDHDTn(1:nspe))
        resultC = dsqrt( 1.d0/(tempDRDP &
		    + 1.d0/resultRho*tempDRDT/tempDHDT*(1.d0-resultRho*tempDHDP)) )
            
    endsubroutine
    
    
    !================== Ideal EOS ===============================
    subroutine getEOSIdealRho(inpP,inpT,eosRhon)
        implicit none
        real(8), intent(inout) :: inpP,inpT,eosRhon
        
        eosRhon = inpP/(eosCp(nq)*(1.d0-1.d0/eosGamma(nq)))/inpT
    
    endsubroutine
    subroutine getEOSIdealHt(inpU,inpT,eosHtn)
        implicit none
        real(8), intent(inout) :: inpU,inpT,eosHtn
        
        eosHtn = eosCp(nq)*cT(i) + 0.5d0*(inpU**2.d0)
    
    endsubroutine
    subroutine getEOSIdealDRDP(inpT,eosDRDPn)
        implicit none
        real(8), intent(inout) :: inpT,eosDRDPn
        
        eosDRDPn = 1.d0/(eosCp(nq)*(1.d0-1.d0/eosGamma(nq)))/inpT
    
    endsubroutine
    subroutine getEOSIdealDRDT(inpP,inpT,eosDRDTn)
        implicit none
        real(8), intent(inout) :: inpP,inpT,eosDRDTn
        
        eosDRDTn = -inpP/(eosCp(nq)*(1.d0-1.d0/eosGamma(nq)))/inpT/inpT
    
    endsubroutine
    subroutine getEOSIdealDHDP(eosDHDPn)
        implicit none
        real(8), intent(inout) :: eosDHDPn
        
        eosDHDPn = 0.d0
    
    endsubroutine
    subroutine getEOSIdealDHDT(eosDHDTn)
        implicit none
        real(8), intent(inout) :: eosDHDTn
        
        eosDHDTn = eosCp(nq)
    
    endsubroutine
    
    !================== NASG EOS ===============================
    subroutine getEOSNASGRho(inpP,inpT,eosRhon)
        implicit none
        real(8), intent(inout) :: inpP,inpT,eosRhon
        
        eosRhon = 1.d0/( (eosGamma(nq)-1.d0)&
            *eosCv(nq)*inpT/(inpP+eosPinf(nq))+eosB(nq) )
    
    endsubroutine
    subroutine getEOSNASGHt(inpP,inpU,inpT,eosHtn)
        implicit none
        real(8), intent(inout) :: inpP,inpU,inpT,eosHtn
        
        eosHtn = eosGamma(nq)*eosCv(nq)*inpT + eosB(nq)*inpP &
            + eosQ(nq) + 0.5d0*(inpU**2.d0)
    
    endsubroutine
    subroutine getEOSNASGDRDP(inpP,inpT,eosDRDPn)
        implicit none
        real(8), intent(inout) :: inpP,inpT,eosDRDPn
        real(8) :: eosrhoii
        
        eosrhoii = 1.d0/( (eosGamma(nq)-1.d0)&
            *eosCv(nq)*inpT/(inpP+eosPinf(nq))+eosB(nq) )
        eosDRDPn = eosrhoii**2.d0*(1.d0/eosrhoii-eosB(nq))/(inpP+eosPinf(nq)) 
    
    endsubroutine
    subroutine getEOSNASGDRDT(inpP,inpT,eosDRDTn)
        implicit none
        real(8), intent(inout) :: inpP,inpT,eosDRDTn
        real(8) :: eosrhoii
        
        eosrhoii = 1.d0/( (eosGamma(nq)-1.d0)&
            *eosCv(nq)*inpT/(inpP+eosPinf(nq))+eosB(nq) )
        eosDRDTn = -eosrhoii**2.d0*(1.d0/eosrhoii-eosB(nq))/inpT
    
    endsubroutine
    subroutine getEOSNASGDHDP(eosDHDPn)
        implicit none
        real(8), intent(inout) :: eosDHDPn
        
        eosDHDPn = eosB(nq)
    
    endsubroutine
    subroutine getEOSNASGDHDT(eosDHDTn)
        implicit none
        real(8), intent(inout) :: eosDHDTn
        
        eosDHDTn = eosGamma(nq)*eosCv(nq)
    
    endsubroutine
    
endmodule
    