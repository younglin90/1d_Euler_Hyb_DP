
module modEOS

    use modValues
    implicit none
    
    contains
    
    !========= Pressure Equations ============
    subroutine getEOS(eosC,eosRho,eosP,eosT,eosEt,eosHt,eosDRDP)
        implicit none
        character(len=*),optional,intent(in) :: eosC,eosRho,eosP,eosT,&
            eosEt,eosHt,eosDRDP
        
        do i=nTCellStr,nTCellEnd
            
            if(present(eosRho)) cRho(i) = cP(i)/(eosCp*(1.d0-1.d0/eosGamma))/cT(i)
            if(present(eosC)) cC(i) = dsqrt(eosGamma*cP(i)/cRho(i))
            if(present(eosP)) cP(i) = cRho(i)*(eosCp*(1.d0-1.d0/eosGamma))*cT(i)
            if(present(eosT)) cT(i) = cP(i)/(eosCp*(1.d0-1.d0/eosGamma))/cRho(i)
            if(present(eosEt)) cEt(i) = eosCp*cT(i) + 0.5d0*(cU(i)**2.d0) - cP(i)/cRho(i)
            if(present(eosHt)) cHt(i) = eosCp*cT(i) + 0.5d0*(cU(i)**2.d0)
            if(present(eosDRDP)) cDRDP(i) = 1.d0/eosCp/(1.d0-1.d0/eosGamma)/cT(i)
            
        enddo
        

    endsubroutine
    

    
endmodule
    