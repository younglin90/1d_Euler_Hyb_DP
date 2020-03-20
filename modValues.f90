
module modValues

    implicit none


    real(8) :: CFL, time_max, eosGamma, eosCp, len, dx, &
     time, charVel, wL(3), wR(3), vecSf, dPN, avgDp, underRelaxFactP, &
     timestep, timestepMax, resiP, fVol, gP, gF, resiU, &
     maxP,minP,maxU,minU,segAlp1,segAlp2,gN,resiR,maxR,minR,&
     saveRho,oldLfUx,oldRfUx,resiT,savecT,savecDummy,&
     underRelaxFactT,underRelaxFactR,underRelaxFactU
    real(8), dimension(:), allocatable :: cX, cRho, cU, cP, cT, &
     cVol, cDt, area, fMdot, fU, fP, fT, cAcRU, cAfRU, cBcRU, &
     cAcUU, cAfUU, cAfUP, cBcUU, cDrDp, cHt, &
     fDrDp, fRho, fFirPE, fAcUU, fSecPE, fDp_fS, &
     cAcPDp, cAfPR, cAfPDp, cBcPDp, cDp, &
     fHt, cAcHH, cAfHU, cBcHH, fThrPE, fDpDn, fDp, cDU, &
     fUx, cEt, fEt, cC, cGECon, cGEMom, cGEEnr, &
     fGEConFlux, fGEMomFlux, fGEEnrFlux, segBmat,cDr, &
     oldcRho,oldcU,oldcP,fgP,fgN,cDEtDT,cDrDT,oldcT,oldcEt,oldcHt,&
    BmatMom1,BmatMom2,BmatPres,coeffMom1,coeffMom2,fcoeffMom1,fcoeffMom2, &
        saveOneValue,BmatEnr,BmatCon
    real(8), dimension(:,:), allocatable :: segAmat,oldCons,oldPrim,&
         AmatPres,AmatMom,AmatEnr,AmatCon
    integer :: nCell, nCellStr, nCellEnd, nTCellStr, nTCellEnd, &
     nFace, nFaceStr, nFaceEnd, i, nstep, lFace, rFace, iterSIMPLE, &
     iterMaxSIMPLE, lCell, rCell, iterPISO, iterMaxPISO, iterDummy
    
    contains
    
    subroutine initValues
    
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
                  cDr(nTCellStr:nTCellEnd), &
                  oldCons(nTCellStr:nTCellEnd,3), &
                  oldPrim(nTCellStr:nTCellEnd,3), &
                  oldcRho(nTCellStr:nTCellEnd), &
                  oldcU(nTCellStr:nTCellEnd), &
                  oldcP(nTCellStr:nTCellEnd), &
                  cDEtDT(nTCellStr:nTCellEnd), &
                  cDrDT(nTCellStr:nTCellEnd), &
                  oldcT(nTCellStr:nTCellEnd), &
                  oldcEt(nTCellStr:nTCellEnd), &
                  oldcHt(nTCellStr:nTCellEnd), &
                  coeffMom1(nTCellStr:nTCellEnd), &
                  coeffMom2(nTCellStr:nTCellEnd), &
                  saveOneValue(nTCellStr:nTCellEnd) )
        cDP(:) = 0.d0
        cDu(:) = 0.d0
    
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
                  fGEEnrFlux(nFaceStr:nFaceEnd), &
                  fgP(nFaceStr:nFaceEnd), &
                  fgN(nFaceStr:nFaceEnd), &
                  fcoeffMom1(nFaceStr:nFaceEnd), &
                  fcoeffMom2(nFaceStr:nFaceEnd)   )
    
        allocate( AmatPres(nCellEnd,nCellEnd), &
                  BmatPres(nCellEnd), &
                  AmatMom(nCellEnd,nCellEnd), &
                  BmatMom1(nCellEnd), &
                  BmatMom2(nCellEnd), &
                  AmatEnr(nCellEnd,nCellEnd), &
                  BmatEnr(nCellEnd), &
                  AmatCon(nCellEnd,nCellEnd), &
                  BmatCon(nCellEnd) )
    endsubroutine
    
endmodule
    