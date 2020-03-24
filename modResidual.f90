
module modResidual

    use modValues
    implicit none
    
    
    contains
    
    subroutine seeMonitor
        implicit none
        
        real(8) :: maxvalue,minvalue
                        
                !resiP = 0.d0
                !resiU = 0.d0
                !resiR = 0.d0
                !resiT = 0.d0
                !do i=nCellStr,nCellEnd
                !    resiP=resiP+cDp(i)**2.d0
                !    resiU=resiU+cDu(i)**2.d0
                !    resiR=resiR+cDr(i)**2.d0
                !    resiT=resiT+cDT(i)**2.d0
                !enddo
                !resiP = dsqrt(resiP)
                !resiU = dsqrt(resiU)
                !resiR = dsqrt(resiR)
                !resiT = dsqrt(resiT)
            
                !> Pressure Residual
                resiP = 0.d0
                maxvalue = -1.d8
                minvalue = 1.d8
                do i=nCellStr,nCellEnd
                    resiP = resiP + cDp(i)**2.d0
                    maxvalue = dmax1(maxvalue,cP(i))
                    minvalue = dmin1(minvalue,cP(i))
                enddo
                resiP = dsqrt(resiP) / (dmax1(maxvalue,0.d0)-dmin1(minvalue,0.d0)+1.d-15)
                    
                !> Velocity Residual
                resiU = 0.d0
                maxvalue = -1.d8
                minvalue = 1.d8
                do i=nCellStr,nCellEnd
                    resiU = resiU + cDU(i)**2.d0
                    maxvalue = dmax1(maxvalue,cU(i))
                    minvalue = dmin1(minvalue,cU(i))
                enddo
                resiU = dsqrt(resiU) / (dmax1(maxvalue,0.d0)-dmin1(minvalue,0.d0)+1.d-15)
                
                !> Density Residual
                resiR = 0.d0
                maxvalue = -1.d8
                minvalue = 1.d8
                do i=nCellStr,nCellEnd
                    resiR = resiR + cDR(i)**2.d0
                    maxvalue = dmax1(maxvalue,cRho(i))
                    minvalue = dmin1(minvalue,cRho(i))
                enddo
                resiR = dsqrt(resiR) / (dmax1(maxvalue,0.d0)-dmin1(minvalue,0.d0)+1.d-15)
                
                !> Temp. Residual
                resiT = 0.d0
                maxvalue = -1.d8
                minvalue = 1.d8
                do i=nCellStr,nCellEnd
                    resiT = resiT + cDT(i)**2.d0
                    maxvalue = dmax1(maxvalue,cT(i))
                    minvalue = dmin1(minvalue,cT(i))
                enddo
                resiT = dsqrt(resiT) / (dmax1(maxvalue,0.d0)-dmin1(minvalue,0.d0)+1.d-15)
            
            !if( mod(nstep,4) == 0 ) then
                write(*,100) nstep,iterSIMPLE,iterPISO,resiP,resiU,resiR,resiT
                !stop
                !/////////////////////////////////
                open(20,file='residual.txt',status='old',position='append')
                WRITE( 20, '(i5,1X,20(F20.9, 1X))' ) nstep,log(resiP),&
                    log(resiU),log(resiR),log(resiT)
                close(20)
                !/////////////////////////////////
            if( mod(nstep,1000) == 0 ) then
                !OPEN(20, FILE = 'residual.txt', STATUS = 'REPLACE')
                close(20)
            endif
            
            !endif
                !=======================================================
    
100 format(i7,1x,i7,1x,i7,1x,6(1pe10.2))
 
                
    endsubroutine
    

    
endmodule
    