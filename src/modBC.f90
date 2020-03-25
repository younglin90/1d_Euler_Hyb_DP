
module modBC

    use modValues
    implicit none
    
    contains
    
    !========= Boundary Conditions ============
    subroutine getBC(leftBC,rightBC)
        implicit none
        character(len=*),optional,intent(in) :: leftBC,rightBC
            
        if(present(leftBC)) then
            if( trim(leftBC) == 'supBC' ) then
                cRho(nCellStr-1) = cRho(nCellStr)
                cU(nCellStr-1)   = cU(nCellStr)
                cP(nCellStr-1)   = cP(nCellStr)
                cT(nCellStr-1)   = cT(nCellStr)
                cY(nCellStr-1,1:nspe)   = cY(nCellStr,1:nspe)
            elseif( trim(leftBC) == 'wallBC' ) then
                
                
            endif
        endif
        
        if(present(leftBC)) then
            if( trim(leftBC) == 'supBC' ) then
                cRho(nCellEnd+1) = cRho(nCellEnd)
                cU(nCellEnd+1)   = cU(nCellEnd)
                cP(nCellEnd+1)   = cP(nCellEnd)
                cT(nCellEnd+1)   = cT(nCellEnd)
                cY(nCellEnd+1,1:nspe)   = cY(nCellEnd,1:nspe)
            elseif( trim(leftBC) == 'wallBC' ) then
                
            endif
        endif

    endsubroutine
    

    
endmodule
    