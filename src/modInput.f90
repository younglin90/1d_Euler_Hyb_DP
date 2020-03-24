
module modInput

    use modValues
    implicit none
    
    contains
    
    subroutine getInput
        implicit none
        character(len=30) :: input_name
        character(80) :: data_name,suffix_num
    
        namelist/input/ &
            nCell,time_max,eosGamma,eosCp,len, &
            iterMaxSIMPLE,iterMaxPISO,underRelaxFactP, &
            underRelaxFactT,underRelaxFactR,underRelaxFactU, &
            timestep,chReconP,chReconU,chReconT,chReconR,&
            chFluxP,chFluxU,chFluxT,chFluxR
    
        input_name = 'config.in'
        open (11,file =input_name,status='old')
        read (11,input)
        close(11)
    endsubroutine
    

    
endmodule
    