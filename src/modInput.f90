
module modInput

    use modValues
    implicit none
    
    contains
    
    !==============================
    subroutine getInput
        implicit none
        character(len=30) :: input_name
        character(80) :: data_name,suffix_num
    
        namelist/input/ &
            nCell,time_max,len, &
            iterMaxSIMPLE,iterMaxPISO,underRelaxFactP, &
            underRelaxFactT,underRelaxFactR,underRelaxFactU, &
            timestep,chReconP,chReconU,chReconT,chReconR,&
            chFluxP,chFluxU,chFluxT,chFluxR,nspe,chReconY
    
        input_name = 'config.in'
        open (11,file =input_name,status='old')
        read (11,input)
        close(11)
    endsubroutine
    
    
    !==============================
    subroutine getInputSpec
        implicit none
        character(len=30) :: input_name
        character(80) :: data_name,suffix_num
    
        namelist/inputSpec/ &
            chEOS,eosCp,eosGamma,eosCv,eosB,eosQ,eosPinf
    
        input_name = 'configSpec.in'
        open (11,file =input_name,status='old')
        read (11,inputSpec)
        close(11)
        
    endsubroutine
    
endmodule
    