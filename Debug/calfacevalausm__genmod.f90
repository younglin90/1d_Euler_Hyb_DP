        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 20 21:45:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALFACEVALAUSM__genmod
          INTERFACE 
            SUBROUTINE CALFACEVALAUSM(WL,WR,GAMMA,MLR,ULR,PLR)
              REAL(KIND=8), INTENT(IN) :: WL(3)
              REAL(KIND=8), INTENT(IN) :: WR(3)
              REAL(KIND=8), INTENT(IN) :: GAMMA
              REAL(KIND=8), INTENT(OUT) :: MLR
              REAL(KIND=8), INTENT(OUT) :: ULR
              REAL(KIND=8), INTENT(OUT) :: PLR
            END SUBROUTINE CALFACEVALAUSM
          END INTERFACE 
        END MODULE CALFACEVALAUSM__genmod
