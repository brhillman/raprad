        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 18 19:23:54 2014
        MODULE RAP_TWOSTR_SW__genmod
          INTERFACE 
            SUBROUTINE RAP_TWOSTR_SW(SOLCONST,NLAYERS,AMU0,W0_TEMP,     &
     &G0_TEMP,SURALB,TAUL_TEMP,GAU_WT,DIS2SUN)
              REAL(KIND=8) :: SOLCONST
              INTEGER(KIND=4) :: NLAYERS
              REAL(KIND=4) :: AMU0
              REAL(KIND=4) :: W0_TEMP(*)
              REAL(KIND=4) :: G0_TEMP(*)
              REAL(KIND=4) :: SURALB
              REAL(KIND=4) :: TAUL_TEMP(*)
              REAL(KIND=8) :: GAU_WT
              REAL(KIND=8) :: DIS2SUN
            END SUBROUTINE RAP_TWOSTR_SW
          END INTERFACE 
        END MODULE RAP_TWOSTR_SW__genmod
