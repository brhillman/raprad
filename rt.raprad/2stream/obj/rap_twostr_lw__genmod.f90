        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 18 18:42:54 2014
        MODULE RAP_TWOSTR_LW__genmod
          INTERFACE 
            SUBROUTINE RAP_TWOSTR_LW(SOLCONST,NLAYERS,AMU0,W0_TEMP,     &
     &G0_TEMP,SURALB,TAUL_TEMP,GAU_WT_TEMP,DIS2SUN,PLANKBND,            &
     &PLANKLAY_TEMP,PLANKLEV_TEMP,BANDWIDTH)
              REAL(KIND=8) :: SOLCONST
              INTEGER(KIND=4) :: NLAYERS
              REAL(KIND=4) :: AMU0
              REAL(KIND=4) :: W0_TEMP(201)
              REAL(KIND=4) :: G0_TEMP(201)
              REAL(KIND=4) :: SURALB
              REAL(KIND=4) :: TAUL_TEMP(201)
              REAL(KIND=8) :: GAU_WT_TEMP(201)
              REAL(KIND=8) :: DIS2SUN
              REAL(KIND=8) :: PLANKBND
              REAL(KIND=8) :: PLANKLAY_TEMP(201)
              REAL(KIND=8) :: PLANKLEV_TEMP(201)
              REAL(KIND=8) :: BANDWIDTH(2)
            END SUBROUTINE RAP_TWOSTR_LW
          END INTERFACE 
        END MODULE RAP_TWOSTR_LW__genmod
