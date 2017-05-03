        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 12:49:18 2014
        MODULE OPPR1__genmod
          INTERFACE 
            SUBROUTINE OPPR1(NLAYERS,PTEMPG,TAULAM,SLOPE,PTEMP,BANDWIDTH&
     &,PLANKBND,GAU_WT,PLANKLEV)
              INTEGER(KIND=4) :: NLAYERS
              REAL(KIND=8) :: PTEMPG
              REAL(KIND=4) :: TAULAM(201)
              REAL(KIND=8) :: SLOPE(201)
              REAL(KIND=8) :: PTEMP(201)
              REAL(KIND=8) :: BANDWIDTH(2)
              REAL(KIND=8) :: PLANKBND
              REAL(KIND=8) :: GAU_WT(201)
              REAL(KIND=8) :: PLANKLEV(201)
            END SUBROUTINE OPPR1
          END INTERFACE 
        END MODULE OPPR1__genmod
