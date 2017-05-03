        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 12:49:43 2014
        MODULE TWOSTRSW__genmod
          INTERFACE 
            SUBROUTINE TWOSTRSW(NLAYER,TAUL,W0,G0,RSFX,B1,B2,EL1,EL2,EM1&
     &,EM2,AF,BF,EF,AK)
              INTEGER(KIND=4) :: NLAYER
              REAL(KIND=4) :: TAUL(*)
              REAL(KIND=4) :: W0(*)
              REAL(KIND=4) :: G0(*)
              REAL(KIND=4) :: RSFX
              REAL(KIND=8) :: B1(201)
              REAL(KIND=8) :: B2(201)
              REAL(KIND=8) :: EL1(201)
              REAL(KIND=8) :: EL2(201)
              REAL(KIND=8) :: EM1(201)
              REAL(KIND=8) :: EM2(201)
              REAL(KIND=8) :: AF(402)
              REAL(KIND=8) :: BF(402)
              REAL(KIND=8) :: EF(402)
              REAL(KIND=8) :: AK(201)
            END SUBROUTINE TWOSTRSW
          END INTERFACE 
        END MODULE TWOSTRSW__genmod
