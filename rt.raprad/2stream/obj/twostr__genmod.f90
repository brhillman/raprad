        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 18 18:42:54 2014
        MODULE TWOSTR__genmod
          INTERFACE 
            SUBROUTINE TWOSTR(NLAYER,IRFLAG,TAUL,W0,G0,RSFX,B1,B2,EL1,  &
     &EL2,EM1,EM2,AF,BF,EF,AK,U1I,U1S,GAMI,EE1)
              INTEGER(KIND=4) :: NLAYER
              INTEGER(KIND=4) :: IRFLAG
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
              REAL(KIND=8) :: U1I
              REAL(KIND=8) :: U1S
              REAL(KIND=8) :: GAMI(201)
              REAL(KIND=8) :: EE1(201)
            END SUBROUTINE TWOSTR
          END INTERFACE 
        END MODULE TWOSTR__genmod
