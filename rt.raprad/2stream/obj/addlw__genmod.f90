        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 12:49:17 2014
        MODULE ADDLW__genmod
          INTERFACE 
            SUBROUTINE ADDLW(NLAYER,TAUL,W0,G0,RSFX,OPD,OPDND,AK,B1,B2, &
     &B3,EM1,EM2,EL1,EL2,AF,BF,EF,U0,SLOPE,PTEMPG,PTEMP,U1I,U1S,FNET,   &
     &SOL_FLUXUP,SOL_FLUXDN,DIRECT_ND,IRFLAG,CK1,CK2)
              INTEGER(KIND=4) :: NLAYER
              REAL(KIND=4) :: TAUL(*)
              REAL(KIND=4) :: W0(*)
              REAL(KIND=4) :: G0(*)
              REAL(KIND=4) :: RSFX
              REAL(KIND=4) :: OPD(201)
              REAL(KIND=4) :: OPDND(201)
              REAL(KIND=8) :: AK(201)
              REAL(KIND=8) :: B1(201)
              REAL(KIND=8) :: B2(201)
              REAL(KIND=8) :: B3(201)
              REAL(KIND=8) :: EM1(201)
              REAL(KIND=8) :: EM2(201)
              REAL(KIND=8) :: EL1(201)
              REAL(KIND=8) :: EL2(201)
              REAL(KIND=8) :: AF(402)
              REAL(KIND=8) :: BF(402)
              REAL(KIND=8) :: EF(402)
              REAL(KIND=4) :: U0
              REAL(KIND=8) :: SLOPE(201)
              REAL(KIND=8) :: PTEMPG
              REAL(KIND=8) :: PTEMP(201)
              REAL(KIND=8) :: U1I
              REAL(KIND=8) :: U1S
              REAL(KIND=8) :: FNET(201)
              REAL(KIND=8) :: SOL_FLUXUP(201)
              REAL(KIND=8) :: SOL_FLUXDN(201)
              REAL(KIND=8) :: DIRECT_ND(201)
              INTEGER(KIND=4) :: IRFLAG
              REAL(KIND=8) :: CK1(201)
              REAL(KIND=8) :: CK2(201)
            END SUBROUTINE ADDLW
          END INTERFACE 
        END MODULE ADDLW__genmod