        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 20 12:49:19 2014
        MODULE NEWFLUX1__genmod
          INTERFACE 
            SUBROUTINE NEWFLUX1(NLAYER,TAUL,PTEMPG,PTEMP,SLOPE,B3,GAMI, &
     &AK,GANGLE,GWEIGHT,Y3,EE1,CK1,CK2,RSFX,U1I,U1S,DIREC,DIRECTU)
              INTEGER(KIND=4) :: NLAYER
              REAL(KIND=4) :: TAUL(201)
              REAL(KIND=8) :: PTEMPG
              REAL(KIND=8) :: PTEMP(201)
              REAL(KIND=8) :: SLOPE(201)
              REAL(KIND=8) :: B3(201)
              REAL(KIND=8) :: GAMI(201)
              REAL(KIND=8) :: AK(201)
              REAL(KIND=8) :: GANGLE(3)
              REAL(KIND=8) :: GWEIGHT(3)
              REAL(KIND=8) :: Y3(3,201)
              REAL(KIND=8) :: EE1(201)
              REAL(KIND=8) :: CK1(201)
              REAL(KIND=8) :: CK2(201)
              REAL(KIND=4) :: RSFX
              REAL(KIND=8) :: U1I
              REAL(KIND=8) :: U1S
              REAL(KIND=8) :: DIREC(201)
              REAL(KIND=8) :: DIRECTU(201)
            END SUBROUTINE NEWFLUX1
          END INTERFACE 
        END MODULE NEWFLUX1__genmod
