
      SUBROUTINE DIFF(NEQN, C, A, B, F, TIN, DT, CONVG)
C
C     Program written by Armistead Russell
C     Revised version by Robert Harley (June 1992)
C     Revised for SAPRC mechanism compiler by Robert Harley (August 1992)
C     Revised to add elevated sources and vert adv by Mike Kleeman (March 2000)
C     Revised 2026-04: interface promoted to REAL*8 to match the upgraded
C                      INTEGR2 driver. The compiled SAPRC mechanism
C                      (DIFUN in saprc14_rev1.f) is auto-generated as
C                      REAL*4, so we bridge through REAL*4 scratch buffers
C                      here. Truncation noise on A, B is bounded by
C                      ~1 ULP of REAL*4 (~1e-7) per call.
C
C...PROGRAM VARIABLES...
C
C     A      - A COEFFICIENT FOR THE DIFFERENTIAL EQUATIONS
C     B      - B COEFFICIENT FOR THE DIFFERENTIAL EQUATIONS
C     C      - CONCENTRATION ARRAY AT TIME T
C     F      - DC/DT AT TIME T
C     TIN    - CURRENT TIME
C     DT     - CURRENT TIME STEP
C     CONVG  - LOGICAL FLAG TRUE IF LAST STEP CONVERGED
C
C...SUBPROGRAMS REQUIRED...
C
C     DIFUN - PROGRAM TO CALCULATE THE CHEMICAL REACTION RATES
C
      IMPLICIT NONE

      INCLUDE 'modlspc.h'
      INCLUDE 'gaskin.h'
      INCLUDE 'common.inc'

      INTEGER I, J
      INTEGER NEQN

      REAL*8  TIN, DT
      REAL*8  C(NEQN), A(NEQN), B(NEQN), F(NEQN)

C     REAL*4 scratch buffers used to interface with the auto-generated
C     REAL*4 DIFUN mechanism callback.
      REAL    C_R4(NEQN), A_R4(NEQN), B_R4(NEQN)
      REAL    RXNRAT(maxrxn)

      logical convg
C
C     COMMON BLOCKS
C
CAliA      COMMON / PARAM / NSPECS, CONST(ncz)
CAliA      COMMON / TRAJT / EMT(LCOLSP)
C
C----------------------------------------------------------------------
C
C     STEP 1:- INITIALIZE THE DIFFERENTIAL TERM
C
C----------------------------------------------------------------------
C
      DO 10 J=1, NEQN
        A_R4(J) = 0.0
        B_R4(J) = 0.0
        C_R4(J) = REAL(C(J))
 10   CONTINUE

C----------------------------------------------------------------------
C
C     STEP 2:- CALCULATE THE CONTRIBUTION DUE TO REACTION
C
C----------------------------------------------------------------------
C
      CALL DIFUN(CONST, C_R4, PSSA, RKZ, RXNRAT,
     &     COEF, A_R4, B_R4)

C--------------------------------------------------------------------
C
C     PHASE 3 :- CALCULATE THE CONTRIBUTION FROM DIRECT SOURCES
C
C---------------------------------------------------------------------
C
CAliA      DO J=1, NEQN
CAliA         A(J) = A(J) + EMT(J)
CAliA      ENDDO
C
C-----------------------------------------------------------------------
C
C     STEP 3 :- promote A, B back to REAL*8 and compute net rate
C
C-----------------------------------------------------------------------
C
      DO 200 I = 1, NEQN
        A(I) = DBLE(A_R4(I))
        B(I) = DBLE(B_R4(I))
        F(I) = A(I) - B(I)*C(I)
 200  CONTINUE
C
C-----END OF THE PROGRAM----
C
      RETURN
      END

