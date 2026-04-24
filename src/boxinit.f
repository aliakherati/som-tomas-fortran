
C     **************************************************
C     *  boxinit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, July 2000

C     Initializes the box model for testing size-resolved microphysics.
C     Starts either with zero concentration or with a restart file if
C     one exists.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE boxinit(time,No,Dp,sigma,INDSP)

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision time
      double precision No
      double precision Dp
      double precision sigma
      INTEGER INDSP

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,l,n

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

CAliA 2    format(200E15.6)

C-----CODE--------------------------------------------------------------

C Initializations

      !Initialize everything to zero at first
      DO N=1,IBINS
         Nk(N) = 0.0
         DO J=1,IORG
            Mk(N,J) = 0.0 
         END DO
      END DO

CAliA      call speciesmap()
      call initbounds()
CAliA      call loginit(No,Dp,sigma,strso4)
      call loginit(No,Dp,sigma,INDSP)
CAliA      call loginit(No,Dp,sigma,IDTSO4)
C       call loginit_2mode(No,Dp,sigma,IDTSO4) !for 2 mode input

C initallize to ammonium sulfate
CCC      DO N=1,IBINS
c         Mk(N,strnh4) = Mk(N,strnh4)/96.0*18.0*2.0
CCC         Mk(N,strnh4) = 0.0
CC      END DO
      
c      do i=1,IM
c         do j=1,JM
c            do l=1,LM
c               T0M(i,j,l,IDTNUMD+3)=10000.d0*boxvol
c               T0M(i,j,l,IDTSO4+3)=10000.d0*boxvol*sqrt(xk(5)*xk(4))
c            enddo
c         enddo
c      enddo   



C Read initial conditions from box.rsf if it exists
CAliA      open(unit=11,file='box.rsf',status='old',err=100)
CAliA      read(11,2) time,(T0M(1,1,1,n),n=1,NTM)
CAliA      close(11)
CAliA      time=time*3600.
CAliA      goto 200

C Otherwise, use these simple initial conditions
CAliA 100  continue
CAliA      T0M(1,1,1,IDTH2SO4)=0.
CAliA      HNO3(1,1,1,1)=100.0e-12    !mixing ratio
CAliA 200  continue

      RETURN
      END
