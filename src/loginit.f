
C     **************************************************
C     *  loginit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Initializes aerosol number and mass distributions with sulfate
C     only and a lognormal size distribution.

C-----INPUTS------------------------------------------------------------

C     No - total number concentration (#/cm3)
C     Dp - lognormal mean diameter (microns)
C     sigma - width parameter

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE loginit(N1, Dp, sigma, idtcomp)

C     INCLUDE FILES

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision N1, Dp, sigma
      integer idtcomp    !indicates chemical component to be added

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer b,n     !bin and species counters
      double precision Dl, Dh     !lower and upper bounds of bin (microns)
      double precision Dk         !mean diameter of current bin (microns)
      double precision np         !number of particles in size bin in GCM cell
      double precision pi
      double precision No
      
C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

C Convert No from #/cm3 to #/box
      No=N1*boxvol

C Loop over grid cells
CAliA      do l=1,lm
CAliA      do j=1,jm
CAliA      do i=1,im

      !Loop over number of size bins
      DO N=1,IBINS

         !Calculate diameter of this size bin
CAliA         Dl=1.0e+6*((6.*xk(n))/(1770.0*3.14))**0.3333
CAliA         Dh=1.0e+6*((6.*xk(n+1))/(1770.0*3.14))**0.3333
         Dl=1.0e+6*((6.*xk(n))/(densorg(1)*3.14))**0.3333
         Dh=1.0e+6*((6.*xk(n+1))/(densorg(1)*3.14))**0.3333
         Dk=sqrt(Dl*Dh)
c         write(*,*) 'Dk(',n,')= ',Dk,' microns'

         !Calculate number concentration
         np=No/(sqrt(2.*pi)*Dk*log(sigma))*exp(-( (log(Dk/Dp))**2./
     &         (2.*(log(sigma))**2.) )) * (Dh-Dl)

         Nk(n) = Nk(n) + np
CAliA         T0M(i,j,l,IDTNUMD+n-1)=T0M(i,j,l,IDTNUMD+n-1)+np
c         print*,'N=',Nk(n)
c         write(*,*)'N=',T0M(i,j,l,IDTNUMD+n-1)

         !Calculate component mass
CCC      AliA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
         Mk(n,idtcomp) = Mk(n,idtcomp) + np*sqrt(xk(n))*sqrt(xk(n+1))
CAliA         Mk(n,strso4) = Mk(n,strso4) + np*sqrt(xk(n))*sqrt(xk(n+1))
CAliA         T0M(i,j,l,idtcomp+n-1)=T0M(i,j,l,idtcomp+n-1)
CAliA     &                         +np*sqrt(xk(n))*sqrt(xk(n+1))
CCC      ////////////////////////////////////////////////////////////////// AliA
c         write(*,*)'M=',T0M(i,j,l,IDTSO4+n-1)

      END DO

      RETURN
      END
