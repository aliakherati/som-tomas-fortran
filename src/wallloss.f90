! function for k_turb implementation
! for internal loop
	subroutine wallloss(k_flat,ke,kw0,pres,T,Gc,Nk,Mk,dt,time,P_temp,P_num,WalVapMas,tau_w,kw_HYBRID) !normal case
!    subroutine wallloss(k_flat,ke,kw0,pres,T,Gc,Nk,Mk,dt,time,P_temp1,P_num,WalVapMas,tau_w,kw_HYBRID)  !for g-p-p on the wall system
	
	implicit none
	  
!-----ARGUMENT DECLARATIONS---------------------------------------------
	  
	integer, parameter :: iorg = 15
!       double precision T0M(1,1,1,446)  
	double precision dt, time
	double precision, parameter :: cdt = 1 !time step, seconds
	double precision T,pres
	double precision k_flat, kw0, ke ! from Jeff A model
	double precision cstar(iorg),cstarT(iorg) ! saturation ratios [ug m3]
!       data cstar /1d-3, 1d-2, 1d-1, 1d0, 1d1, 1d2, 1d3, 1d4/     ! for POA wall loss
	data cstar /1d-3,1d-2,1d-1, 1d0, 1d1, 1d2, 1d3,1d4,1d5, 1d6, 1d7, 1d8, 1d9, 1d10, 1d11/ 
   
	double precision Hvap(iorg) ! heat of vaporization [kJ mol^-1] 
!       data Hvap /230., 210., 189.3, 169.5, 150.8, 133.2, 116.6, 101.13/ !original 	 
!       data Hvap /97., 93.,89.,85.,81.,77.,73.,69./   !Andy's calculation
	data Hvap /97.,93.,89.,85.,81.,77.,73.,69.,65.,61.,57.,53.,49.,45.,41/ 

	double precision WalVapMas(iorg)
	double precision tau_w(iorg)	  
!       double precision diorg(iorg), ditemporg(iorg)
      ! Diffusion coefficients of the surrogate compounds in air at T_ref [m2/s]
      ! and the temperature-dependent factor
!       data diorg /iorg*5.0e-6/
!       data ditemporg /iorg*1.75/	  

!-----VARIABLE DECLARATIONS---------------------------------------------	  

	real (kind=8) rho, dil,rhoSI ! readin parameters
	real (kind=8), parameter :: pi = 3.1415926
	integer  i,j,k,n	!counters
	integer, parameter :: ibins = 36
	integer, parameter :: NBINS = 36
	integer, parameter :: IDTNUMD = 15
	integer, parameter :: icomp = 18
	!integer num
!       integer, parameter :: num = nint(600/cdt)
	integer, parameter :: num = 10

	real (kind=8) kB,g,R
	double precision xk(ibins+1),Mo,Dl,Dh,Dk(ibins),kw_turb(ibins)
	double precision kw_HYBRID(ibins)
	double precision Nk(ibins), Mk(ibins,icomp), Gc(icomp-1)
	real (kind=8) lambda, mu,  Cc(ibins),B(ibins),D(ibins),vT(ibins),x(ibins),D1(ibins)
	real (kind=8) V_bag, R_bag,StoV, S_bag, height, S1, S2
	real (kind=8) R_bag1,StoV1, S_bag1, factor ! for FSL volume test
	real (kind=8) Mw_p(iorg) ! mmolecular weight of organic species
	real (kind=8) gamma_p, alpha_w(iorg)
!       data alpha_w /1.80E-6, 1.16E-6,7.45E-7,4.79E-7,3.08E-7,1.98E-7,1.27E-7,8.17E-8/ ! Zhang et al. ACP, 2015
!       data alpha_w /1.0E-5,1.0E-5,1.0E-5,1.0E-5,1.0E-5,1.0E-5,1.0E-5,1.0E-5/  ! Ziemann et al., AST, 2008
!       data alpha_w /1.0E-0,1.0E-0,1.0E-0,1.0E-0,1.0E-0,1.0E-0,1.0E-0,1.0E-0/  
	data alpha_w /15*1.0E-5/
!       data alpha_w /15*1.0/
!       data alpha_w /1.80E-6,1.16E-6,7.45E-7,4.79E-7,3.08E-7,1.98E-7,1.27E-7,8.17E-8, &
!	& 5.25E-8,3.38E-8,2.17E-8,1.40E-8,8.97E-9,5.77E-9,3.71E-9/	  
	real (kind=8), external  :: debye1
	  
	real (kind=8) Na, Ra, MWair, morg(iorg), mair, dorg, dair
	real (kind=8) cmean(iorg),kw_gas(iorg), Dg(iorg)
!       real (kind=8), parameter :: wlc = 120.0
	real (kind=8) wlc(iorg)
!       data wlc /1.85E-6,1.30E-5,9.16E-5,6.44E-4,4.53E-3,3.18E-2,2.24E-1,1.57/  !Zhang et al., ACP, 2015, first round
!       data wlc /3.26E-6,2.29E-5,1.61E-4,1.13E-3,7.96E-3,5.60E-2,3.94E-1,2.77/  !Zhang et al., ACP, 2015, second round
!       data wlc /3.26E-6,2.29E-5,1.61E-4,1.13E-3,7.96E-3,5.60E-2,3.94E-1,2.77,1.95E1,1.37E2, &
!	& 9.62E2,6.76E3,4.75E4,3.34E5,2.35E6/
!       data wlc /1.30E-5,9.16E-5,6.44E-4,4.53E-3,3.18E-2,2.24E-1,1.57E0,11.1,7.78E1,5.47E2, &
!	& 3.85E3,2.70E4,1.90E5,1.34E6,9.40E6/ ! corrected Zhang data...

!       data wlc /0.064,0.064,0.064,0.064,0.255,1.014,4.038,16.076,80.,80., &
!	& 80.,80.,80.,80.,80./ ! Krechmer ES&T 2016	
	data wlc /15*120.0/	! in the unit of umol/m3
	real (kind=8) Gc_temp(iorg,num), P_temp(icomp),WalVap(iorg,num),Gci(iorg,num)
	real (kind=8) kw_gas_off(iorg),before(icomp),after(icomp)
	real (kind=8) P_temp1(ibins,icomp), before1(ibins,icomp) ! for g-p-p on the wall
	real (kind=8) P_num(ibins), before2(ibins)	  
	
	  
     !num = nint(dt/cdt)
!-----chamber properties------------------------------------------
	lambda = 65.2		! mean free path of particles in air,  nm
     mu = 1.458e-6*T**1.5/(110.4 + T) ! viscosity of air, kg/m-s	
     V_bag = 7. ! volume of filled bag, m3	
!     V_bag = 3014.	 
!      V_bag = 13.1e-3 !m3
!----for spherical chamber-----------------------------------------	 
!     R_bag = (V_bag*3./(4.*pi))**(1./3.) ! radius of bag, m
!     S_bag = 4.*pi*R_bag**2.	 
!-----for cubic chamber--------------------------------------------
      R_bag = 7.**(1./3.)
      S_bag = 6.*R_bag**2.
     StoV = S_bag/V_bag ! surface to volume ratio for smoke chamber
!------for FSL chamber --------------------------------
     ! R_bag1 = 7.**(1./3.)
     ! S_bag1 = 6.*R_bag1**2.
     ! StoV1 = S_bag1/7. ! surface to volume ratio for smoke chamber
     ! StoV = (12.4*12.4+12.4*19.6*2.)*2./V_bag	 !for FSL chamber
	  ! R_bag = (12.4*12.4*19.6)**(1./3.) !for FSL chamber
	  ! factor = StoV/StoV1
!-----for cylinder / rectangle reactor--------------------------------------------
!     height = V_bag/(14.*14.*1e-4) !m
!	 S1 = 14.0**2*1e-4 !m
!	 S2 = height*14.0*1e-2 !m
!     StoV = (S1*2.+S2*4.)/V_bag ! surface to volume ratio for smoke chamber  
	  
      factor = 1.0
	 
	 kB = 1.38e-23 ! Boltzmann constant, kg-m2/s2-K
     g = 9.81 ! gravity, m/s2
     R = 8.314 ! universal gas constant, kg-m3/s2-mol-K
	 rho = 1.4
	 rhoSI = 1000*rho !density of particles, kg/m3
	 gamma_p = 1.0
	 Na = 6.023e23 ! Avogadro's number, molec/mol


      do j=1,iorg
         cstarT(j)=cstar(j)*(298.0/T)*exp(Hvap(j)*1.d3*(1.d0/298.0-1.d0/T)/R)
      enddo	
  
	 
!  molecular weight of hypothetical compound, kg/mol, assumed from n-alkanes	 
	 do j=1,9
	  MW_p(j) = 0.434 - 0.045*log10(cstarT(j)) 
!      MW_p(j) = 0.200
	  morg(j) = MW_p(j)/Na	  
!	  print*, 'morg'
!	  print*, MW_p(j)*1000.0
	 end do
	 
    do j=10,15
	  MW_p(j) = 0.200
	  morg(j) = MW_p(j)/Na	  
	 end do     	
	 
	 
	 ! Size calculation
      Mo = 1.0d-21*2.d0**(-6.)

      do i=1,ibins+1
         xk(i)=Mo*2.d0**(i-1)		 
      end do
      
      do i=1,ibins
         Dl=1.0e+6*((6.*xk(i))/(1770.0*pi))**0.3333
         Dh=1.0e+6*((6.*xk(i+1))/(1770.0*pi))**0.3333
         Dk(i)=sqrt(Dl*Dh)*1000.
		 !print *,Dk(nn)
      end do
	  
	   ! calculated aerosol properties
	do i=1,ibins
	    kw_HYBRID(i) = 0.0
    end do
	 
	do i=1,ibins
		Cc(i) = 1. + (2*lambda/Dk(i))*(1.257 + 0.4*exp(-(0.55*Dk(i)/lambda)))
		B(i) = Cc(i)/(3*pi*mu*Dk(i)*1e-9) ! particle mobility, kg/s
		D(i) = kB*T*B(i) ! diffusivity, m2/s
		vT(i) = rhoSI*(Dk(i)*1e-9)**2*g/(18.*mu)
		x(i) = pi*vT(i)/(2.*sqrt(ke*D(i))) !  skip the fitting
		D1(i) = debye1(x(i)) 
 		kw_TURB(i) = 6.*sqrt(ke*D(i))/(pi*R_bag)*D1(i) + vT(i)/(4.*R_bag/3.) ! for cubic / spheric
!        kw_TURB(i) = 2.0*sqrt(ke*D(i))/pi*2*(S1+S2) + vT(i) / tanh(pi*vT(i)/4./sqrt(ke*D(i)))  ! for rectangle, cubic
        
		kw_HYBRID(i) = kw0*factor + kw_TURB(i)
!		write(*,*) 'kw_HYBRID'
!			print*, StoV1,StoV,kw0*factor*3600.,kw_TURB(i)*3600.,kw_HYBRID(i)*3600.
	end do


!	loss particle to the wall

	  !for bulky non-organic species, use k_flat
!   	  do n=1,IDTNUMD-1 !check the position

    
    do j=1,icomp
        before(j) = 0.0
        after(j) = 0.0        
    end do	

	  do n=1,ibins
	  	do j=1,icomp
	    before1(n,j) = 0.0
        end do
		before2(n) = 0.0
	  end do
	
!	do j=1,icomp+1
!   	    do n=1,NBINS !for size-resolved species, use k_turb
!          dilt2 = 1./kw_turb(n)
!      print *,kw_turb(n)*3600.
!          before(j) = before(j) + T0M(1,1,1,IDTNUMD-1+n+(j-1)*NBINS)
!	      T0M(1,1,1,IDTNUMD-1+n+(j-1)*NBINS)=  &
!     &                 T0M(1,1,1,IDTNUMD-1+n+(j-1)*NBINS)*  &
!     &              	 exp(-dt*kw_HYBRID(n))
!	     after(j) = after(j) + T0M(1,1,1,IDTNUMD-1+n+(j-1)*NBINS)
!   	    end do
!		P_temp(j) = P_temp(j)+before(j)-after(j)
!    end do

!!!!!!!!!!!!!!!!!!!!!! particle loss calculation, at 2014/12/11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
 	    do n=1,NBINS !for size-resolved species, use k_turb
		 before2(n) = Nk(n)
         Nk(n) = Nk(n)*exp(-dt*kw_HYBRID(n))
		 P_num(n) = P_num(n) + before2(n)- Nk(n)
		end do
		
		do j=1,icomp
          do n = 1, NBINS		! normal cases
!            if ((n .GE. 7) .AND. (n .LE. 22)) then ! for comparison with APE model
!             if (n .LE. 26) then !for coarse mode calculation
               before(j) = before(j) + Mk(n,j) !for normal case
!			   before1(n,j) = Mk(n,j)  !! for g-p-p on the wall
!		    endif
	      Mk(n,j)= Mk(n,j)* exp(-dt*kw_HYBRID(n))
!		    if ((n .GE. 7) .AND. (n .LE. 22)) then  ! for comparison with APE model
!            if (N .LE. 26) then ! for coarse mode calculation
		       after(j) = after(j) + Mk(n,j)
!			   P_temp1(n,j) = P_temp1(n,j) + before1(n,j) - Mk(n,j) ! for g-p-p on the wall
!		    endif
		  end do
	      P_temp(j) = P_temp(j)+before(j)-after(j)
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	

! Question: ke can be applied on both of the number and mass size distribution?
!   organic gas properties
!       alpha_w = 1.e-5 ! mass accomodation coefficient of organics onto chamber wall
	   !(i.e., sorption; Matsunaga and Ziemann, 2010)
       MWair = (0.79*28. + 0.21*32.)/1000.  ! molecular weight of air, kg/mol
       
       mair = MWair/Na;
        dorg = 10e-10  ! diameter of organic molecule, m
       dair = (0.79*1.09e-10 + 0.21*1.21e-10) ! diameter of air molecule, m
	   

	do j=1,iorg
	   cmean(j) = sqrt(3.*R*T/MW_p(j)) ! root mean square speed of gas, m/s
	   Dg(j) = (2./3.)*sqrt(kB*T/pi*0.5*(1./morg(j)+1/mair))*1./ &
	       (pi*(0.5*(dorg + dair))**2.)/Na*(R*T/pres) ! diffusivity of organic vapors, m2/s
!       Dg(j) = diorg(j)*(T/298.0)**ditemporg(j)
       kw_gas(j) = StoV*(alpha_w(j)*cmean(j)/4.)/(1.+(pi/2.)*(alpha_w(j)*cmean(j)/(4.*sqrt(ke*Dg(j))))) ! Zhang et al., 2014
!       solubility(j,m,time) = (cmean(j,m)*sqrt(pi*time)/(4*10^logHb(j,m)*Ra*T(j)*sqrt(Dg(j,m))))^-1; % inverse uptake based on Henry's law
!        gamma(j) = (1./alpha_w + cmean(j)*sqrt(pi*time)/(4.*10.^logHb(j)*Ra*T*sqrt(Dg(j))))^-1. ! uptake accounting for mass accommodation and solubility, uptake coeffcient
!        kw_gas(j) = StoV*(gamma(j)*cmean(j)/4.)/(1.+(pi/2.)*(gamma(j)*cmean(j)/(4.*sqrt(ke*Dg(j))))) ! gas-phase wall loss rate, 1/s
!       kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc/200./gamma_p
!       kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc/(MW_p(j)*1000.0)/gamma_p
		kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc(j)/(MW_p(j)*1000.0)/gamma_p
	   tau_w(j) = 1./(kw_gas(j)+kw_gas_off(j))
! 	  T0M(1,1,1,j+6)= T0M(1,1,1,j+6)*exp(-kw_gas(j)*t)
!     write(*,*) 'vapor loss rate constant'
!     print *,'kw_gas=',kw_gas(j),kw_gas_off(j) !,Dg(j)	 
     end do	
	 
! swap T0M array into Gc_temp array dependent on time step
     do j=1,iorg
	   Gc_temp(j,1) =  Gc(j+1)
!	   write(*,*) 'initial vapor concentration'
!	   print*, Gc_temp(j,1)/7000000.*1.E15
!	   Gc_temp2(j,1) = 0.0
	   WalVap(j,1) = WalVapMas(j)
	   Gci(j,1) = 0.0
	 end do
	 
!	 ! Euler's method
!      do k=2,num   ! time step
!         do j=1,iorg !species
!		 if ((WalVap(j,k-1)+Gc_temp(j,k-1)*kw_gas(j)-WalVap(j,k-1)*kw_gas_off(j)) .GE. 0.0) then		 ! added in vbs extension
!		  Gc_temp(j,k) = Gc_temp(j,k-1)+WalVap(j,k-1)*kw_gas_off(j)-Gc_temp(j,k-1)*kw_gas(j)
!          WalVap(j,k) = WalVap(j,k-1)+Gc_temp(j,k-1)*kw_gas(j)-WalVap(j,k-1)*kw_gas_off(j)
!		  else
!         Gc_temp(j,k) = Gc_temp(j,k-1) + WalVap(j,k-1)
!		 WalVap(j,k) = 0.0
!		 endif
!         end do
!      end do		

     ! integrate method
      do k=2,num   !time step
          do j=1,iorg	  !species
		  Gc_temp(j,k) = Gc_temp(j,k-1)*exp(-kw_gas(j)*cdt)+WalVap(j,k-1)*(1-exp(-kw_gas_off(j)*cdt))
		  WalVap(j,k) = WalVap(j,k-1)*exp(-kw_gas_off(j)*cdt)+Gc_temp(j,k-1)*(1-exp(-kw_gas(j)*cdt))
		  end do
	  end do
 

!       print *, num,WalVap(8,1),WalVap(8,num),WalVap(8,1)/WalVap(8,num)
!   Swap Gc_temp back into T0M array and Swap out vapour wall loss
      do j=1,iorg
	    Gc(j+1) = Gc_temp(j,num)
		if (WalVap(j,num) .LT. 1e-30) then  !add in vbs extension
		   WalVapMas(j) = 0.0
		else
	    WalVapMas(j) = WalVap(j,num)
		endif
	  end do
	     
	
100 Format(1X,E9.3)	
	END subroutine
	
 
	
!debye1 function, taken from MISFUN	
	
	function debye1 ( xvalue )

!*****************************************************************************80
!
!! DEBYE1 calculates the Debye function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE1, the value of the function.
!
  implicit none

  real ( kind = 8 ) adeb1(0:18)
  real ( kind = 8 ) cheval
  real ( kind = 8 ) debye1
  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nexp
  integer ( kind = 4 ), parameter :: nterms = 15
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: quart = 0.25D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) debinf,expmx, &
       nine,rk,sum1,t,thirt6,xk,xlim,xlow, &
       xupper
  data nine,thirt6 /9.0d0, 36.0d0 /
  data debinf/0.60792710185402662866d0/
  data adeb1/2.40065971903814101941d0, &
             0.19372130421893600885d0, &
            -0.623291245548957703d-2, &
             0.35111747702064800d-3, &
            -0.2282224667012310d-4, &
             0.158054678750300d-5, &
            -0.11353781970719d-6, &
             0.835833611875d-8, &
            -0.62644247872d-9, &
             0.4760334890d-10, &
            -0.365741540d-11, &
             0.28354310d-12, &
            -0.2214729d-13, &
             0.174092d-14, &
            -0.13759d-15, &
             0.1093d-16, &
            -0.87d-18, &
             0.7d-19, &
            -0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
  data xlow,xupper,xlim/0.298023d-7,35.35051d0,708.39642d0/

  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE1 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye1 = zero
  else if ( x < xlow ) then
    debye1 = ( ( x - nine ) * x + thirt6 ) / thirt6
  else if ( x <= four ) then
    t = ( ( x * x / eight ) - half ) - half
    debye1 = cheval ( nterms, adeb1, t ) - quart * x
  else

    debye1 = one / ( x * debinf )
    if ( x < xlim ) then
      expmx = exp ( -x )
      if ( xupper < x ) then
        debye1 = debye1 - expmx * ( one + one / x )
      else
        sum1 = zero
        rk = aint ( xlim / x )
        nexp = int ( rk )
        xk = rk * x
        do i = nexp, 1, -1
          t = ( one + one / xk ) / rk
          sum1 = sum1 * expmx + t
          rk = rk - one
          xk = xk - x
        end do
        debye1 = debye1 - sum1 * expmx
      end if
    end if
  end if

  return
end
subroutine debye1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE1_VALUES returns some values of Debye's function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
            0.99951182471380889183D+00, &
            0.99221462647120597836D+00, &
            0.96918395997895308324D+00, &
            0.88192715679060552968D+00, &
            0.77750463411224827642D+00, &
            0.68614531078940204342D+00, &
            0.60694728460981007205D+00, &
            0.53878956907785587703D+00, &
            0.48043521957304283829D+00, &
            0.38814802129793784501D+00, &
            0.36930802829242526815D+00, &
            0.32087619770014612104D+00, &
            0.29423996623154246701D+00, &
            0.27126046678502189985D+00, &
            0.20523930310221503723D+00, &
            0.16444346567994602563D+00, &
            0.10966194482735821276D+00, &
            0.82246701178200016086D-01, &
            0.54831135561510852445D-01, &
            0.32898681336964528729D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
             0.0019531250D+00, &
             0.0312500000D+00, &
             0.1250000000D+00, &
             0.5000000000D+00, &
             1.0000000000D+00, &
             1.5000000000D+00, &
             2.0000000000D+00, &
             2.5000000000D+00, &
             3.0000000000D+00, &
             4.0000000000D+00, &
             4.2500000000D+00, &
             5.0000000000D+00, &
             5.5000000000D+00, &
             6.0000000000D+00, &
             8.0000000000D+00, &
            10.0000000000D+00, &
            15.0000000000D+00, &
            20.0000000000D+00, &
            30.0000000000D+00, &
            50.0000000000D+00 /)
!
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x  = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return

end

function cheval ( n, a, t )

!*****************************************************************************80
!
!! CHEVAL evaluates a Chebyshev series.
!
!  Discussion:
!
!    This function evaluates a Chebyshev series, using the
!    Clenshaw method with Reinsch modification, as analysed
!    in the paper by Oliver.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    J Oliver,
!    An error analysis of the modified Clenshaw method for
!    evaluating Chebyshev and Fourier series,
!    Journal of the IMA,
!    Volume 20, 1977, pages 379-391.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of terms in the sequence.
!
!    Input, real ( kind = 8 ) A(0:N), the coefficients of the Chebyshev series.
!
!    Input, real ( kind = 8 ) T, the value at which the series is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CHEVAL, the value of the Chebyshev series at T.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: a(0:n)
  real ( kind = 8 ) :: cheval
  real ( kind = 8 ) :: d1
  real ( kind = 8 ) :: d2
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: t
  real ( kind = 8 ), parameter :: test = 0.6D+00
  real ( kind = 8 ) :: tt
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) :: u0
  real ( kind = 8 ) :: u1
  real ( kind = 8 ) :: u2
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  u1 = zero
!
!  T <= -0.6, Reinsch modification.
!
  if ( t <= -test ) then

    d1 = zero
    tt = ( t + half ) + half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) - d2
      u1 = d1 - u2
    end do

    cheval = ( d1 - d2 ) / two
!
!  -0.6 < T < 0.6, Standard Clenshaw method.
!
  else if ( t < test ) then

    u0 = zero
    tt = t + t

    do i = n, 0, -1
      u2 = u1
      u1 = u0
      u0 = tt * u1 + a(i) - u2
    end do

    cheval = ( u0 - u2 ) / two
!
!  0.6 <= T, Reinsch modification.
!
  else

    d1 = zero
    tt = ( t - half ) - half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) + d2
      u1 = d1 + u2
    end do

    cheval = ( d1 + d2 ) / two

  end if

  return
end
	  
