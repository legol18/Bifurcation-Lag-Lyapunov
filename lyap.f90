   MODULE CLIMATE
   implicit none
   integer, parameter:: id=7, mlyap=2, NEQ=(mlyap+1)*id, idx1=2, idx2=3, idx3=4, idx4=5, idx5=6!! id is the system dimension, idx1, idx2 are signals to compare for lag
   integer, parameter:: npar=17, ntrans=100000, niter=100000, ntot=ntrans+niter, nhist=50000, imaxlag=niter/10 !! npar is number of system parameters, ntrans, niter are transient and iteration loop sizes
   real(kind=8), parameter:: h=0.1d0  !! time step for iteration
   real(kind=8):: par(npar) !! parameter array
   
  CONTAINS
  
     SUBROUTINE similarity(xx,yy,simfx)!! similarity function calculation, DO NOT TOUCH unless you REALLY know what you are doing!
     implicit none
     real(kind=8), intent(in):: xx(niter),yy(niter)
     real(kind=8), intent(out):: simfx(imaxlag)
     real(kind=8):: simsum, sqm1,sqm2
     integer:: i,j
     
     simfx=0.d0

     do i=1,imaxlag
       simsum=0.d0; sqm1=0.d0; sqm2=0.d0
       do j=1,niter-i
        simsum=simsum+(xx(j)-yy(j+i))**2
        sqm1=sqm1+xx(j)**2/dfloat(niter-i)
        sqm2=sqm2+yy(j)**2/dfloat(niter-i)
       enddo
      simfx(i)=(simsum/dfloat(niter-i))/sqrt(sqm1*sqm2)
     enddo
     
     RETURN
     end SUBROUTINE similarity

        SUBROUTINE DERIVS(NEQ,T,Y,YDOT) !! system definition below. Only to change the system if need be.
        IMPLICIT NONE
        INTEGER, intent(in):: NEQ
        real(kind=8), intent(in):: T, Y(NEQ) 
        real(kind=8), intent(out):: YDOT(NEQ)
        real(kind=8):: jac(id,id),jsum
        real(kind=8):: D,S0,gamma,n,rpb,rbb,rmb,Hpb,Hbb,Hmb,ae,gpp,gbp,Kpp,Kbp,epp,ebp,w,V,c!!!parameters
        integer:: itr,j,ix,iy

    ! Parameterset
        ! diluation rate
        D = par(1); 
        ! nutrient availability
        S0 = par(2)!; %0.08; %3
        ! Physiological prey parameters
        !gamma = par(4);! %proportion of C:N->4.5 in marine bacteria
        n = par(3); !% n<1 increase grazing preasure on G. ecological cost of phenot. plat. (n=1 linear trade off, n=0.6 strong trade off, n=0.2 very strong trade off)e.g. increased overall grazing preasure n<1
        rpb = par(4)
        rbb = par(5)
        rmb = par(6)
        Hpb = par(7)
        Hbb = par(8)
        Hmb = par(9)
        ae = par(10)  !! assimilation efficieny
        ! Physiological predator parameters
        gpp = par(11);! %1.2
        gbp = par(12);
        Kpp = par(13);! %0.1
        Kbp = par(14);
        !epp = par(15);! %0.3
        !ebp = par(16);! %0.3
        
        ! Generalist parameters
        !%c = 0.01; %1 in the paper Gaedke et al.% reduction of max growth compared with the specialist (direct costs)e.g. 0.2 means 20% less groth rate than specialist
        c = par(15);! %1 in the paper Gaedke et al.% reduction of max growth compared with the specialist (direct costs)e.g. 0.2 means 20% less groth rate than specialist
        w = par(16);! % Scaling of non adaptative trait adjustment
        V = par(17) ! (speed of adaptation (0.01 slow, 0.1 intemediate, 1 fast trait changes)
        
		YDOT(1)=D*(S0-y(1))-(rpb*y(1)/(y(1)+Hpb)*y(2)+rbb*y(1)/(y(1)+Hbb)*y(3)+rmb*(1-c)*y(1)/(y(1)+Hmb)*y(4))  							!!Nutrients															!! Nutrients

        YDOT(2)=(ae*rpb*(y(1)/(y(1)+Hpb))-D)*y(2)-gpp*(y(2)/(y(2)+(1-y(7))**n*y(4)+Kpp))*y(5)                                          !! Planktonic bacteria

        YDOT(3)=(ae*rbb*(y(1)/(y(1)+Hbb))-D)*y(3)-gbp*(y(3)/(y(3)+(y(7))**n*y(4)+Kbp))*y(6)            									!! Biofilm bacteria

        YDOT(4)=(ae*rmb*(1-c)*y(1)/(y(1)+Hmb)-D)*y(4)-gpp*(1-y(7))**n*y(4)*y(5)/(y(2)+(1-y(7))**n*y(4)+Kpp)-gbp*(y(7)**n*y(4)) & 
                *y(6)/(y(3)+(y(7)**n*y(4))+Kbp)                                                                                        !! Generalist bacteria

        YDOT(5)=(ae*gpp*(y(2)+(1-y(7))**n*y(4))/(y(2)+(1-y(7))**n*y(4)+Kpp)-D)*y(5)                                                     !! Planktonic Predator

        YDOT(6)=(ae*gbp*(y(3)+y(7)**n*y(4))/(y(3)+y(7)**n*y(4)+Kbp)-D)*y(6) 															!! Biofilm Predator
		
		YDOT(7)=V*(-gpp*y(5)*(-n*(1-y(7))**(n-1)*(y(2)+Kpp)/(y(2)+(1-y(7))**n*y(4)+Kpp)**2)-gbp*y(6)*(n*y(7)**(n-1)	&
				*(y(3)+Kbp)/(y(3)+y(7)**n*y(4)+Kbp)**2))+(w/(y(7))**2)-(w/(1-y(7))**2) 													!Individual change in alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JACOBIAN   
        jac(1,1)=(rbb*y(1)*y(3))/(Hbb + y(1))**2 - (rbb*y(3))/(Hbb + y(1)) - &
        &(rpb*y(2))/(Hpb + y(1)) - D + (rpb*y(1)*y(2))/(Hpb + y(1))**2 + &
        &(rmb*y(4)*(c-1.d0))/(Hmb + y(1)) - (rmb*y(1)*y(4)*(c-1.d0))/(Hmb + y(1))**2

        jac(1,2)=-(rpb*y(1))/(Hpb + y(1))

        jac(1,3)=-(rbb*y(1))/(Hbb + y(1))

        jac(1,4)=(rmb*y(1)*(c-1.d0))/(Hmb + y(1))

        jac(1,5)=0.d0

        jac(1,6)=0.d0

        jac(1,7)=0.d0

        jac(2,1)=y(2)*((ae*rpb)/(Hpb + y(1)) - (ae*rpb*y(1))/(Hpb + y(1))**2)

        jac(2,2)=(gpp*y(2)*y(5))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2 - &
        &(gpp*y(5))/(Kpp + y(2) + y(4)*(1 - y(7))**n) - D + (ae*rpb*y(1))/(Hpb + y(1))

        jac(2,3)=0.d0

        jac(2,4)=(gpp*y(2)*y(5)*(1.d0 - y(7))**n)/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n)**2

        jac(2,5)=-(gpp*y(2))/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n)

        jac(2,6)=0.d0

        jac(2,7)=-(gpp*n*y(2)*y(4)*y(5)*(1.d0-y(7))**(n-1))/(Kpp + y(2) + y(4)*(1.d0-y(7))**n)**2

        jac(3,1)=y(3)*((ae*rbb)/(Hbb + y(1)) - (ae*rbb*y(1))/(Hbb + y(1))**2)

        jac(3,2)=0.d0

        jac(3,3)=(ae*rbb*y(1))/(Hbb + y(1)) - (gbp*y(6))/(Kbp + y(3) + &
        &y(4)*y(7)**n) - D + (gbp*y(3)*y(6))/(Kbp + y(3) + y(4)*y(7)**n)**2

        jac(3,4)=(gbp*y(3)*y(6)*y(7)**n)/(Kbp + y(3) + y(4)*y(7)**n)**2

        jac(3,5)=0.d0

        jac(3,6)=-(gbp*y(3))/(Kbp + y(3) + y(4)*y(7)**n)

        jac(3,7)=(gbp*n*y(3)*y(4)*y(6)*y(7)**(n - 1))/(Kbp + y(3) + y(4)*y(7)**n)**2

        jac(4,1)=-y(4)*((ae*rmb*(c-1.d0))/(Hmb + y(1)) - (ae*rmb*y(1)*(c-1.d0))/(Hmb + y(1))**2)

        jac(4,2)=(gpp*y(4)*y(5)*(1 - y(7))**n)/(Kpp + y(2) + y(4)*(1.d0-y(7))**n)**2

        jac(4,3)=(gbp*y(4)*y(6)*y(7)**n)/(Kbp + y(3) + y(4)*y(7)**n)**2

        jac(4,4)=(gbp*y(4)*y(6)*y(7)**(2*n))/(Kbp + y(3) + y(4)*y(7)**n)**2 - &
        &(gbp*y(6)*y(7)**n)/(Kbp + y(3) + y(4)*y(7)**n) - &
        &(gpp*y(5)*(1.d0-y(7))**n)/(Kpp + y(2) + y(4)*(1.d0-y(7))**n) - D + &
        &(gpp*y(4)*y(5)*(1.d0-y(7))**(2*n))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2 - &
        &(ae*rmb*y(1)*(c-1.d0))/(Hmb + y(1))

        jac(4,5)=-(gpp*y(4)*(1.d0 - y(7))**n)/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n)

        jac(4,6)=-(gbp*y(4)*y(7)**n)/(Kbp + y(3) + y(4)*y(7)**n)

        jac(4,7)=(gpp*n*y(4)*y(5)*(1.d0-y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1 - y(7))**n) &
        &-(gbp*n*y(4)*y(6)*y(7)**(n - 1))/(Kbp + y(3) + y(4)*y(7)**n) +&
        &(gbp*n*y(4)**2*y(6)*y(7)**n*y(7)**(n - 1))/(Kbp + y(3) + y(4)*y(7)**n)**2 -&
        &(gpp*n*y(4)**2*y(5)*(1 - y(7))**n*(1 - y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2

        jac(5,1)=0.d0

        jac(5,2)=y(5)*((ae*gpp)/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n) -&
        &(ae*gpp*(y(2) + y(4)*(1.d0 - y(7))**n))/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n)**2)

        jac(5,3)=0.d0

        jac(5,4)=y(5)*((ae*gpp*(1.d0 - y(7))**n)/(Kpp + y(2) + y(4)*(1 - y(7))**n) -&
        &(ae*gpp*(y(2) + y(4)*(1.d0 - y(7))**n)*(1.d0 - y(7))**n)/(Kpp + y(2) +&
        &y(4)*(1.d0 - y(7))**n)**2)

        jac(5,5)=(ae*gpp*(y(2) + y(4)*(1.d0 - y(7))**n))/(Kpp + y(2) + y(4)*(1.d0 - y(7))**n) - D

        jac(5,6)=0.d0

        jac(5,7)=-y(5)*((ae*gpp*n*y(4)*(1.d0 - y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1 - y(7))**n) -&
        &(ae*gpp*n*y(4)*(y(2) + y(4)*(1.d0 - y(7))**n)*(1.d0 - y(7))**(n - 1))/(Kpp + y(2) +&
        &y(4)*(1 - y(7))**n)**2)

        jac(6,1)=0.d0

        jac(6,2)=0.d0

        jac(6,3)=y(6)*((ae*gbp)/(Kbp + y(3) + y(4)*y(7)**n) -&
        &(ae*gbp*(y(3) + y(4)*y(7)**n))/(Kbp + y(3) + y(4)*y(7)**n)**2)

        jac(6,4)=y(6)*((ae*gbp*y(7)**n)/(Kbp + y(3) + y(4)*y(7)**n) - &
        &(ae*gbp*y(7)**n*(y(3) + y(4)*y(7)**n))/(Kbp + y(3) + y(4)*y(7)**n)**2)

        jac(6,5)=0.d0

        jac(6,6)=(ae*gbp*(y(3) + y(4)*y(7)**n))/(Kbp + y(3) + y(4)*y(7)**n) - D

        jac(6,7)=y(6)*((ae*gbp*n*y(4)*y(7)**(n - 1))/(Kbp + y(3) + y(4)*y(7)**n) -&
        &(ae*gbp*n*y(4)*y(7)**(n - 1)*(y(3) + y(4)*y(7)**n))/(Kbp + y(3) + y(4)*y(7)**n)**2)

        jac(7,1)=0.d0

        jac(7,2)=V*((gpp*n*y(5)*(1 - y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2 - &
        &(2*gpp*n*y(5)*(Kpp + y(2))*(1 - y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1.d0-y(7))**n)**3)

        jac(7,3)=-V*((gbp*n*y(6)*y(7)**(n - 1))/(Kbp + y(3) + y(4)*y(7)**n)**2 - &
        &(2*gbp*n*y(6)*y(7)**(n - 1)*(Kbp + y(3)))/(Kbp + y(3) + y(4)*y(7)**n)**3)

        jac(7,4)=V*((2*gbp*n*y(6)*y(7)**n*y(7)**(n - 1)*(Kbp + y(3)))/(Kbp + y(3) + y(4)*y(7)**n)**3 &
        &-(2*gpp*n*y(5)*(Kpp + y(2))*(1 - y(7))**n*(1 - y(7))**(n - 1))/(Kpp + y(2) + &
        &y(4)*(1 - y(7))**n)**3)

        jac(7,5)=(V*gpp*n*(Kpp + y(2))*(1 - y(7))**(n - 1))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2

        jac(7,6)=-(V*gbp*n*y(7)**(n - 1)*(Kbp + y(3)))/(Kbp + y(3) + y(4)*y(7)**n)**2

        jac(7,7)=(2*w)/(y(7) - 1)**3 + V*((2*gbp*n**2*y(4)*y(6)*y(7)**(2*n - 2)*(Kbp & 
        &+y(3)))/(Kbp + y(3) + y(4)*y(7)**n)**3 - (gpp*n*y(5)*(Kpp + y(2))*(n - 1)*(1 - &
        &y(7))**(n - 2))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**2 + (2*gpp*n**2*y(4)*y(5)*(Kpp + &
        &y(2))*(1 - y(7))**(2*n - 2))/(Kpp + y(2) + y(4)*(1 - y(7))**n)**3 - &
        &(gbp*n*y(6)*y(7)**(n - 2)*(Kbp + y(3))*(n - 1))/(Kbp + y(3) + y(4)*y(7)**n)**2) - &
        &(2*w)/y(7)**3
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        do itr=1,mlyap
         j=itr*id
           do ix=1,id
            jsum=0.d0
             do iy=1,id
              jsum=jsum+jac(ix,iy)*Y(iy+j)
             enddo
            YDOT(ix+j)=jsum
           enddo
         enddo
         
      RETURN
      END SUBROUTINE DERIVS

         SUBROUTINE gs(y,bnorm)
         implicit none
         real(kind=8),intent(inout):: y(NEQ), bnorm(mlyap)
         real(kind=8):: a(mlyap)
         INTEGER:: itt,it,iii
         a=0.d0; bnorm=0.d0
         gs1:do itt=1,mlyap
         if(itt==1) then
         do it=1,id
             bnorm(itt)=bnorm(itt)+y(it+id)**2
         enddo
       !bnorm(itt)=sum(y(id+1:2*id))
         do it=1,id
              y(it+id)=y(it+id)/dsqrt(bnorm(itt))
         enddo
         !y(id+1:2*id)=y(id+1:2*id)/dsqrt(bnorm(itt))
         else
         gs2:do iii=1,(itt-1)
             a(iii)=0.d0
         do it=1,id
             a(iii)=a(iii)+y(it+iii*id)*y(it+itt*id)
         enddo
         do it=1,id
             y(it+itt*id)=y(it+itt*id)-a(iii)*y(it+iii*id)
         enddo
         enddo gs2
         do it=1,id
             bnorm(itt)=bnorm(itt)+y(it+itt*id)**2
         enddo
         do it=1,id
             y(it+itt*id)=y(it+itt*id)/dsqrt(bnorm(itt))
         enddo
         endif
         enddo gs1

         return
         end subroutine gs

        
!      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!     Subroutine to define the exact Jacobian for this problem
!        IMPLICIT NONE
!        INTEGER NEQ, ML, MU, NROWPD
!        DOUBLE PRECISION T, Y, PD
!        DIMENSION Y(NEQ), PD(NROWPD,NEQ)
!        PD=0.d0
        !N = Y(1)
        !PHI = Y(2)
        !PD(1,1) = -(ALPHA*PHI+BETA)
        !PD(2,1) = PHI*RHO + TAU
        !PD(1,2) = -N*ALPHA
        !PD(2,2) = RHO*N - SIGMA
!        RETURN

!      END SUBROUTINE JACD

    END MODULE CLIMATE



!******************************************************************

    PROGRAM DEMOCLIMATE
      USE CLIMATE
      USE DVODE_F90_M
!     Type declarations:

      IMPLICIT NONE
      real(kind=8):: t,TOUT,xe1, twopi, xlag1, xlag11, xlag2, xlag21, xlag3, xlag31, xlag4, xlag41, xperh, pstep
      integer(kind=8):: icount,j,looparm,i,inic,xint
      real(kind=8),allocatable:: Y(:),asum(:),bnorm(:),Ysave(:,:),Ysafe(:,:),xmax(:,:),xmin(:,:),xhistmx(:,:),xhistmn(:,:)
      real(kind=8),allocatable:: histmx_min(:),histmx_max(:),histmx_h(:),histmn_min(:),histmn_max(:),histmn_h(:)!,asum(:),bnorm(:)
      real(kind=8), allocatable :: simfx(:),xper(:),simfx1(:),simfx11(:), simfx2(:), simfx21(:)
      real(kind=8), allocatable :: simfx3(:),simfx31(:), simfx4(:), simfx41(:)
!!!!!!!!!!!!!! INTEGRATOR ENTRIES; USUALLY NOT TO BE TOUCHED !!!!!!     

      INTEGER ITASK, ISTATE, ISTATS, IOUT, SWITCH
      DOUBLE PRECISION RSTATS, RTOL, ATOL
      DIMENSION RSTATS(22), ISTATS(31)
      TYPE (VODE_OPTS) :: OPTIONS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      twopi=6.28318530718d0 !! value of 2pi

      allocate(Y(NEQ),xmax(id,3),xmin(id,3),xhistmx(id,nhist),xhistmn(id,nhist),histmx_min(id),histmx_max(id),histmx_h(id), &
      histmn_min(id),histmn_max(id),histmn_h(id),Ysave(id,niter),Ysafe(id,niter),xper(id),asum(mlyap),bnorm(mlyap))
      
      histmx_min=-0.2d0; histmx_max=5.d0
      histmn_min=-0.2d0; histmn_max=5.d0
      
      histmx_h=(histmx_max-histmx_min)/dfloat(nhist)
      histmn_h=(histmn_max-histmn_min)/dfloat(nhist)

!     Set  the problem parameters, you can change parameter values here:
       par(1) = 0.1d0 !! D  
        ! nutrient availability
        par(2) = 1.d0  !! S0   %0.08; %3
        ! Physiological prey parameters
        !gamma = par(4) = ;! %proportion of C:N->4.5 in marine bacteria
        par(3) = 0.8d0 !! n % n<1 increase grazing preasure on G. ecological cost of phenot. plat. (n=1 linear trade off, n=0.6 strong trade off, n=0.2 very strong trade off)e.g. increased overall grazing preasure n<1
        par(4) = 3.33d0 !! rpb
        par(5) = 3.33d0 !! rbb
        par(6) = 3.33d0 !! rmb
        par(7) = 0.01d0 !! Hpb
        par(8) = 0.01d0 !! Hbb
        par(9) = 0.01d0 !! Hmb
        par(10) =  0.3d0 !! ae assimilation efficieny
        ! Physiological predator parameters
        par(11) = 1.2d0  !! gpp %1.2
        par(12) = 1.2d0  !! gbp
        par(13) = 0.1d0; !! Kpp %0.1
        par(14) = 0.1d0; !! Kbp
        ! Generalist parameters
        !%c = 0.01; %1 in the paper Gaedke et al.% reduction of max growth compared with the specialist (direct costs)e.g. 0.2 means 20% less groth rate than specialist
        par(15) = 0.01d0  !! c %1 in the paper Gaedke et al.% reduction of max growth compared with the specialist (direct costs)e.g. 0.2 means 20% less groth rate than specialist
        par(16) = 0.0001d0 !! w  % Scaling of non adaptative trait adjustment
        par(17) =  0.7d0   !! V  (speed of adaptation (0.01 slow, 0.1 intemediate, 1 fast trait changes)

	   xint=200   !! TO CHANGE THE RESOLUTION OF THE CALCULATION
       pstep=1.d0/dfloat(xint)
    do looparm=1,xint+1; print*,looparm   !! parameter loop 
	xe1=pstep+(1.d0-pstep)/dfloat(xint)*dfloat(looparm-1); par(15)=xe1
!!!!!!!! Set the integration parameters:
      RTOL = 1.0D-6
      ATOL = 1.0D-9
      IOUT = 1
      ITASK = 1
      ISTATE = 1
     ! WRITE (6,*) RTOL, ATOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Set the initial conditions:
     outer: do inic=1,20
      call random_number(Y) !! random initial conditions in [0, 1]
      Y(1:id)=0.5d0*Y(1:id) !! random initial conditions in [0, 0.5]
      if (mlyap>0) call gs(Y,bnorm) 
      ISTATE = 1
!     Output points:
      T=0.d0
      TOUT=0.d0
!     Start with nonstiff method:
      OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL)
!     Perform the integration:
       asum=0.d0
      time: DO IOUT = 1, ntot
       TOUT = TOUT + h
       CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!!!!! ERROR ANALYSIS: DO NOT TOUCH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gather and write the integration statistics for this problem:
        CALL GET_STATS(RSTATS,ISTATS)
!        Stop the integration if an error occurred:
        IF (ISTATE<0) THEN
          WRITE (6,*) ISTATE
          print*, 'Convergence issue, trying with new initial conditions'
          exit time
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Print the solution and write the plot data file:
         if(IOUT>ntrans) then
               do i=1,id
                if (Y(i)<0.00000001d0) Y(i)=Y(i)+0.0000001d0
               enddo
           if (mlyap>0) then 
             call gs(Y,bnorm)
             asum=asum+0.5d0*dlog(bnorm)
           endif
           ISTATE=1
           Ysave(:,IOUT-ntrans)=Y(1:id) !! saving the signals in a big array for later manipulation
           !write(33,24)T,Y(1:id)
          else
          if (mlyap>0) call gs(Y,bnorm)
         endif  
      END DO time ! IOUT
      !print*,ISTATE
      IF (ISTATE>0) THEN
       !write(*,24) par(15),asum/(dfloat(niter)*h) !!! LEs
       if (mlyap>0) write(32,24) par(15),asum/(dfloat(niter)*h) !!! LEs
!!! MAX/MIN Estimation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     xhistmx=0.d0; xhistmn=0.d0
      xmax(:,1)=Ysave(1:id,1); xmax(:,2)=Ysave(1:id,2)
      xmin(:,1)=Ysave(1:id,1); xmin(:,2)=Ysave(1:id,2)
      do IOUT=3,niter
       xmax(:,3)=Ysave(1:id,IOUT)
       xmin(:,3)=Ysave(1:id,IOUT)
       do i=1,id
        if(xmax(i,2)>=xmax(i,1).and.xmax(i,2)>=xmax(i,3)) then !! maxima histogram entries
         do j=1,nhist
          if(xmax(i,2)>=histmx_min(i)+dfloat(j-1)*histmx_h(i).and.xmax(i,2)<=histmx_min(i)+dfloat(j)*histmx_h(i)) then
          xhistmx(i,j)=xhistmx(i,j)+1.d0
          endif 
         enddo 
        endif
        if(xmin(i,2)<=xmin(i,1).and.xmin(i,2)<=xmin(i,3)) then !! minima histogram entries
         do j=1,nhist
          if(xmin(i,2)>=histmn_min(i)+dfloat(j-1)*histmn_h(i).and.xmin(i,2)<=histmn_min(i)+dfloat(j)*histmn_h(i)) then
          xhistmn(i,j)=xhistmn(i,j)+1.d0
          endif 
         enddo
        endif
     enddo 
      xmax(:,1)=xmax(:,2); xmax(:,2)=xmax(:,3)
      xmin(:,1)=xmin(:,2); xmin(:,2)=xmin(:,3)
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!! detrending all the signals for lag calculation
    Ysafe=Ysave
    do i=1,id
     Ysave(i,:)=Ysave(i,:)-sum(Ysave(i,:))/dfloat(niter)
    enddo
!!!!! SIMILARITY FUNCTION CALCULATION for any 2 signals, lets say Y(idx1) and Y(idx2) here and sampling rate is 'h'
!! calculating signal periodicity
  allocate(simfx(imaxlag), simfx1(imaxlag), simfx11(imaxlag), simfx2(imaxlag), simfx21(imaxlag)) !! ALLOCATE DYNAMIC ARRAY IN THE CALCULATION
  allocate(simfx3(imaxlag), simfx31(imaxlag), simfx4(imaxlag), simfx41(imaxlag)) !!
  
   xper=0.d0
   do i=1,id
     simfx=0.d0 
     call similarity(Ysave(i,:),Ysave(i,:),simfx)
    do j=2,imaxlag
     if(simfx(j)<simfx(j-1).and.simfx(j)<simfx(j+1)) then
       xper(i)=dfloat(j)*h; exit
      endif
     enddo
   enddo

   do i=1,id
    xperh=sum(dabs(Ysave(i,1:int(xper(i)/h))))/dfloat(int(xper(i)/h))
    if(xperh<0.000001d0.or.sum(Ysafe(i,:))/dfloat(niter)<0.000001d0) xper(i)=0.d0
   enddo

   write(1,24) par(15), xper(1:id)  !! writing the period data for all the signals   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lag calculation  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lag2-5 calculation

 simfx1=0.d0; simfx11=0.d0; xlag1=0.d0; xlag11=0.d0
  call similarity(Ysave(idx1,:),Ysave(idx4,:),simfx1)
  
   if(xper(idx1)>h.and.xper(idx4)>h) then
    do j=2,imaxlag
     if(simfx1(j)<simfx1(j-1).and.simfx1(j)<simfx1(j+1)) then
       xlag1=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag1=dfloat(niter)*h
   endif
   
   call similarity(Ysave(idx4,:),Ysave(idx1,:),simfx11)
   
   if(xper(idx4)>h.and.xper(idx1)>h) then
    do j=2,imaxlag
     if(simfx11(j)<simfx11(j-1).and.simfx11(j)<simfx11(j+1)) then
       xlag11=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag11=dfloat(niter)*h; print*, "non-oscillatory signal encountered"
   endif
   if(xlag11<xlag1) xlag1=xlag11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lag3-6 calculation

 simfx2=0.d0; simfx21=0.d0; xlag2=0.d0; xlag21=0.d0
  call similarity(Ysave(idx2,:),Ysave(idx5,:),simfx2)
  
   if(xper(idx2)>h.and.xper(idx5)>h) then
    do j=2,imaxlag
     if(simfx2(j)<simfx2(j-1).and.simfx2(j)<simfx2(j+1)) then
       xlag2=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag2=dfloat(niter)*h
   endif
   
   call similarity(Ysave(idx5,:),Ysave(idx2,:),simfx21)
   
   if(xper(idx5)>h.and.xper(idx2)>h) then
    do j=2,imaxlag
     if(simfx21(j)<simfx21(j-1).and.simfx21(j)<simfx21(j+1)) then
       xlag21=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag21=dfloat(niter)*h; print*, "non-oscillatory signal encountered"
   endif
   if(xlag21<xlag2) xlag2=xlag21      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lag4-5 calculation

 simfx3=0.d0; simfx31=0.d0; xlag3=0.d0; xlag31=0.d0
  call similarity(Ysave(idx3,:),Ysave(idx4,:),simfx3)
  
   if(xper(idx3)>h.and.xper(idx4)>h) then
    do j=2,imaxlag
     if(simfx3(j)<simfx3(j-1).and.simfx3(j)<simfx3(j+1)) then
       xlag3=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag3=dfloat(niter)*h
   endif
   
   call similarity(Ysave(idx4,:),Ysave(idx3,:),simfx31)
   
   if(xper(idx4)>h.and.xper(idx3)>h) then
    do j=2,imaxlag
     if(simfx31(j)<simfx31(j-1).and.simfx31(j)<simfx31(j+1)) then
       xlag31=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag31=dfloat(niter)*h; print*, "non-oscillatory signal encountered"
   endif
   if(xlag31<xlag3) xlag3=xlag31
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lag4-6 calculation

 simfx4=0.d0; simfx41=0.d0; xlag4=0.d0; xlag41=0.d0
  call similarity(Ysave(idx3,:),Ysave(idx5,:),simfx4)
  
   if(xper(idx3)>h.and.xper(idx5)>h) then
    do j=2,imaxlag
     if(simfx4(j)<simfx4(j-1).and.simfx4(j)<simfx4(j+1)) then
       xlag4=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag4=dfloat(niter)*h
   endif
   
   call similarity(Ysave(idx5,:),Ysave(idx3,:),simfx41)
   
   if(xper(idx5)>h.and.xper(idx3)>h) then
    do j=2,imaxlag
     if(simfx41(j)<simfx41(j-1).and.simfx41(j)<simfx41(j+1)) then
       xlag41=dfloat(j)*h; exit
      endif
     enddo
   else
    xlag41=dfloat(niter)*h; print*, "non-oscillatory signal encountered"
   endif
   if(xlag41<xlag4) xlag4=xlag41      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(2,24) par(15), xlag1/xper(idx1), xlag2/xper(idx2), xlag3/xper(idx3), xlag4/xper(idx3) !! writing the lag data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!! writing the bifurcation data from histograms
    do i=1,id
      do j=1,nhist
       if(xhistmx(i,j)>0.0000001d0) then
        write(10+i,24)par(15), histmx_min(i)+dfloat(j)*histmx_h(i)
       endif
       if(xhistmn(i,j)>0.0000001d0) then
        write(20+i,24)par(15), histmn_min(i)+dfloat(j)*histmn_h(i)
       endif
      enddo
     enddo 
       deallocate(simfx,simfx1,simfx11,simfx2,simfx21,simfx3,simfx31,simfx4,simfx41)
      END IF
     enddo outer
    enddo   !!! looparm   

      24  Format(10F20.10)
!     Write the integration final integration statistics:
      WRITE (6,*) ISTATS(11), ISTATS(12), ISTATS(13)
!     Format statements for this problem:

END PROGRAM DEMOCLIMATE

