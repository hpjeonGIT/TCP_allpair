!
! ####### Pure Verlet time integration routine with No Thermostat #########
! 
! Purely integration routine with no thermostat. Basically, this routine is
! applied when new potential is implemented - Any potential should not increase
! the given temperature or increase temperature as iteration goes. 
! Also for NVE (microcanonical ensemble) simulation.
!
SUBROUTINE VVerletStoch1(NS, qel, qion, param, dt, gsend)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     USE DATASTR
     IMPLICIT NONE
     REAL(KIND=DP):: fluct, x, r, v1, v2
     REAL(KIND=DP):: rand1, rand2, ran2
     REAL:: ranf
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(GH)::gsend(NS%Ng)
TYPE(PM)::param
REAL(KIND=DP)  ::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL(KIND=DP) :: x, xm
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k)*(1.D0-0.5D0*param%alpha*dt/xm) + &
           & 0.5D0*dt*qel(i)%ff(k)/xm
      x = qel(i)%xx(k) + dt*qel(i)%xv(k)
      qel(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      gsend(i)%xx(k) = qel(i)%xx(k)
   END DO
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k)*(1.D0-0.5D0*param%alpha*dt/xm) + &
           & 0.5D0*dt*qion(i)%ff(k)/xm
      x = qion(i)%xx(k) + dt*qion(i)%xv(k)
      qion(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      gsend(i+NS%Nel)%xx(k) = qion(i)%xx(k)
   END DO
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch1
!
! 2nd Verlet routine for stochastic thermostat ###############################
SUBROUTINE VVerletStoch2(NS, qel, qion, param, sys, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     USE DATASTR
     IMPLICIT NONE
     REAL(KIND=DP):: fluct, x, r, v1, v2
     REAL(KIND=DP):: rand1, rand2, ran2
     REAL:: ranf
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
REAL(KIND=DP)  ::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, id
REAL(KIND=DP) :: xm, sigma, beta, eta, mv2
sigma = 1.D0       ! Deviation for Gaussian distribution
sys%mv2el = 0.0D0
sys%mv2ion = 0.0D0
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   beta = DSQRT(2.D0*param%alpha*param%Te/dt) ! Temperature unit is eV
   DO k = 1, 3
      eta = fluct(sigma)
      qel(i)%ff(k) = qel(i)%ff(k) + eta*beta
      qel(i)%xv(k) = (qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm) / &
           & (1.D0 + 0.5D0*param%alpha*dt/xm)
   END DO
   mv2 = (qel(i)%xv(1)**2 + qel(i)%xv(2)**2 + qel(i)%xv(3)**2)*xm
   sys%mv2el = sys%mv2el + mv2
   id = INT(mv2/3./NS%dv_e) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_e(id) = NS%VDF_e(id) + 1.D0
   END IF
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   beta = DSQRT(2.*param%alpha*param%Ti/dt) ! Temperature unit is eV
   DO k = 1, 3
      eta = fluct(sigma)
      qion(i)%ff(k) = qion(i)%ff(k) + eta*beta
      qion(i)%xv(k) = (qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm) / &
           & (1.D0 + 0.5D0*param%alpha*dt/xm)
   END DO
   mv2 = (qion(i)%xv(1)**2 + qion(i)%xv(2)**2 + qion(i)%xv(3)**2)*xm
   sys%mv2ion = sys%mv2ion + mv2
   id = INT(mv2/3./NS%dv_i) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_i(id) = NS%VDF_i(id) + 1.D0
   END IF
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch2
!
! Verlet routine without thermostat - for microcanonical ensemble (NVE) ######
! ############################################################################
! Also can be used for isokinetic/berendsen thermostat as 1st verlet routine
SUBROUTINE VVerletNotemp1(NS, qel, qion, param, dt, gsend)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(GH)::gsend(NS%Ng)
TYPE(PM)::param
REAL(KIND=DP)  ::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL(KIND=DP) :: x, xm
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
      x = qel(i)%xx(k) + dt*qel(i)%xv(k)
      qel(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      gsend(i)%xx(k) = qel(i)%xx(k)
   END DO
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
      x = qion(i)%xx(k) + dt*qion(i)%xv(k)
      qion(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))     
      gsend(i+NS%Nel)%xx(k) = qion(i)%xx(k)
   END DO
END DO
!
RETURN
!
END SUBROUTINE VVerletNotemp1
!
! 2nd routine for Verlet routine - without thermostat #########################
SUBROUTINE VVerletNotemp2(NS, qel, qion, param, sys, dt)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
REAL(KIND=DP)  ::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, id
REAL(KIND=DP) :: xm, mv2
!
!
sys%mv2el  = 0.0D0
sys%mv2ion = 0.0D0
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
   END DO
   mv2 = (qel(i)%xv(1)**2 + qel(i)%xv(2)**2 + qel(i)%xv(3)**2)*xm
   sys%mv2el = sys%mv2el + mv2
   id = INT(mv2/3.D0/NS%dv_e) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_e(id) = NS%VDF_e(id) + 1.D0
   END IF
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
   END DO
   mv2 = (qion(i)%xv(1)**2 + qion(i)%xv(2)**2 + qion(i)%xv(3)**2)*xm
   sys%mv2ion = sys%mv2ion + mv2
   id = INT(mv2/3.D0/NS%dv_i) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_i(id) = NS%VDF_i(id) + 1.D0
   END IF
END DO
!
RETURN
!
END SUBROUTINE VVerletNotemp2
!
! Berendsen thermostat ###################################################
! 1st routine is shared with no thermostat verlet integration routine 
SUBROUTINE VVerletBerend2(NS, qel, qion, param, sys, dt, COMM)
USE DataStr
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
REAL(KIND=DP)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, IERR, id
REAL(KIND=DP) :: xm, lambda, T0, mv2
!
!
sys%mv2el  = 0.0D0
sys%mv2ion = 0.0D0
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
      sys%mv2el = sys%mv2el + xm*qel(i)%xv(k)**2
   END DO
END DO
CALL MPI_ALLREDUCE(sys%mv2el, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(3*NS%NTel)
lambda = DSQRT(1.D0+ dt*(param%Te/T0 - 1.D0)/param%tau)
sys%mv2el  = 0.0D0
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k)*lambda
   END DO
   mv2 = (qel(i)%xv(1)**2 + qel(i)%xv(2)**2 + qel(i)%xv(3)**2)*xm
   sys%mv2el = sys%mv2el + mv2
   id = INT(mv2/3.D0/NS%dv_e) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_e(id) = NS%VDF_e(id) + 1.D0
   END IF
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
      sys%mv2ion = sys%mv2ion + xm*qion(i)%xv(k)**2
   END DO
END DO
CALL MPI_ALLREDUCE(sys%mv2ion, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(3*NS%NTion)
lambda = DSQRT(1.D0+ dt*(param%Ti/T0 - 1.D0)/param%tau)
sys%mv2ion = 0.0D0
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k)*lambda
   END DO
   mv2 = (qion(i)%xv(1)**2 + qion(i)%xv(2)**2 + qion(i)%xv(3)**2)*xm
   sys%mv2ion = sys%mv2ion + mv2
   id = INT(mv2/3.D0/NS%dv_i) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_i(id) = NS%VDF_i(id) + 1.
   END IF
END DO
!
RETURN
!
END SUBROUTINE VVerletBerend2
!
! Isokinetic thermostat - linear velocity scaling thermostat ##############
SUBROUTINE VVerletIsokin2(NS, qel, qion, param, sys, dt, COMM)
USE DataStr
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
REAL(KIND=DP)  ::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k, IERR, id
REAL(KIND=DP) :: xm, lambda, T0, mv2
!
!
sys%mv2el  = 0.0D0
sys%mv2ion = 0.0D0
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
      sys%mv2el = sys%mv2el + xm*qel(i)%xv(k)**2
   END DO
END DO
CALL MPI_ALLREDUCE(sys%mv2el, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(3*NS%NTel)
lambda = DSQRT(param%Te/T0)
sys%mv2el  = 0.0D0
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k)*lambda
   END DO
   mv2 = (qel(i)%xv(1)**2 + qel(i)%xv(2)**2 + qel(i)%xv(3)**2)*xm
   sys%mv2el = sys%mv2el + mv2
   id = INT(mv2/3.D0/NS%dv_e) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_e(id) = NS%VDF_e(id) + 1.
   END IF
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
      sys%mv2ion = sys%mv2ion + xm*qion(i)%xv(k)**2
   END DO
END DO
CALL MPI_ALLREDUCE(sys%mv2ion, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(3*NS%NTion)
lambda = DSQRT(param%Ti/T0)
sys%mv2ion = 0.0D0
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k)*lambda
   END DO
   mv2 = (qion(i)%xv(1)**2 + qion(i)%xv(2)**2 + qion(i)%xv(3)**2)*xm
   sys%mv2ion = sys%mv2ion + mv2
   id = INT(mv2/3./NS%dv_i) + 1
   IF ( id <= Ndist + 1 ) THEN
      NS%VDF_i(id) = NS%VDF_i(id) + 1.D0
   END IF
END DO
!
RETURN
!
END SUBROUTINE VVerletIsokin2
!
! ############## Random Gaussian(Normal) Distribution Function ################
! 
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the 
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x 
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
! 
FUNCTION fluct(x)
  USE DATASTR
  IMPLICIT NONE
  REAL(KIND=DP):: fluct, x, r, v1, v2
  REAL(KIND=DP):: rand1, rand2
  !
  ! Initialization
  r=1.
  DO WHILE (r.ge.1.)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.D0*rand1 - 1.D0
     v2 = 2.D0*rand2 - 1.D0
     r = v1*v1+v2*v2
  END DO
  fluct = v1*DSQRT(-2.*DLOG(r)/r)*x
  RETURN
END FUNCTION fluct
!
! Velocity autocorrelation function and diffusion coefficient estimation #####
! ############################################################################
SUBROUTINE DIFFUSION(NS, qel, qion)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL(KIND=DP) :: mean_D
!
! Electron VAF
mean_D = 0.0D0
DO i=1, NS%Nel
   DO k=1, 3      
      mean_D = mean_D + qel(i)%xx(k)**2
   END DO
END DO
!sys%Del = mean_D
!
! Ion VAF
mean_D = 0.0D0
DO i=1, NS%Nion
   DO k=1, 3
      mean_D = mean_D + qion(i)%xx(k)**2
   END DO
END DO
!sys%Dion = mean_D
!
RETURN
END SUBROUTINE DIFFUSION
!
! Velocity autocorrelation function and diffusion coefficient estimation #####
! ############################################################################
SUBROUTINE VAF(NS, qel, qion, sys)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(ST)::sys
!
! INTERNAL VARIABLES
INTEGER:: i
REAL(KIND=DP) :: ve(3), vi(3)
!
! Electron VAF
sys%cel = 0.0D0
ve = 0.0D0
vi = 0.0D0
DO i=1, NS%Nel
   sys%cel = sys%cel + qel(i)%xv(1)*qel(i)%vi(1) + &
        & qel(i)%xv(2)*qel(i)%vi(2) + qel(i)%xv(3)*qel(i)%vi(3)
   ve(:) = ve(:) + qel(i)%xv(:)
END DO
!
! Ion VAF
sys%cion = 0.0D0
DO i=1, NS%Nion
   sys%cion = sys%cion + qion(i)%xv(1)*qion(i)%vi(1) + &
        & qion(i)%xv(2)*qion(i)%vi(2) + qion(i)%xv(3)*qion(i)%vi(3)
   vi(:) = vi(:) + qion(i)%xv(:)
END DO
!
! Mutual VAF
sys%Jei(:) = (ve(:) - vi(:))*0.5D0
!sys%Jie(:) = (vi(:) - ve(:))*0.5
!
! For diffusion coefficient, they will be implemented later - 122006
RETURN
END SUBROUTINE VAF
