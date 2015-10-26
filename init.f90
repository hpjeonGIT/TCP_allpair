!
! Particle trace initialization ##############################################
SUBROUTINE Rinit(NS, qel, qion, param)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
!
INTEGER:: i
!
DO i = 1, NS%Nel
   qel(i)%xx(:)  = qel(i)%xx(:)-param%box(:)*DNINT(qel(i)%xx(:)/param%box(:))
END DO
DO i = 1, NS%Nion
   qion(i)%xx(:) = qion(i)%xx(:)-param%box(:)*DNINT(qion(i)%xx(:)/param%box(:))
END DO
!
RETURN
END SUBROUTINE Rinit
!
! ################ From Frenkel and Smit ####################
! 
SUBROUTINE Vinit(NS, qel, qion, param, sys, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
!
INTEGER:: i, j, IERR
REAL(KIND=DP) :: xm, lambda, T0, xv(NS%Npt,3), rand, sumv(3), mv2, ve(3), vi(3), J0(3)
!
! Electron velocity initialization
sys%mv2el = 0.0D0
xm = param%xm(1)
sumv = 0.0D0
DO i = 1, NS%Nel
   DO j=1, 3
      CALL RANDOM_NUMBER(rand)
      xv(i,j) = rand - 0.5D0
      sumv(j) = sumv(j) + xv(i,j)
   END DO
END DO
!
sumv = sumv/DBLE(NS%Nel)
DO i=1, NS%Nel
   xv(i,:)= xv(i,:) - sumv(:)
   sys%mv2el = sys%mv2el + xm*(xv(i,1)**2 + xv(i,2)**2 + xv(i,3)**2)
END DO
CALL MPI_ALLREDUCE(sys%mv2el, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(NS%NTel*3)
lambda = DSQRT(param%Te/T0)
sys%mv2el = 0.0D0
DO i = 1, NS%Nel
   DO j = 1,3
      qel(i)%xv(j) = xv(i,j)*lambda
      qel(i)%ff(j) = 0.0D0
      sys%mv2el = sys%mv2el + xm*qel(i)%xv(j)**2
   END DO
END DO
!
! Ion velocity initialization
sys%mv2ion = 0.0D0
xm = param%xm(2)
sumv = 0.0D0
DO i = 1, NS%Nion
   DO j=1, 3
      CALL RANDOM_NUMBER(rand)
      xv(i,j) = rand - 0.5D0
      sumv(j) = sumv(j) + xv(i,j)
   END DO
END DO
!
sumv = sumv/DBLE(NS%Nion)
DO i=1, NS%Nion
   xv(i,:)= xv(i,:) - sumv(:)
   sys%mv2ion = sys%mv2ion + xm*(xv(i,1)**2 + xv(i,2)**2 + xv(i,3)**2)
END DO
CALL MPI_ALLREDUCE(sys%mv2ion, mv2, 1, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
T0 = mv2/DBLE(NS%NTion*3)
lambda = SQRT(param%Ti/T0)
sys%mv2ion = 0.0D0
DO i = 1, NS%Nion
   DO j = 1,3
      qion(i)%xv(j) = xv(i,j)*lambda
      qion(i)%ff(j) = 0.0D0
      sys%mv2ion = sys%mv2ion + xm*qion(i)%xv(j)**2
   END DO
END DO
!
! Initial velocity flux
ve = 0.0D0
vi = 0.0D0
DO i=1, NS%Nel
   ve(:) = ve(:) + qel(i)%xv(:)
END DO
DO i=1, NS%Nion
   vi(:) = vi(:) + qion(i)%xv(:)
END DO
J0(:) = (ve(:) - vi(:))*0.5D0
CALL MPI_ALLREDUCE(J0, sys%J0, 3, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
!
RETURN
END SUBROUTINE Vinit
!
! Particle initial velocity initialization ###################################
SUBROUTINE Vinit_Again(NS, qel, qion, sys, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(ST)::sys
TYPE(MP)::COMM
!
INTEGER:: i, IERR
REAL(KIND=DP) :: ve(3), vi(3), J(3)
!
DO i = 1, NS%Nel
   qel(i)%vi(:)  = qel(i)%xv(:)
END DO
DO i = 1, NS%Nion
   qion(i)%vi(:) = qion(i)%xv(:)
END DO
!
! Initial velocity flux
ve = 0.0
vi = 0.0
DO i=1, NS%Nel
   ve(:) = ve(:) + qel(i)%xv(:)
END DO
DO i=1, NS%Nion
   vi(:) = vi(:) + qion(i)%xv(:)
END DO
J(:) = (ve(:) - vi(:))*0.5D0
CALL MPI_ALLREDUCE(J, sys%J0, 3, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
!
RETURN
END SUBROUTINE Vinit_Again
!
! Particle initial force initialization ###################################
SUBROUTINE Finit(NS, qel, qion, param, sys, COMM, gsend, gemit, grec, pme, &
     & PSI, dt, PRE_data)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(ST)::sys
TYPE(MP)::COMM
TYPE(GH)::gsend(NS%Ng), gemit(NS%Ng), grec(NS%Ng), pme(NS%NLpt)
REAL(KIND=DP):: PSI(param%K(1), param%K(2), param%K(3))
REAL(KIND=DP):: dt, PRE_data(Ngrid,4)
!
INTEGER:: i, k, Ne, Ni, NUMITER, IERR, ISTATUS(MPI_STATUS_SIZE)
REAL(KIND=DP) :: x, Ur
!
DO i = 1, NS%Nel
   DO k = 1, 3
      x = qel(i)%xx(k) - dt*qel(i)%xv(k)
      qel(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      gsend(i)%xx(k) = qel(i)%xx(k)
   END DO
END DO
!
DO i = 1, NS%Nion
   DO k = 1, 3
      x = qion(i)%xx(k) - dt*qion(i)%xv(k)
      qion(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))     
      gsend(i+NS%Nel)%xx(k) = qion(i)%xx(k)
   END DO
END DO

IF (COMM%ITAG) THEN
   NUMITER = INT(COMM%NUMPROC/2) - 1
ELSE 
   NUMITER = INT(COMM%NUMPROC/2)
END IF

CALL Force_SPME4(NS, param, sys, COMM, PSI, gsend, pme)
IF (COMM%TAG) THEN
   CALL Force_direct_self(NS, qel, qion, param, sys, PRE_data)
   DO i=1, NUMITER
      Ne  = NS%NOPT(COMM%SRC(i)+1,1)
      Ni = NS%NOPT(COMM%SRC(i)+1,2)
      CALL MPI_SENDRECV(gsend, NS%Npt, COMM%PTCL, COMM%DEST(i), 1, &
           & grec, Ne+Ni, COMM%PTCL, COMM%SRC(i), 1, &
           & COMM%COMM1D, ISTATUS, IERR)
      CALL Force_direct_ext(NS, qel, qion, param, grec, gemit, &
           & Ur, Ne, Ni, PRE_data)
      CALL MPI_SENDRECV(gemit, Ne+Ni, COMM%PTCL, COMM%SRC(i), 0, &
           & grec, NS%Npt, COMM%PTCL, COMM%DEST(i), 0, &
           & COMM%COMM1D, ISTATUS, IERR)
      CALL Force_direct_sum(NS, qel, qion, grec)
   END DO
   IF (COMM%ITAG) THEN
      Ne  = NS%NOPT(COMM%SRC(i)+1,1)
      Ni = NS%NOPT(COMM%SRC(i)+1,2)
      CALL MPI_SENDRECV(gsend, NS%Npt, COMM%PTCL, COMM%DEST(i), 1, &
           & grec, Ne+Ni, COMM%PTCL, COMM%SRC(i), 1, &
           & COMM%COMM1D, ISTATUS, IERR)
      CALL Force_direct_ext(NS, qel, qion, param, grec, gemit, &
           & Ur, Ne, Ni, PRE_data)
   END IF
END IF
   CALL SPME_SUM(NS, qel, qion, sys, COMM, pme, grec)
   IF (COMM%TAG) CALL REMOVE_RIGIDMTN(NS, qel, qion, COMM)


DO i = 1, NS%Nel
   DO k = 1, 3
      x = qel(i)%xx(k) + dt*qel(i)%xv(k)
      qel(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))
      gsend(i)%xx(k) = qel(i)%xx(k)
   END DO
END DO
!
DO i = 1, NS%Nion
   DO k = 1, 3
      x = qion(i)%xx(k) + dt*qion(i)%xv(k)
      qion(i)%xx(k) = x - param%box(k)*DNINT(x/param%box(k))     
      gsend(i+NS%Nel)%xx(k) = qion(i)%xx(k)
   END DO
END DO

END SUBROUTINE Finit
