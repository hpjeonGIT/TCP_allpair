!
! ################## Self-interaction term of EWALD sum #######################
! 
SUBROUTINE EWALD_SELF(NS, param, sys, PSI, COMM)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
REAL(KIND=DP)  :: PSI(param%K(1), param%K(2), param%K(3))
!
! Internal variables
INTEGER:: i, j, k, ix, iy, iz, m
REAL(KIND=DP) :: m2, rc_cut
!
!
rc_cut = MIN(MIN(param%box(1),param%box(2)),param%box(3))/2.0D0
IF (param%rc > rc_cut) THEN
   CALL MPI_Finalize(COMM%rc)
   STOP " *** error at cut-off radius *** "
END IF
param%a = 3.1D0/param%rc
param%rc2 = param%rc**2
!
!
IF (.NOT.COMM%TAG) THEN
   DO i=1, param%K(1)
      ix = MOD(i-1+param%K(1)/2, param%K(1)) - param%K(1)/2
      DO j=1, param%K(2)
         iy = MOD(j-1+param%K(2)/2, param%K(2)) - param%K(2)/2
         DO k=1, param%K(3)
            iz = MOD(k-1+param%K(3)/2, param%K(3)) - param%K(3)/2
            m = ix**2 + iy**2 + iz**2
            IF( m /= 0) THEN
               m2 = DBLE(ix)**2/param%box(1)**2 + DBLE(iy)**2/param%box(2)**2+&
                    & DBLE(iz)**2/param%box(3)**2
               PSI(i,j,k) = EXP(-m2*pi**2/param%a**2)/m2
            END IF
         END DO
      END DO
   END DO
ELSE
   sys%Uo = - DBLE(NS%Nel)*param%q(1)**2  - DBLE(NS%Nion)*param%q(2)**2
   sys%Uo = sys%Uo * param%a/sqrtpi
END IF
!
RETURN
END SUBROUTINE EWALD_SELF
!
! ################# Kelbg potential for coulomb interaction ###################
SUBROUTINE Force_direct_self(NS, qel, qion, param, sys)
USE DATASTR
IMPLICIT NONE
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg 
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(PM):: param
TYPE(ST):: sys
!
!
INTEGER:: i, j, k, idx
REAL(KIND=DP) :: r2, xr(3), df, r, q1, q2, x, y
!
!
DO i=1, NS%Nel
   qel(i)%ff(:) = 0.0D0
END DO
DO i=1, NS%Nion
   qion(i)%ff(:) = 0.0D0
END DO
sys%Ur = 0.0D0
!
!
! Direct - real space sum for EWALD energy and force 
!
! Direct sum for electron-electron
q1 = param%q(1)
q2 = param%q(1)
DO i=1, NS%Nel-1
   DO j=i+1, NS%Nel
      r2 = 0.0D0
      DO k=1,3
         x = qel(i)%xx(k) - qel(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      IF (r2 < param%rc2) THEN
         ! RDF
         !
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ee)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ee(idx) = NS%RDF_ee(idx) + 2.D0
         END IF
         x = 1.D0 - DERF(param%a*r)
         sys%Ur = sys%Ur + q1*q2*x/r
         df = eps*q1*q2*(2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + x)/r2/r
         DO k=1,3
            qel(i)%ff(k) = qel(i)%ff(k) + df*xr(k)
            qel(j)%ff(k) = qel(j)%ff(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
! 
! Direct sum for electron-ion
q1 = param%q(1)
q2 = param%q(2)
DO i=1, NS%Nel
   DO j=1, NS%Nion
      r2 = 0.0D0
      DO k=1,3
         x = qel(i)%xx(k) - qion(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      !
      IF (r2 < param%rc2) THEN
         ! RDF
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ei)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ei(idx) = NS%RDF_ei(idx) + 2.D0
         END IF
         x = 1.D0 - DERF(param%a*r)
         y = DEXP(-r/param%rs)
         sys%Ur = sys%Ur + q1*q2*(x-y)/r
         df = eps*q1*q2*(( 2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + &
              & x - y)/r2/r - y/r2/param%rs)

         DO k=1,3
            qel(i)%ff(k)  = qel(i)%ff(k)  + df*xr(k)
            qion(j)%ff(k) = qion(j)%ff(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
! Direct sum for ion-ion
q1 = param%q(2)
q2 = param%q(2)
DO i=1, NS%Nion-1
   DO j=i+1, NS%Nion
      r2 = 0.0D0
      DO k=1,3
         x = qion(i)%xx(k) - qion(j)%xx(k)
         xr(k) = x - param%box(k)*ANINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      !
      IF (r2 < param%rc2) THEN
         ! RDF
         r = SQRT(r2)
         idx = INT(r/NS%dx_ii)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ii(idx) = NS%RDF_ii(idx) + 2.D0
         END IF
         x = 1.D0 - DERF(param%a*r)
         sys%Ur = sys%Ur + q1*q2*x/r 
         df = eps*q1*q2*(2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + x)/r2/r

         DO k=1,3
            qion(i)%ff(k) = qion(i)%ff(k) + df*xr(k)
            qion(j)%ff(k) = qion(j)%ff(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
RETURN
END SUBROUTINE Force_direct_self
!
! ############## Routine for direct sum with external particles ###############
SUBROUTINE Force_direct_ext(NS, qel, qion, param, grec, gemit, &
     & Ur, Ne, Ni)
USE DATASTR
IMPLICIT NONE
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg 
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(GH):: grec(NS%Ng), gemit(NS%Ng)
TYPE(PM):: param
REAL(KIND=DP):: Ur
INTEGER :: Ne, Ni
!
!
INTEGER:: i, j, k, idx
REAL(KIND=DP) :: r2, r, xr(3), df, q1, q2, x, y
Ur = 0.0D0

gemit(:)%xx(1) = 0.0D0
gemit(:)%xx(2) = 0.0D0
gemit(:)%xx(3) = 0.0D0
!
!
! Direct sum for electron-electron_g
q1 = param%q(1)
q2 = param%q(1)
DO i=1, NS%Nel
   DO j=1, Ne
      r2 = 0.0D0
      DO k=1,3
         x = qel(i)%xx(k) - grec(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      IF (r2 < param%rc2) THEN
         !
         ! RDF
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ee)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ee(idx) = NS%RDF_ee(idx) + 2.0D0
         END IF
         x  = 1.0D0 - DERF(param%a*r)
         Ur = Ur + q1*q2*x/r
         df = eps*q1*q2*(2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + x)/r2/r
         DO k=1,3
            qel(i)%ff(k)   = qel(i)%ff(k)   + df*xr(k)
            gemit(j)%xx(k) = gemit(j)%xx(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
! 
! Direct sum for electron-ion_g
q1 = param%q(1)
q2 = param%q(2)
DO i=1, NS%Nel
   DO j=Ne+1, Ne+Ni
      r2 = 0.0D0
      DO k=1,3
         x = qel(i)%xx(k) - grec(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      IF (r2 < param%rc2) THEN
         !
         ! RDF
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ei)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ei(idx) = NS%RDF_ei(idx) + 2.D0
         END IF
         x  = 1.D0 - DERF(param%a*r)
         y  = DEXP(-r/param%rs)
         Ur = Ur + q1*q2*(x - y)/r
         df = eps*q1*q2*(( 2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + &
              & x - y)/r2/r - y/r2/param%rs)

         DO k=1,3
            qel(i)%ff(k)   = qel(i)%ff(k)  + df*xr(k)
            gemit(j)%xx(k) = gemit(j)%xx(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
! 
!Direct sum for ion-electron_g
q1 = param%q(2)
q2 = param%q(1)
DO i=1, NS%Nion
   DO j=1, Ne
      r2 = 0.0D0
      DO k=1,3
         x = qion(i)%xx(k) - grec(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      IF (r2 < param%rc2) THEN
         !
         ! RDF
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ei)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ei(idx) = NS%RDF_ei(idx) + 2.D0
         END IF
         x  = 1.D0 - DERF(param%a*r)
         y  = DEXP(-r/param%rs)
         Ur = Ur + q1*q2*(x - y)/r
         df = eps*q1*q2*(( 2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + &
              & x - y)/r2/r - y/r2/param%rs)

         DO k=1,3
            qion(i)%ff(k)  = qion(i)%ff(k)  + df*xr(k)
            gemit(j)%xx(k) = gemit(j)%xx(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
!
!Direct sum for ion-ion_g
q1 = param%q(2)
q2 = param%q(2)
DO i=1, NS%Nion
   DO j=Ne+1, Ne+Ni
      r2 = 0.0D0
      DO k=1,3
         x = qion(i)%xx(k) - grec(j)%xx(k)
         xr(k) = x - param%box(k)*DNINT(x/param%box(k))
         r2 = r2 + xr(k)**2
      END DO
      !
      IF (r2 < param%rc2) THEN
         !
         ! RDF
         r = DSQRT(r2)
         idx = INT(r/NS%dx_ii)+1
         IF ( idx <= Ndist+1) THEN
            NS%RDF_ii(idx) = NS%RDF_ii(idx) + 2.D0
         END IF
         x  = 1.D0 - DERF(param%a*r)
         Ur = Ur + q1*q2*x/r 
         df = eps*q1*q2*(2.D0*DEXP(-param%a**2*r2)*param%a*r/sqrtpi + x)/r2/r

         DO k=1,3
            qion(i)%ff(k)  = qion(i)%ff(k)  + df*xr(k)
            gemit(j)%xx(k) = gemit(j)%xx(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
RETURN
END SUBROUTINE Force_direct_ext
!
! ################ Routine for external direct force sum #####################
SUBROUTINE Force_direct_sum(NS, qel, qion, grec)
USE DATASTR
IMPLICIT NONE
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(GH):: grec(NS%Ng)
!
! INTERNAL VARIABLES
INTEGER:: i, k
!
! Electron force sum
DO i=1, NS%Nel
   DO k=1, 3
      qel(i)%ff(k) = qel(i)%ff(k) + grec(i)%xx(k)
   END DO
END DO
!
! Ion force sum
DO i=1, NS%Nion
   DO k=1, 3
      qion(i)%ff(k) = qion(i)%ff(k) + grec(i+NS%Nel)%xx(k)
   END DO
END DO
!
RETURN
END SUBROUTINE Force_direct_sum
!
! ################ Routine for SPME ###########################################
! ################ Using FFTE library #########################################
!
SUBROUTINE Force_SPME(NS, param, sys, COMM, PSI, gsend, pme)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg 
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(GH):: gsend(NS%Ng), pme(NS%NLpt)
REAL(KIND=DP)  :: PSI(param%K(1), param%K(2), param%K(3))
!
!
INTEGER :: i, j, k, l, m, n, x, ix, iy, iz, IERR
REAL(KIND=DP)  :: ux, uy, uz, V, Epme, AX(4), BX(4), CX(4), q1, z(3), &
     AXX(4), BXX(4), CXX(4)
REAL(KIND=DP),  ALLOCATABLE:: dQdx(:,:,:,:,:), Qtmp(:,:,:), Qsum(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE::Q(:,:,:), Qinv(:,:,:), Qfin(:,:,:), Qalg(:,:,:)
!
CALL MPI_GATHERV(gsend, NS%CNT(COMM%LID+1), COMM%PTCL, pme, NS%CNT, &
     & NS%DPL, COMM%PTCL, 0, COMM%LOCAL, IERR)
!
IF (.NOT. COMM%TAG) THEN
   !
   ! Initialization and build Q matrix
   ALLOCATE(dQdx(4, 4, 4, NS%NTpt,3), &
        & Qtmp(param%K(1), param%K(2), param%K(3)), &
        & Qsum(param%K(1), param%K(2), param%K(3)), &
        & Q(param%K(1), param%K(2), param%K(3)), &
        & Qinv(param%K(1), param%K(2), param%K(3)), &
        & Qfin(param%K(1), param%K(2), param%K(3)), &
        & Qalg(param%K(1), param%K(2), param%K(3)))
   V = param%box(1)*param%box(2)*param%box(3)
   Q = 0.0D0
   Qtmp = 0.0D0
   m = 0
   q1 = param%q(1)
   DO l=2, COMM%Nlsize
      n = COMM%iy*(COMM%Ncolmn-1) + l - 1
      DO x=1, NS%NOPT(n,1)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         z(1) = DINT(ux) - ux + 1.D0
         z(2) = DINT(uy) - uy + 1.D0
         z(3) = DINT(uz) - uz + 1.D0
         AX(1) =      z(1)**3                                /6.D0
         AXX(1)=              -  3.*z(1)**2                  /6.D0
         z(1) = z(1) + 1.D0
         AX(2) = (-3.*z(1)**3 + 12.*z(1)**2 - 12.*z(1) +  4.)/6.D0
         AXX(2)= (               9.*z(1)**2 - 24.*z(1) + 12.)/6.D0
         z(1) = z(1) + 1.D0
         AX(3) = ( 3.*z(1)**3 - 24.*z(1)**2 + 60.*z(1) - 44.)/6.D0
         AXX(3)= (            -  9.*z(1)**2 + 48.*z(1) - 60.)/6.D0
         z(1) = z(1) + 1.D0
         AX(4) = (   -z(1)**3 + 12.*z(1)**2 - 48.*z(1) + 64.)/6.D0
         AXX(4)= (               3.*z(1)**2 - 24.*z(1) + 48.)/6.D0
         BX(1) =      z(2)**3                                /6.D0
         BXX(1)=              -  3.*z(2)**2                  /6.D0
         z(2) = z(2) + 1.D0
         BX(2) = (-3.*z(2)**3 + 12.*z(2)**2 - 12.*z(2) +  4.)/6.D0
         BXX(2)= (               9.*z(2)**2 - 24.*z(2) + 12.)/6.D0
         z(2) = z(2) + 1.D0
         BX(3) = ( 3.*z(2)**3 - 24.*z(2)**2 + 60.*z(2) - 44.)/6.D0
         BXX(3)= (            -  9.*z(2)**2 + 48.*z(2) - 60.)/6.D0
         z(2) = z(2) + 1.D0
         BX(4) = (   -z(2)**3 + 12.*z(2)**2 - 48.*z(2) + 64.)/6.D0
         BXX(4)= (               3.*z(2)**2 - 24.*z(2) + 48.)/6.D0
         CX(1) =      z(3)**3                                /6.D0
         CXX(1)=              -  3.*z(3)**2                  /6.D0
         z(3) = z(3) + 1.D0
         CX(2) = (-3.*z(3)**3 + 12.*z(3)**2 - 12.*z(3) +  4.)/6.D0
         CXX(2)= (               9.*z(3)**2 - 24.*z(3) + 12.)/6.D0
         z(3) = z(3) + 1.D0
         CX(3) = ( 3.*z(3)**3 - 24.*z(3)**2 + 60.*z(3) - 44.)/6.D0
         CXX(3)= (            -  9.*z(3)**2 + 48.*z(3) - 60.)/6.D0
         z(3) = z(3) + 1.D0
         CX(4) = (   -z(3)**3 + 12.*z(3)**2 - 48.*z(3) + 64.)/6.D0
         CXX(4)= (               3.*z(3)**2 - 24.*z(3) + 48.)/6.D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  Qtmp(ix,iy,iz) = Qtmp(ix,iy,iz) + q1*AX(i)*BX(j)*CX(k)
                  dQdx(i,j,k,m,1) = &
                       & q1*AXX(i)*BX(j)*CX(k)*DBLE(param%K(1))/param%box(1)
                  dQdx(i,j,k,m,2) = &
                       & q1*AX(i)*BXX(j)*CX(k)*DBLE(param%K(2))/param%box(2)
                  dQdx(i,j,k,m,3) = &
                       & q1*AX(i)*BX(j)*CXX(k)*DBLE(param%K(3))/param%box(3)
               END DO
            END DO
         END DO
      END DO
      m = m + NS%NOPT(n,2)
   END DO
   m = 0
   q1 = param%q(2)
   DO l=2, COMM%Nlsize
      n = COMM%iy*(COMM%Ncolmn-1) + l - 1
      m = m + NS%NOPT(n,1)
      DO x=1, NS%NOPT(n,2)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         z(1) = DINT(ux) - ux + 1.D0
         z(2) = DINT(uy) - uy + 1.D0
         z(3) = DINT(uz) - uz + 1.D0
         AX(1) =      z(1)**3                                /6.D0
         AXX(1)=              -  3.*z(1)**2                  /6.D0
         z(1) = z(1) + 1.D0
         AX(2) = (-3.*z(1)**3 + 12.*z(1)**2 - 12.*z(1) +  4.)/6.D0
         AXX(2)= (               9.*z(1)**2 - 24.*z(1) + 12.)/6.D0
         z(1) = z(1) + 1.D0
         AX(3) = ( 3.*z(1)**3 - 24.*z(1)**2 + 60.*z(1) - 44.)/6.D0
         AXX(3)= (            -  9.*z(1)**2 + 48.*z(1) - 60.)/6.D0
         z(1) = z(1) + 1.D0
         AX(4) = (   -z(1)**3 + 12.*z(1)**2 - 48.*z(1) + 64.)/6.D0
         AXX(4)= (               3.*z(1)**2 - 24.*z(1) + 48.)/6.D0
         BX(1) =      z(2)**3                                /6.D0
         BXX(1)=              -  3.*z(2)**2                  /6.D0
         z(2) = z(2) + 1.D0
         BX(2) = (-3.*z(2)**3 + 12.*z(2)**2 - 12.*z(2) +  4.)/6.D0
         BXX(2)= (               9.*z(2)**2 - 24.*z(2) + 12.)/6.D0
         z(2) = z(2) + 1.D0
         BX(3) = ( 3.*z(2)**3 - 24.*z(2)**2 + 60.*z(2) - 44.)/6.D0
         BXX(3)= (            -  9.*z(2)**2 + 48.*z(2) - 60.)/6.D0
         z(2) = z(2) + 1.D0
         BX(4) = (   -z(2)**3 + 12.*z(2)**2 - 48.*z(2) + 64.)/6.D0
         BXX(4)= (               3.*z(2)**2 - 24.*z(2) + 48.)/6.D0
         CX(1) =      z(3)**3                                /6.D0
         CXX(1)=              -  3.*z(3)**2                  /6.D0
         z(3) = z(3) + 1.D0
         CX(2) = (-3.*z(3)**3 + 12.*z(3)**2 - 12.*z(3) +  4.)/6.D0
         CXX(2)= (               9.*z(3)**2 - 24.*z(3) + 12.)/6.D0
         z(3) = z(3) + 1.D0
         CX(3) = ( 3.*z(3)**3 - 24.*z(3)**2 + 60.*z(3) - 44.)/6.D0
         CXX(3)= (            -  9.*z(3)**2 + 48.*z(3) - 60.)/6.D0
         z(3) = z(3) + 1.D0
         CX(4) = (   -z(3)**3 + 12.*z(3)**2 - 48.*z(3) + 64.)/6.D0
         CXX(4)= (               3.*z(3)**2 - 24.*z(3) + 48.)/6.D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  Qtmp(ix,iy,iz) = Qtmp(ix,iy,iz) + q1*AX(i)*BX(j)*CX(k)
                  dQdx(i,j,k,m,1) = &
                       & q1*AXX(i)*BX(j)*CX(k)*DBLE(param%K(1))/param%box(1)
                  dQdx(i,j,k,m,2) = &
                       & q1*AX(i)*BXX(j)*CX(k)*DBLE(param%K(2))/param%box(2)
                  dQdx(i,j,k,m,3) = &
                       & q1*AX(i)*BX(j)*CXX(k)*DBLE(param%K(3))/param%box(3)
               END DO
            END DO
         END DO
      END DO
   END DO
   CALL MPI_ALLREDUCE(Qtmp, Qsum, param%K(1)*param%K(2)*param%K(3), &
        & MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
   Q(:,:,:) = DCMPLX(Qsum(:,:,:))
   Qinv = Q
   CALL ZFFT3D(Qinv,param%K(1),param%K(2),param%K(3),0)
   CALL ZFFT3D(Qinv,param%K(1),param%K(2),param%K(3),+1)
   !
   ! Convolution
   DO i=1, param%K(1)
      z(1) = (4.D0 + 2.D0*DCOS(2.D0*pi*(i-1)/param%K(1)))**2
      DO j=1, param%K(2)
         z(2) = (4.D0 + 2.D0*DCOS(2.D0*pi*(j-1)/param%K(2)))**2
         DO k=1, param%K(3)
            z(3) = (4.D0 + 2.D0*DCOS(2.D0*pi*(k-1)/param%K(3)))**2
            Qinv(i,j,k) = Qinv(i,j,k)*PSI(i,j,k)*36.D0**3/z(1)/z(2)/z(3)
         END DO
      END DO
   END DO
   !
   !
   Qfin = Qinv
   CALL ZFFT3D(Qfin, param%K(1), param%K(2), param%K(3), 0)
   CALL ZFFT3D(Qfin, param%K(1), param%K(2), param%K(3), -1)
   !
   ! Sum of PME energy
   Epme = 0.
   DO i=1, param%K(1)
      DO j=1, param%K(2)
         DO k=1, param%K(3)
            Epme = Epme + Q(i,j,k)*Qfin(i,j,k)
         END DO
      END DO
   END DO
   Epme = Epme*DBLE(param%K(1)*param%K(2)*param%K(3))
   q1 = eps*DBLE(param%K(1)*param%K(2)*param%K(3))
   !
   ! Force estimation
   m = 0 
   DO l=2, COMM%Nlsize
      n = COMM%iy*(COMM%Ncolmn-1) + l - 1
      DO x=1, NS%NOPT(n,1)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         pme(m)%xx(:) = 0.0D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  pme(m)%xx(1) = pme(m)%xx(1) - &
                       & q1*dQdx(i,j,k,m,1)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(2) = pme(m)%xx(2) - &
                       & q1*dQdx(i,j,k,m,2)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(3) = pme(m)%xx(3) - &
                       & q1*dQdx(i,j,k,m,3)*Qfin(ix,iy,iz)/V/pi
               END DO
            END DO
         END DO
      END DO
      m = m + NS%NOPT(n,2)
   END DO
   m = 0
   DO l=2, COMM%Nlsize
      n = COMM%iy*(COMM%Ncolmn-1) + l - 1
      m = m + NS%NOPT(n,1)
      DO x=1, NS%NOPT(n,2)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         pme(m)%xx(:) = 0.0D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  pme(m)%xx(1) = pme(m)%xx(1) - &
                       & q1*dQdx(i,j,k,m,1)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(2) = pme(m)%xx(2) - &
                       & q1*dQdx(i,j,k,m,2)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(3) = pme(m)%xx(3) - &
                       & q1*dQdx(i,j,k,m,3)*Qfin(ix,iy,iz)/V/pi
               END DO
            END DO
         END DO
      END DO
   END DO
   !
   sys%Um = eps*Epme/V/pi/2.
!   sys%Epot = sys%Epot + (sys%Ur + sys%Uo + sys%Um)*eps
   DEALLOCATE(dQdx, Qtmp, Qsum, Q, Qinv, Qfin, Qalg)

END IF
!
RETURN
!
END SUBROUTINE FORCE_SPME
!
! ################ Routine for SPME ###########################################
! ################ Using fftw3 library ###################################
!
SUBROUTINE Force_SPME4(NS, param, sys, COMM, PSI, gsend, pme)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!INCLUDE 'fftw3.f'
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg 
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PM):: param
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(GH):: gsend(NS%Ng), pme(NS%NLpt)
REAL(KIND=DP)  :: PSI(param%K(1), param%K(2), param%K(3))
!
!
INTEGER  :: i, j, k, l, m, x, ix, iy, iz, IERR
REAL(KIND=DP):: V, ux, uy, uz, Epme, AX(4), BX(4), CX(4), q1, z(3), &
     AXX(4), BXX(4), CXX(4)
REAL(KIND=DP),  ALLOCATABLE:: dQdx(:,:,:,:,:), Qtmp(:,:,:), Qsum(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE::Q(:,:,:), Qinv(:,:,:), Qfin(:,:,:), Qalg(:,:,:)
!
CALL MPI_GATHERV(gsend, NS%CNT(COMM%LID+1), COMM%PTCL, pme, NS%CNT, &
     NS%DPL, COMM%PTCL, 0, COMM%LOCAL, IERR)
!
IF (.NOT. COMM%TAG) THEN
   !
   ! Initialization and build Q matrix
   ALLOCATE(dQdx(4, 4, 4, NS%NTpt,3), &
        & Qtmp(param%K(1), param%K(2), param%K(3)), &
        & Qsum(param%K(1), param%K(2), param%K(3)), &
        & Q(param%K(1), param%K(2), param%K(3)), &
        & Qinv(param%K(1), param%K(2), param%K(3)), &
        & Qfin(param%K(1), param%K(2), param%K(3)), &
        & Qalg(param%K(1), param%K(2), param%K(3)))
   V = param%box(1)*param%box(2)*param%box(3)
   Q = 0.0
   Qtmp = 0.0
   m = 0
   q1 = param%q(1)
   DO l=1, COMM%Ncolmn-1
      DO x=1, NS%NGHB(l,1)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         z(1) = DINT(ux) - ux + 1.D0
         z(2) = DINT(uy) - uy + 1.D0
         z(3) = DINT(uz) - uz + 1.D0
         AX(1) =      z(1)**3                                /6.D0
         AXX(1)=              -  3.*z(1)**2                  /6.D0
         z(1) = z(1) + 1.D0
         AX(2) = (-3.*z(1)**3 + 12.*z(1)**2 - 12.*z(1) +  4.)/6.D0
         AXX(2)= (               9.*z(1)**2 - 24.*z(1) + 12.)/6.D0
         z(1) = z(1) + 1.D0
         AX(3) = ( 3.*z(1)**3 - 24.*z(1)**2 + 60.*z(1) - 44.)/6.D0
         AXX(3)= (            -  9.*z(1)**2 + 48.*z(1) - 60.)/6.D0
         z(1) = z(1) + 1.D0
         AX(4) = (   -z(1)**3 + 12.*z(1)**2 - 48.*z(1) + 64.)/6.D0
         AXX(4)= (               3.*z(1)**2 - 24.*z(1) + 48.)/6.D0
         BX(1) =      z(2)**3                                /6.D0
         BXX(1)=              -  3.*z(2)**2                  /6.D0
         z(2) = z(2) + 1.D0
         BX(2) = (-3.*z(2)**3 + 12.*z(2)**2 - 12.*z(2) +  4.)/6.D0
         BXX(2)= (               9.*z(2)**2 - 24.*z(2) + 12.)/6.D0
         z(2) = z(2) + 1.D0
         BX(3) = ( 3.*z(2)**3 - 24.*z(2)**2 + 60.*z(2) - 44.)/6.D0
         BXX(3)= (            -  9.*z(2)**2 + 48.*z(2) - 60.)/6.D0
         z(2) = z(2) + 1.D0
         BX(4) = (   -z(2)**3 + 12.*z(2)**2 - 48.*z(2) + 64.)/6.D0
         BXX(4)= (               3.*z(2)**2 - 24.*z(2) + 48.)/6.D0
         CX(1) =      z(3)**3                                /6.D0
         CXX(1)=              -  3.*z(3)**2                  /6.D0
         z(3) = z(3) + 1.D0
         CX(2) = (-3.*z(3)**3 + 12.*z(3)**2 - 12.*z(3) +  4.)/6.D0
         CXX(2)= (               9.*z(3)**2 - 24.*z(3) + 12.)/6.D0
         z(3) = z(3) + 1.D0
         CX(3) = ( 3.*z(3)**3 - 24.*z(3)**2 + 60.*z(3) - 44.)/6.D0
         CXX(3)= (            -  9.*z(3)**2 + 48.*z(3) - 60.)/6.D0
         z(3) = z(3) + 1.D0
         CX(4) = (   -z(3)**3 + 12.*z(3)**2 - 48.*z(3) + 64.)/6.D0
         CXX(4)= (               3.*z(3)**2 - 24.*z(3) + 48.)/6.D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  Qtmp(ix,iy,iz) = Qtmp(ix,iy,iz) + q1*AX(i)*BX(j)*CX(k)
                  dQdx(i,j,k,m,1) = &
                       & q1*AXX(i)*BX(j)*CX(k)*DBLE(param%K(1))/param%box(1)
                  dQdx(i,j,k,m,2) = &
                       & q1*AX(i)*BXX(j)*CX(k)*DBLE(param%K(2))/param%box(2)
                  dQdx(i,j,k,m,3) = &
                       & q1*AX(i)*BX(j)*CXX(k)*DBLE(param%K(3))/param%box(3)
               END DO
            END DO
         END DO
      END DO
      m = m + NS%NGHB(l,2)
   END DO
   m = 0
   q1 = param%q(2)
   DO l=1, COMM%Ncolmn-1
      m = m + NS%NGHB(l,1)
      DO x=1, NS%NGHB(l,2)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5)/param%box(3)
         z(1) = DINT(ux) - ux + 1.D0
         z(2) = DINT(uy) - uy + 1.D0
         z(3) = DINT(uz) - uz + 1.D0
         AX(1) =      z(1)**3                                /6.D0
         AXX(1)=              -  3.*z(1)**2                  /6.D0
         z(1) = z(1) + 1.D0
         AX(2) = (-3.*z(1)**3 + 12.*z(1)**2 - 12.*z(1) +  4.)/6.D0
         AXX(2)= (               9.*z(1)**2 - 24.*z(1) + 12.)/6.D0
         z(1) = z(1) + 1.D0
         AX(3) = ( 3.*z(1)**3 - 24.*z(1)**2 + 60.*z(1) - 44.)/6.D0
         AXX(3)= (            -  9.*z(1)**2 + 48.*z(1) - 60.)/6.D0
         z(1) = z(1) + 1.D0
         AX(4) = (   -z(1)**3 + 12.*z(1)**2 - 48.*z(1) + 64.)/6.D0
         AXX(4)= (               3.*z(1)**2 - 24.*z(1) + 48.)/6.D0
         BX(1) =      z(2)**3                                /6.D0
         BXX(1)=              -  3.*z(2)**2                  /6.D0
         z(2) = z(2) + 1.D0
         BX(2) = (-3.*z(2)**3 + 12.*z(2)**2 - 12.*z(2) +  4.)/6.D0
         BXX(2)= (               9.*z(2)**2 - 24.*z(2) + 12.)/6.D0
         z(2) = z(2) + 1.D0
         BX(3) = ( 3.*z(2)**3 - 24.*z(2)**2 + 60.*z(2) - 44.)/6.D0
         BXX(3)= (            -  9.*z(2)**2 + 48.*z(2) - 60.)/6.D0
         z(2) = z(2) + 1.D0
         BX(4) = (   -z(2)**3 + 12.*z(2)**2 - 48.*z(2) + 64.)/6.D0
         BXX(4)= (               3.*z(2)**2 - 24.*z(2) + 48.)/6.D0
         CX(1) =      z(3)**3                                /6.D0
         CXX(1)=              -  3.*z(3)**2                  /6.D0
         z(3) = z(3) + 1.D0
         CX(2) = (-3.*z(3)**3 + 12.*z(3)**2 - 12.*z(3) +  4.)/6.D0
         CXX(2)= (               9.*z(3)**2 - 24.*z(3) + 12.)/6.D0
         z(3) = z(3) + 1.D0
         CX(3) = ( 3.*z(3)**3 - 24.*z(3)**2 + 60.*z(3) - 44.)/6.D0
         CXX(3)= (            -  9.*z(3)**2 + 48.*z(3) - 60.)/6.D0
         z(3) = z(3) + 1.D0
         CX(4) = (   -z(3)**3 + 12.*z(3)**2 - 48.*z(3) + 64.)/6.D0
         CXX(4)= (               3.*z(3)**2 - 24.*z(3) + 48.)/6.D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  Qtmp(ix,iy,iz) = Qtmp(ix,iy,iz) + q1*AX(i)*BX(j)*CX(k)
                  dQdx(i,j,k,m,1) = &
                       & q1*AXX(i)*BX(j)*CX(k)*DBLE(param%K(1))/param%box(1)
                  dQdx(i,j,k,m,2) = &
                       & q1*AX(i)*BXX(j)*CX(k)*DBLE(param%K(2))/param%box(2)
                  dQdx(i,j,k,m,3) = &
                       & q1*AX(i)*BX(j)*CXX(k)*DBLE(param%K(3))/param%box(3)
               END DO
            END DO
         END DO
      END DO
   END DO
   CALL MPI_ALLREDUCE(Qtmp, Qsum, param%K(1)*param%K(2)*param%K(3), &
        & MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
   Q(:,:,:) = DCMPLX(Qsum(:,:,:))
   Qinv = Q
!   CALL dfftw_plan_dft_3d(plan,param%K(1),param%K(2),param%K(3), &
!        Q, Qinv, FFTW_BACKWARD, FFTW_ESTIMATE)
!   CALL dfftw_execute(plan)
!   CALL dfftw_destroy_plan(plan)
   CALL ZFFT3D(Qinv,param%K(1),param%K(2),param%K(3),0)
   CALL ZFFT3D(Qinv,param%K(1),param%K(2),param%K(3),+1)
   !
   ! Convolution
   DO i=1, param%K(1)
      z(1) = (4.D0 + 2.D0*DCOS(2.D0*pi*(i-1)/param%K(1)))**2
      DO j=1, param%K(2)
         z(2) = (4.D0 + 2.D0*DCOS(2.D0*pi*(j-1)/param%K(2)))**2
         DO k=1, param%K(3)
            z(3) = (4.D0 + 2.D0*DCOS(2.D0*pi*(k-1)/param%K(3)))**2
            Qinv(i,j,k) = Qinv(i,j,k)*PSI(i,j,k)*36.D0**3/z(1)/z(2)/z(3)
         END DO
      END DO
   END DO
   !
   !
   Qfin = Qinv
!   CALL dfftw_plan_dft_3d(plan, param%K(1), param%K(2), param%K(3), &
!        Qinv, Qfin, FFTW_FORWARD, FFTW_ESTIMATE)
!   CALL dfftw_execute(plan)
!   CALL dfftw_destroy_plan(plan)
   CALL ZFFT3D(Qfin, param%K(1), param%K(2), param%K(3), 0)
   CALL ZFFT3D(Qfin, param%K(1), param%K(2), param%K(3), -1)
   !
   ! Sum of PME energy
   Epme = 0.D0
   DO i=1, param%K(1)
      DO j=1, param%K(2)
         DO k=1, param%K(3)
            Epme = Epme + Q(i,j,k)*Qfin(i,j,k)
         END DO
      END DO
   END DO
   !Epme = Epme !*DBLE(param%K(1)*param%K(2)*param%K(3))
   !q1 = eps !*DBLE(param%K(1)*param%K(2)*param%K(3))
   Epme = Epme*DBLE(param%K(1)*param%K(2)*param%K(3))
   q1 = eps*DBLE(param%K(1)*param%K(2)*param%K(3))
   !
   ! Force estimation
   m = 0 
   DO l=1, COMM%Ncolmn-1
      DO x=1, NS%NGHB(l,1)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5D0)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5D0)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5D0)/param%box(3)
         pme(m)%xx(:) = 0.0D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  pme(m)%xx(1) = pme(m)%xx(1) - &
                       & q1*dQdx(i,j,k,m,1)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(2) = pme(m)%xx(2) - &
                       & q1*dQdx(i,j,k,m,2)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(3) = pme(m)%xx(3) - &
                       & q1*dQdx(i,j,k,m,3)*Qfin(ix,iy,iz)/V/pi
               END DO
            END DO
         END DO
      END DO
      m = m + NS%NGHB(l,2)
   END DO
   m = 0
   DO l=1, COMM%Ncolmn-1
      m = m + NS%NGHB(l,1)
      DO x=1, NS%NGHB(l,2)
         m = m + 1
         ux = DBLE(param%K(1))*(pme(m)%xx(1)+param%box(1)*.5D0)/param%box(1)
         uy = DBLE(param%K(2))*(pme(m)%xx(2)+param%box(2)*.5D0)/param%box(2)
         uz = DBLE(param%K(3))*(pme(m)%xx(3)+param%box(3)*.5D0)/param%box(3)
         pme(m)%xx(:) = 0.0D0
         DO i=1,4
            DO j=1,4
               DO k=1,4
                  ix = MOD(INT(ux)+i-2+param%K(1),param%K(1))+1
                  iy = MOD(INT(uy)+j-2+param%K(2),param%K(2))+1
                  iz = MOD(INT(uz)+k-2+param%K(3),param%K(3))+1
                  pme(m)%xx(1) = pme(m)%xx(1) - &
                       & q1*dQdx(i,j,k,m,1)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(2) = pme(m)%xx(2) - &
                       & q1*dQdx(i,j,k,m,2)*Qfin(ix,iy,iz)/V/pi
                  pme(m)%xx(3) = pme(m)%xx(3) - &
                       & q1*dQdx(i,j,k,m,3)*Qfin(ix,iy,iz)/V/pi
               END DO
            END DO
         END DO
      END DO
   END DO
   !
   sys%Um = eps*Epme/V/pi/2.D0
!   sys%Epot = sys%Epot + (sys%Ur + sys%Uo + sys%Um)*eps
   DEALLOCATE(dQdx, Qtmp, Qsum, Q, Qinv, Qfin, Qalg)

END IF
!
RETURN
!
END SUBROUTINE FORCE_SPME4
!
! Scatter PME force onto direct PEs ###########################################
SUBROUTINE SPME_SUM(NS, qel, qion, sys, COMM, pme, grec)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(GH):: grec(NS%Ng), pme(NS%NLpt)
!
!
INTEGER :: i, k, IERR
!
CALL MPI_BCAST(sys%Um, 1, MPI_REAL8, 0, COMM%LOCAL, IERR)
CALL MPI_SCATTERV(pme, NS%CNT, NS%DPL, COMM%PTCL, grec, NS%CNT(COMM%LID+1), &
     & COMM%PTCL, 0, COMM%LOCAL, IERR)
!
IF (COMM%TAG) THEN
   !
   ! Electron force sum
   DO i=1, NS%Nel
      DO k=1, 3
         qel(i)%ff(k) = qel(i)%ff(k) + grec(i)%xx(k)
      END DO
   END DO
   !
   ! Ion force sum
   DO i=1, NS%Nion
      DO k=1, 3
         qion(i)%ff(k) = qion(i)%ff(k) + grec(i+NS%Nel)%xx(k)
      END DO
   END DO
   !
END IF
!
RETURN 
END SUBROUTINE SPME_SUM

!
! Remove rigid body motion ####################################################
! #############################################################################
SUBROUTINE REMOVE_RIGIDMTN(NS, qel, qion, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(MP):: COMM
!!
INTEGER:: i, IERR
REAL(KIND=DP) :: Fx(3), ifx(3)
!
Fx = 0.0
DO i=1,NS%Nel
   Fx(:) = Fx(:) + qel(i)%ff(:)
END DO
DO i=1, NS%Nion
   Fx(:) = Fx(:) + qion(i)%ff(:)
END DO
CALL MPI_ALLREDUCE(Fx, ifx, 3, MPI_REAL8, MPI_SUM, COMM%COMM1D, IERR)
!
! Find average rigid force
ifx(:) = ifx(:)/DBLE(NS%NTpt)
!
! Subtraction for all of the particles
DO i=1,NS%Nel
   qel(i)%ff(:) = qel(i)%ff(:) - ifx(:)
END DO
DO i=1, NS%Nion
   qion(i)%ff(:) = qion(i)%ff(:) - ifx(:)
END DO
!
RETURN
END SUBROUTINE REMOVE_RIGIDMTN
!
! Pressure estimation  ########################################################
! #############################################################################
SUBROUTINE BAROSTAT(NS, qel, qion, sys, param)
USE DATASTR
IMPLICIT NONE
!
!
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(PM):: param
TYPE(ST):: sys
!
!
INTEGER:: i, j
REAL(KIND=DP) :: Press, V
!
! Electron pressure
V = param%box(1)*param%box(2)*param%box(3)
Press = 0.0
DO i=1, NS%Nel
   DO j=1,3
      Press = Press + qel(i)%xx(j)*qel(i)%ff(j)
   END DO
END DO
sys%Pel = Press
!
! Ion pressure
Press = 0.0
DO i=1, NS%Nion
   DO j=1,3
      Press = Press + qion(i)%xx(j)*qion(i)%ff(j)
   END DO
END DO
sys%Pion = Press
!
RETURN
END SUBROUTINE BAROSTAT
!

