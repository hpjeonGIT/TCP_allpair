!
! ############ binary restart file generation ###################
SUBROUTINE Restart(NS, qel, qion, ttime, param, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(TM):: ttime
TYPE(PM):: param
TYPE(MP):: COMM
!
!
INTEGER:: i, j, k, Nfreq, ccnt(COMM%NUMPROC), displ(COMM%NUMPROC), IERR
INTEGER:: eid(NS%Nel), iid(NS%Nion), id(NS%NTpt)
REAL(KIND=DP) :: ixx(NS%Nion*3), ixv(NS%Nion*3), iff(NS%Nion*3), &
     & exx(NS%Nel*3), exv(NS%Nel*3), eff(NS%Nel*3), &
     & xx(NS%NTpt*3), xv(NS%NTpt*3), ff(NS%NTpt*3)
CHARACTER(LEN=20):: FILENAME, FILE2
!
! Dump position and velocity into local variables
k = 0
DO i=1, NS%Nel
   DO j=1, 3
      k = k + 1
      exx(k) = qel(i)%xx(j)
      exv(k) = qel(i)%xv(j)
      eff(k) = qel(i)%ff(j)
   END DO
   eid(i) = qel(i)%id
END DO
k = 0
DO i=1, NS%Nion
   DO j=1, 3
      k = k + 1
      ixx(k) = qion(i)%xx(j)
      ixv(k) = qion(i)%xv(j)
      iff(k) = qion(i)%ff(j)
   END DO
   iid(i) = qion(i)%id
END DO
!
! Electron data dump
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,1)*3
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(exx, NS%Nel*3, MPI_REAL8, xx, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_GatherV(exv, NS%Nel*3, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_GatherV(eff, NS%Nel*3, MPI_REAL8, ff, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
!
! id set
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,1)
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(eid, NS%Nel, MPI_INTEGER, id, ccnt, displ, MPI_INTEGER, 0, &
     & COMM%COMM1D, IERR)
IF (COMM%MYID2 == 0) THEN
   !
   ! File name decision
   Nfreq = INT(ttime%Nloop/ttime%Nrest)
   WRITE(FILENAME,50) Nfreq
   WRITE(FILE2,   60) Nfreq
   !
   OPEN(UNIT=30,file=FILENAME)
   WRITE(30,*) NS%NTel+ NS%NTion
   WRITE(30,*) "frame = ", Nfreq, "energy = 0"
   k = 1
   DO i=1, NS%NTel
      WRITE(30,200) "Cu ", xx(k), xx(k+1), xx(k+2), xv(k), xv(k+1), xv(k+2), ff(k), ff(k+1), ff(k+2), id(i)
      k = k + 3
   END DO
END IF
!
! Ion data dump
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,2)*3
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(ixx, NS%Nion*3, MPI_REAL8, xx, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_GatherV(ixv, NS%Nion*3, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_GatherV(iff, NS%Nion*3, MPI_REAL8, ff, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
!
! id set
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,2)
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(iid, NS%Nion, MPI_INTEGER, id, ccnt, displ, MPI_INTEGER, 0, &
     & COMM%COMM1D, IERR)


IF (COMM%MYID2 == 0) THEN
   k = 1
   DO i=1, NS%NTion
      WRITE(30,200) "C ", xx(k), xx(k+1), xx(k+2), xv(k), xv(k+1), xv(k+2), ff(k), ff(k+1), ff(k+2), id(i)
     k = k + 3
   END DO
   CLOSE(30)
   !
   ! Sampling file close and reopen
   CLOSE(55)
   OPEN(UNIT=55,file=FILE2)
   WRITE(55,100)
END IF
50  FORMAT("REST", I2.2, ".xyz")
60  FORMAT("ener", I2.2, ".dat")
100 FORMAT("# time(fs),       T_el(eV),      T_ion(eV), Potential Energy, &
         & Total energy, Electron VAF, Ion VAF,  mutual VAF")
200 FORMAT(A3, 9(1X, ES14.6), 1X, I6)
!
RETURN
END SUBROUTINE Restart
!
!
!
! ############ binary initial velocity file generation ###################
SUBROUTINE  Vinit_write(NS, qel, qion, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(MP):: COMM
!
!
INTEGER:: i, j, k, ccnt(COMM%NUMPROC), displ(COMM%NUMPROC), IERR
REAL(KIND=DP) :: ixv(NS%Nion*3), exv(NS%Nel*3), xv(NS%NTpt*3)
!
! Electron data dump
k = 0
DO i=1, NS%Nel
   DO j=1, 3
      k = k + 1
      exv(k) = qel(i)%xv(j)
   END DO
END DO
k = 0
DO i=1, NS%Nion
   DO j=1, 3
      k = k + 1
      ixv(k) = qion(i)%xv(j)
   END DO
END DO
!
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,1)*3
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(exv, NS%Nel*3, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
IF (COMM%MYID2 == 0) THEN
   OPEN(UNIT=30,file="init_v.dat")
   k = 1
   DO i=1, NS%NTel
      WRITE(30,*) xv(k), xv(k+1), xv(k+2)
      k = k + 3
   END DO
END IF
!
! Ion data dump
DO i=1, COMM%NUMPROC
   ccnt(i) = NS%NOPT(i,2)*3
END DO
displ(1) = 0
DO i=1, COMM%NUMPROC-1
   displ(i+1) = ccnt(i) + displ(i)
END DO
CALL MPI_GatherV(ixv, NS%Nion*3, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     & COMM%COMM1D, IERR)
IF (COMM%MYID2 == 0) THEN
   k = 1
   DO i=1, NS%NTion
      WRITE(30,*) xv(k), xv(k+1), xv(k+2)
     k = k + 3
   END DO
   CLOSE(30)
END IF
RETURN
END SUBROUTINE Vinit_write
!
! ################### RDF print ####################
SUBROUTINE DUMP(NS, ttime, param, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(TM):: ttime
TYPE(PM):: param
TYPE(MP):: COMM
!
!
INTEGER:: i, Nfreq, IERR
REAL(KIND=DP) :: Ntau, V, dV_ee, dV_ei, dV_ii, &
     & Nee(Ndist+1), Nei(Ndist+1), Nii(Ndist+1), VDe(Ndist+1), VDi(Ndist+1)
CHARACTER(LEN=20):: FILERDF, FILEVDF
!
CALL MPI_REDUCE(NS%RDF_ee, Nee, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%RDF_ei, Nei, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%RDF_ii, Nii, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%VDF_e, VDe, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     & COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%VDF_i, VDi, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     & COMM%COMM1D, IERR)
!
IF (COMM%MYID2 == 0) THEN
   Ntau = DBLE(ttime%Ndump)
   Nfreq = INT(ttime%Nloop/ttime%Ndump)
   WRITE(FILERDF,20) Nfreq
   WRITE(FILEVDF,30) Nfreq
   !
   OPEN(UNIT=30,file=FILERDF)
   V = param%box(1)*param%box(2)*param%box(3)
   WRITE(30,100)
   dV_ee=DBLE(NS%NTel)**2 *Ntau*NS%dx_ee**3*4.*pi/3.D0/V
   dV_ei=DBLE(NS%NTpt)**2 *Ntau*NS%dx_ei**3*4.*pi/3.D0/V
   dV_ii=DBLE(NS%NTion)**2*Ntau*NS%dx_ii**3*4.*pi/3.D0/V   
   DO i=1, Ndist
      WRITE(30,200) DBLE(i)*NS%dx_ee, Nee(i)/dV_ee/DBLE(i**3-(i-1)**3), &
           & DBLE(i)*NS%dx_ei, Nei(i)/dV_ei/DBLE(i**3-(i-1)**3), &
           & DBLE(i)*NS%dx_ii, Nii(i)/dV_ii/DBLE(i**3-(i-1)**3)
   END DO
   CLOSE(30)
   !
   ! ==================== VDF printout #####################################
   OPEN(UNIT=31, file = FILEVDF)
   WRITE(31,300)
   DO i=1, Ndist
      WRITE(31, 400) NS%dv_e*DBLE(i), VDe(i)/Ntau, NS%dv_i*DBLE(i), VDi(i)/Ntau
   END DO
   CLOSE(31)
END IF
NS%RDF_ee = 0.D0
NS%RDF_ei = 0.D0
NS%RDF_ii = 0.D0
NS%VDF_e  = 0.D0
NS%VDF_i  = 0.D0
!
20  FORMAT("RDF", I3.3, ".dat")
30  FORMAT("VDF", I3.3, ".dat")
100 FORMAT("#  x_ee,      RDF_ee,      x_ei    RDF_ei,       x_ii     RDF_ii")
200 FORMAT(6(ES14.6, 1X))
300 FORMAT("#  T_e        VDF_e            T_i      VDF_i")
400 FORMAT(4(ES14.6, 1X))
!
RETURN
END SUBROUTINE DUMP
!
! ################ Intermediate temperature and energy print #################
SUBROUTINE SAMPLING(NS, ttime, sys, COMM, param)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(TM):: ttime
TYPE(ST):: sys
TYPE(MP):: COMM
TYPE(PM):: param
!
!
INTEGER:: IERR
REAL(KIND=DP) :: Rel, Rion, mv2el, mv2ion, Epot, cel, cion, V, Jei(3), cei
!
! Reference number
V = param%box(1)*param%box(2)*param%box(3)
Rel = DBLE(NS%NTel*3 - Nref)
Rion = DBLE(NS%NTion*3 - Nref)
!
! Collect data
CALL MPI_REDUCE(sys%mv2el,  mv2el, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%mv2ion, mv2ion,1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%Epot,   Epot,  1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
!CALL MPI_REDUCE(sys%cel,    cel,   1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
!CALL MPI_REDUCE(sys%cion,   cion,  1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
!CALL MPI_REDUCE(sys%Jei,    Jei,   3, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
!
! Print data
!IF (COMM%MYID2 == 0) PRINT *, "direct + const", Epot
Epot = sys%Um + Epot
!cei  = (Jei(1)*sys%J0(1) + Jei(2)*sys%J0(2) + Jei(3)*sys%J0(3)) / &
!     & (sys%J0(1)**2 + sys%J0(2)**2 + sys%J0(3)**2)
IF (COMM%MYID2 == 0) THEN
!   PRINT *, Epot-sys%Um, "direct", sys%Um, "pme"
   cel = cel/DBLE(NS%NTel)
   cion = cion/DBLE(NS%NTion)
   WRITE(55, 200) ttime%tnow*tps, mv2el/Rel, mv2ion/Rion, Epot, &
        & 0.5D0*(mv2el+mv2ion)+Epot
!, cel, cion
   !& cei*(1./param%xm(1) + 1./param%xm(2))*0.5D0
END IF
200 FORMAT(5(ES15.8, 1X))
!
RETURN
END SUBROUTINE SAMPLING
