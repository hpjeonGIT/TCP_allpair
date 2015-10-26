PROGRAM TCP_SPME_VAF
!
! Two component plasma - electron and ion - analysis code for PBC
! This code is specified for very high density plasma
!
! Aug. 03. 2006,
! Theoretical division, Los Alamos National Laboratory
! Department of Applied Science, University of California, Davis
! Byoungseon Jeon
!
! Sept. 11. 2006,
! Smooth particle mesh Ewald imlemented
!
! Sept. 16. 2006,
! FFTE library applied
!
! Oct. 03. 2006,
! DPME implemented - split PME and direct sum nodes
!
! Oct. 30. 2006,
! Diffusion coefficient and particle tracing implemented.
! Position array %xx save position as absolute ones - not trucated by PBC
!
! Dec. 20. 2006.
! VAF implemented. Diffusion coefficient canceled.
!
! Feb. 20. 2007.
! Implement 5th order interpolation for SPME
!
! July 21. 2009.
! Correct double precision functions and force routine at the initial step
!
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
!
TYPE(NM):: NS
TYPE(PT), POINTER:: qel(:), qion(:)
REAL(KIND=DP) ,  POINTER:: PSI(:,:,:)
TYPE(GH), ALLOCATABLE:: gsend(:), gemit(:), grec(:), pme(:)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(TM):: ttime
TYPE(MP):: COMM
REAL(KIND=DP):: dt
!
!
! INTERNAL VARIABLES
REAL:: time0, time1, time2, secnds
REAL(KIND=DP):: a_time, r1_time,r2_time, d_time,o_time
REAL(KIND=DP):: Ur
INTEGER:: Nloop_max, NUMITER, i, Ne, Ni, ISTATUS(MPI_STATUS_SIZE), IERR
LOGICAL:: TAG, TAG_REF
time1 = secnds(time0)
!
! Main node decision
CALL MPI_INITIALIZE(COMM, NS)
IF (COMM%GID == 1) THEN
   TAG = .TRUE.
ELSE
   TAG = .FALSE.
END IF
!
! Open XYZ file for output
IF (TAG) THEN
   OPEN(UNIT=55, file="ener00.dat")
   WRITE(55,100)
END IF
100 FORMAT("# time(fs),       T_el(eV),      T_ion(eV), Potential Energy, &
         & Total energy, Electron VAF, Ion VAF,  mutual VAF")
!
! Parsing input data and allocate pointer variables
CALL Init(NS, qel, qion, ttime, param, sys, dt, PSI, COMM)
ALLOCATE(gsend(NS%Ng), gemit(NS%Ng), grec(NS%Ng), pme(NS%NLpt))
CALL EWALD_SELF(NS, param, sys, PSI, COMM)
!
! Velocity initialization
IF (param%rest == 'OFF') THEN
   ! Initial velocity from Normal distribution
   IF (COMM%TAG) CALL Vinit(NS, qel, qion, param, sys, COMM)
   CALL Finit(NS, qel, qion, param, sys, COMM, gsend, gemit, grec, pme, PSI, &
        & dt)
END IF
!
! Initialization
ttime%Ndump = NINT(ttime%tdump/dt)
ttime%Nrest = NINT(ttime%trest/dt)
ttime%Nsamp = NINT(ttime%tsamp/dt)
ttime%Nloop = 0
ttime%tnow = 0.0D0
Nloop_max = NINT(ttime%tmax/dt)
NS%RDF_ee = 0.D0
NS%RDF_ei = 0.D0
NS%RDF_ii = 0.D0
NS%VDF_e  = 0.D0
NS%VDF_i  = 0.D0
TAG_REF = .FALSE.
r1_time = 0.D0
r2_time = 0.D0
d_time = 0.D0
o_time = 0.D0
!
!
IF (COMM%ITAG) THEN 
   NUMITER = INT(COMM%NUMPROC/2) - 1
ELSE 
   NUMITER = INT(COMM%NUMPROC/2)
END IF

DO WHILE (ttime%Nloop < Nloop_max)
   a_time = MPI_Wtime()
   ttime%tnow = ttime%tnow + dt
   ttime%Nloop = ttime%Nloop + 1
   !
   ! L.2 Time integration 
   IF (COMM%TAG) THEN
      IF  (param%thermo == 'STOCH ')  THEN   
         CALL VVerletStoch1(NS, qel, qion, param, dt, gsend)
      ELSE 
         CALL VVerletNotemp1(NS, qel, qion, param, dt, gsend)
      END IF
   ELSE
      gsend(1)%xx(:) = 0.0D0
   END IF
   o_time = o_time + MPI_Wtime() - a_time
   !
   ! Allocate particle information onto PME PE
   a_time = MPI_Wtime()
   CALL Force_SPME4(NS, param, sys, COMM, PSI, gsend, pme)
   r1_time = r1_time + MPI_Wtime() - a_time
   !
   ! Direct sum only on direct PEs
   a_time = MPI_Wtime()
   IF (COMM%TAG) THEN
      CALL Force_direct_self(NS, qel, qion, param, sys)
      DO i=1, NUMITER
         Ne  = NS%NOPT(COMM%SRC(i)+1,1)
         Ni = NS%NOPT(COMM%SRC(i)+1,2)
         CALL MPI_SENDRECV(gsend, NS%Npt, COMM%PTCL, COMM%DEST(i), 1, &
              & grec, Ne+Ni, COMM%PTCL, COMM%SRC(i), 1, &
              & COMM%COMM1D, ISTATUS, IERR)
         CALL Force_direct_ext(NS, qel, qion, param, grec, gemit, &
              & Ur, Ne, Ni)
         CALL MPI_SENDRECV(gemit, Ne+Ni, COMM%PTCL, COMM%SRC(i), 0, &
              & grec, NS%Npt, COMM%PTCL, COMM%DEST(i), 0, &
              & COMM%COMM1D, ISTATUS, IERR)
         CALL Force_direct_sum(NS, qel, qion, grec)
         sys%Ur = sys%Ur + Ur
      END DO
      IF (COMM%ITAG) THEN
         Ne  = NS%NOPT(COMM%SRC(i)+1,1)
         Ni = NS%NOPT(COMM%SRC(i)+1,2)
         CALL MPI_SENDRECV(gsend, NS%Npt, COMM%PTCL, COMM%DEST(i), 1, &
              & grec, Ne+Ni, COMM%PTCL, COMM%SRC(i), 1, &
              & COMM%COMM1D, ISTATUS, IERR)
         CALL Force_direct_ext(NS, qel, qion, param, grec, gemit, &
              & Ur, Ne, Ni)
         sys%Ur = sys%Ur + Ur*0.5D0
      END IF
   END IF
   d_time = d_time + MPI_Wtime() - a_time
   !
   ! PME force scattering on direct PEs
   a_time = MPI_Wtime()
   CALL SPME_SUM(NS, qel, qion, sys, COMM, pme, grec)
   r2_time = r2_time + MPI_Wtime() - a_time
   sys%Epot = (sys%Ur + sys%Uo)*eps
   !
   ! Remove rigid body motion for direct PEs
   IF (COMM%TAG) CALL REMOVE_RIGIDMTN(NS, qel, qion, COMM)
   !
   ! Verlet routine for second step
   a_time = MPI_Wtime()
   IF  (param%thermo == 'STOCH ')  THEN   
      IF (COMM%TAG) CALL VVerletStoch2(NS, qel, qion, param, sys, dt)
   ELSE IF (param%thermo == 'BEREND')  THEN
      IF (COMM%TAG) CALL VVerletBerend2(NS, qel, qion, param, sys, dt, COMM)
   ELSE IF (param%thermo == 'ISOKIN')  THEN
      IF (COMM%TAG) CALL VVerletIsokin2(NS, qel, qion, param, sys, dt, COMM)
   ELSE IF (param%thermo == 'NOTEMP')  THEN
      IF (COMM%TAG) CALL VVerletNotemp2(NS, qel, qion, param, sys, dt)
   END IF
   !
   IF (MOD(ttime%Nloop,ttime%Ndump) == 0) THEN
      IF (COMM%TAG) CALL DUMP(NS, ttime, param, COMM)
   END IF
   IF (MOD(ttime%Nloop,ttime%Nsamp) == 0) THEN
      IF (COMM%TAG) THEN
         !CALL BAROSTAT(NS, qel, qion, sys, param)
         !CALL DIFFUSION(NS, qel, qion)
         !CALL VAF(NS, qel, qion, sys)
         CALL SAMPLING(NS, ttime, sys, COMM, param)
      END IF
   END IF
   IF (MOD(ttime%Nloop,ttime%Nrest) == 0) THEN
      IF (COMM%TAG) CALL Restart(NS, qel, qion, ttime, param, COMM)
   END IF
   o_time = o_time + MPI_Wtime() - a_time
END DO
!
IF (TAG) THEN
   CLOSE(55)
END IF

IF (COMM%TAG) THEN
   DEALLOCATE(qel, qion, NS%NOPT, COMM%DEST, COMM%SRC, gsend, gemit, grec, &
        & NS%DPL, NS%CNT, NS%NGHB, PSI, pme)
ELSE 
   DEALLOCATE(NS%NOPT, COMM%DEST, COMM%SRC, gsend, gemit, grec, PSI, &
        & NS%DPL, NS%CNT, NS%NGHB, pme)
END IF
time2 = secnds(time1)
IF (TAG) PRINT '(A15, F8.2)', "Wall time is ", time2
IF (TAG) PRINT '(A20,F8.2,1X,A16,1X,F8.2,A14,F8.2)', "direct sum cost = ", &
     & d_time, " idling time = ", o_time, "comm. time = ", r1_time + r2_time
IF (COMM%GID==0) PRINT '(A10, F8.2, A14, F8.2)', "PME cost =", r1_time, &
     & "idling time = ", r2_time
CALL MPI_Finalize(COMM%rc)
!
STOP
!
CONTAINS
!
! ################### Data parsing and memory allocation ####################
! Unit normalization
! Length: 1. means 1Angstrom = 10E-10 m
! Mass: 1. means 1.6605E-27 kg (a.m.u)
! Energy: 1. means 1 eV = 1.6022E-19 J
! Time: 1. means 10.18 fs = 1.018E-14 sec
SUBROUTINE INIT(NS, qel, qion, ttime, param, sys, dt, PSI, COMM)
USE DATASTR
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(PT), POINTER:: qel(:), qion(:)
REAL(KIND=DP) ,  POINTER:: PSI(:,:,:)
TYPE(TM):: ttime
TYPE(PM):: param
TYPE(ST):: sys
TYPE(NM):: NS
TYPE(MP):: COMM
REAL(KIND=DP)  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, Ntemp1, Ntemp2, NCPUS, OpenStatus, id
CHARACTER(len=5):: dummy
REAL(KIND=DP) :: x, y, z, vx, vy, vz, fx, fy, fz, ix, iy, iz, xm
!
!
![[[[[[[[[[[[[[[[[ "input.dat" parsing and data arrangement ]]]]]]]]]]]]]]]]]
OPEN (UNIT=11, file="tcp.prm", IOSTAT= OpenStatus)
IF (OpenStatus > 0) THEN
   CALL MPI_Finalize(COMM%rc)
   STOP " *** CANNOT OPEN tcp.prm *** "
END IF

READ(11,*) dummy
!
! Time parameter
! 1.0 = 10.18fs
READ(11,*) dummy
READ(11,*) ttime%tmax, ttime%trest, ttime%tdump, ttime%tsamp, dt
ttime%tmax  = ttime%tmax  / tps
ttime%trest = ttime%trest / tps
ttime%tdump = ttime%tdump / tps
ttime%tsamp = ttime%tsamp / tps
dt = dt / tps
!
! Box size
READ(11,*) dummy
READ(11,*) (param%box(j), j=1,3)
!
! Initial temperature
READ(11,*) dummy
READ(11,*) param%Te
READ(11,*) param%Ti
!
! i.1 Number of all particles
READ(11,*) dummy
READ(11,*) NS%NTel
READ(11,*) NS%NTion
NS%NTpt = NS%NTel+NS%NTion
!
! cut-off radius for direct sum
READ(11,*) dummy
READ(11,*) param%rc
!
! De broglie wave length
READ(11,*) dummy
READ(11,*) param%rs
!
! p.5 Mass for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%xm(i)
END DO
!
! p.5 Charge for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%q(i)
END DO
!
! p.7 Berendsen thermostat/damping constant
READ(11,*) dummy
READ(11,*) param%thermo
READ(11,*) param%tau, param%alpha
!
! Restart option
READ(11,*) dummy
READ(11,*) param%rest, param%trace
!
! RDF parameter
READ(11,*) dummy
READ(11,*) NS%dx_ee, NS%dx_ei, NS%dx_ii
NS%dv_e = param%Te*4.D0/DBLE(Ndist)
NS%dv_i = param%Ti*4.D0/DBLE(Ndist)
!
! FFT parameter
READ(11,*) dummy
READ(11,*) param%K(1), param%K(2), param%K(3)
!
! File close
CLOSE(11)
!
!#################### READ particle position and id data ######################
!
! Number of particle distribution
NCPUS = NS%NCPUS - NS%Npe
NS%NOPT(:,1) = INT(NS%NTel/NCPUS)
Ntemp1 = MOD(NS%NTel, NCPUS)
DO i=1, Ntemp1
   NS%NOPT(i,1) = NS%NOPT(i,1) + 1
END DO
NS%NOPT(:,2) = INT(NS%NTion/NCPUS)
Ntemp1 = MOD(NS%NTion, NCPUS)
DO i=1, Ntemp1
   NS%NOPT(i,2) = NS%NOPT(i,2) + 1
END DO
!
IF (COMM%TAG) THEN
   NS%Nel  = NS%NOPT(COMM%MYID2 + 1, 1)
   NS%Nion = NS%NOPT(COMM%MYID2 + 1, 2)
   NS%Npt  = NS%Nel + NS%Nion
   NS%Ng   = MAXVAL(NS%NOPT(:,1))+MAXVAL(NS%NOPT(:,2))
ELSE
   NS%Nel  = 1
   NS%Nion = 1
   NS%Npt  = 1
   NS%Ng   = NS%NTpt
END IF
NS%CNT(1) = 0
DO i=2, COMM%Nlsize
   j = COMM%iy*(COMM%Ncolmn-1) + i - 1
   NS%CNT(i) = NS%NOPT(j,1)+NS%NOPT(j,2)
   NS%NGHB(i-1,1) = NS%NOPT(j,1)
   NS%NGHB(i-1,2) = NS%NOPT(j,2)
END DO
NS%DPL = 0
DO i=1, COMM%Nlsize-1
   NS%DPL(i+1) = NS%CNT(i) + NS%DPL(i)
END DO
NS%NLpt = NS%DPL(COMM%Nlsize) + NS%CNT(COMM%Nlsize)
!
ALLOCATE(qel(NS%Nel), qion(NS%Nion), PSI(param%K(1),param%K(2),param%K(3)))
!
sys%mv2el = 0.0D0
sys%mv2ion = 0.0D0
IF (COMM%TAG) THEN
   Ntemp1 = 0
   Ntemp2 = 0
   DO i=1, COMM%MYID2
      Ntemp1 = Ntemp1 + NS%NOPT(i,1)
      Ntemp2 = Ntemp2 + NS%NOPT(i,2)
   END DO
   IF (param%rest == 'OFF') THEN
      OPEN (UNIT=15, file="input.xyz", STATUS = "OLD", IOSTAT= OpenStatus)
      IF (OpenStatus > 0) THEN
         CALL MPI_Finalize(COMM%rc)
         STOP " *** CANNOT OPEN input.xyz *** "
      END IF
      
      READ(15,*) dummy
      READ(15,*) dummy
      !
      ! Read only positions
      j = 0
      DO i=1, NS%NTel      
         READ(15,*) dummy, x, y, z, id
         IF (i > Ntemp1 .AND. j < NS%Nel) THEN
            j = j + 1
            qel(j)%xx(1) = x
            qel(j)%xx(2) = y
            qel(j)%xx(3) = z
            qel(j)%id = id
         END IF
      END DO
      j = 0
      DO i=1, NS%NTion
         READ(15,*) dummy, x, y, z, id
         IF (i > Ntemp2 .AND. j < NS%Nion) THEN
            j = j + 1
            qion(j)%xx(1) = x
            qion(j)%xx(2) = y
            qion(j)%xx(3) = z
            qion(j)%id = id
         END IF
      END DO
      CLOSE(15)
   ELSE
      !
      ! Restart - read previous velocity
      OPEN (UNIT=15, file="restart.xyz", STATUS = "OLD", IOSTAT= OpenStatus)
      IF (OpenStatus > 0) THEN
         CALL MPI_Finalize(COMM%rc)
         STOP " *** CANNOT OPEN restart.xyz *** "
      END IF

      READ(15,*) dummy
      READ(15,*) dummy
      xm = param%xm(1)
      j = 0
      DO i=1, NS%NTel
         READ(15,*) dummy, x, y, z, vx, vy, vz, fx, fy, fz, id
         IF (i > Ntemp1 .AND. j < NS%Nel) THEN
            j = j + 1
            qel(j)%xx(1) = x
            qel(j)%xx(2) = y
            qel(j)%xx(3) = z

            qel(j)%xv(1) = vx
            qel(j)%xv(2) = vy
            qel(j)%xv(3) = vz

            qel(j)%ff(1) = fx
            qel(j)%ff(2) = fy
            qel(j)%ff(3) = fz
            
            qel(j)%id = id

            sys%mv2el = sys%mv2el + &
                 xm*(qel(j)%xv(1)**2 + qel(j)%xv(2)**2 + qel(j)%xv(3)**2)
         END IF
      END DO
      j = 0
      xm = param%xm(2)
      DO i=1, NS%NTion
         READ(15,*) dummy, x, y, z, vx, vy, vz, fx, fy, fz, id
         IF (i > Ntemp2 .AND. j < NS%Nion) THEN
            j = j + 1
            qion(j)%xx(1) = x
            qion(j)%xx(2) = y
            qion(j)%xx(3) = z

            qion(j)%xv(1) = vx
            qion(j)%xv(2) = vy
            qion(j)%xv(3) = vz

            qion(j)%ff(1) = fx
            qion(j)%ff(2) = fy
            qion(j)%ff(3) = fz
            
            qion(j)%id = id

            sys%mv2ion = sys%mv2ion + &
                 xm*(qion(j)%xv(1)**2 + qion(j)%xv(2)**2 + qion(j)%xv(3)**2)
         END IF
      END DO
      CLOSE(15)

    END IF
ELSE
   !
   ! Dummy data
   qel(1)%xx(:) = 0.0D0
   qion(1)%xx(:) = 0.0D0
END IF
!
! temperature estimation for restart
IF (param%rest /= 'OFF') THEN
   CALL MPI_REDUCE(sys%mv2el, x, 1, MPI_REAL8, MPI_SUM, 0, &
        COMM%COMM1D, IERR)
   CALL MPI_REDUCE(sys%mv2ion, y, 1, MPI_REAL8, MPI_SUM, 0, &
        COMM%COMM1D, IERR)
   sys%temp = x + y
   CALL MPI_BCAST(sys%temp, 1, MPI_REAL8, 0, COMM%COMM1D, IERR)
END IF
!
! Number of particle check
CALL MPI_REDUCE(NS%Nel, Ntemp1,1,MPI_INTEGER, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%Nion, Ntemp2,1,MPI_INTEGER,MPI_SUM, 0, COMM%COMM1D, IERR)

IF (COMM%GID /= 0 .AND. COMM%COMM1D == 0) THEN
   IF (Ntemp1 /= NS%NTel .OR. Ntemp2 /= NS%NTion) THEN
      PRINT *, "parsing error", Ntemp1, Ntemp2, "while given ones are", &
           NS%NTel, NS%NTion
      STOP
   END IF
END IF
!
!STOP
!
!
END SUBROUTINE INIT
!
END PROGRAM TCP_SPME_VAF
