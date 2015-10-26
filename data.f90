MODULE DATASTR
IMPLICIT NONE
!
! PRECISION
! SI = single precision for integer
! DI = double precision for integer
! SP = single precision for real
! DP = double precision for real
!INTEGER, PARAMETER:: SI = SELECTED_INT_KIND(7),  DI = SELECTED_INT_KIND(15)
!INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(7), DP = SELECTED_REAL_KIND(14)
INTEGER, PARAMETER:: SI = 4, DI = 8
INTEGER, PARAMETER:: SP = 4, DP = 8
! Nparam = Number of particle kinds
! Ndist = RDF index
! Nref = 3 for PBC, 6 for vacuum
INTEGER, PARAMETER:: Nparam = 2, Ndist = 500, Nref = 0
INTEGER, PARAMETER:: Nfine = 10   ! Refining size
INTEGER, PARAMETER:: Ngmax = 8   ! Maximum size of local group communicator
INTEGER, PARAMETER:: Nrefn_limit = 4 ! Maximum number of refining
INTEGER, PARAMETER:: Ngrid = 5000 ! Grid for force/energy interpolation
REAL(KIND=DP), PARAMETER:: eps = 14.399644154020082D0 
REAL(KIND=DP), PARAMETER:: tps = 10.180505306645389D0
REAL(KIND=DP), PARAMETER:: pi  = 3.1415926535897931D0
REAL(KIND=DP), PARAMETER:: sqrtpi = 1.7724538509055159D0
!
! Particle data type
! xx = current position
! xv = current velocity
! ff = current force
! xo = old position
! vo = old velocity
TYPE PT
   INTEGER:: id
   REAL(KIND=DP):: xx(3), xv(3), ff(3), vi(3)
END TYPE PT
!
! Ghost particle data type
! xx = position for initial communication/force for later communication
TYPE GH
   REAL(KIND=DP) :: xx(3)
END TYPE GH
!
! Time data type
! Nloop = Number of loops
! Ndump = Number of dumps
! Nsamp = Number of samplings
! tmax = Maximum simulation time
! tdump = xyz file dump time
! tsamp = sampling time
! tnow = Current time
TYPE TM
   INTEGER:: Nloop, Ndump, Nsamp, Nrest
   REAL(KIND=DP) :: tmax, tdump, tsamp, trest, tnow
END TYPE TM
!
! Parameter data type
! xm = Mass of particle
! tau = Berendsen thermostat parameter
! alpha = Statistical thermostat damping parameter
! rs = thermal de Broglie wave length
! a
! box = Simulation box size
! Te = Initial electron temperature
! Ti = Initial ion temperature
! thermo = Thermostat kind (Berend/Stoch/Notemp) 
! Nmax = Maximum iteration for Fourier sum
! N2max = Maximum square of iteration for Fourier sum
TYPE PM
   REAL(KIND=DP):: xm(Nparam), q(Nparam), tau, alpha, a, box(3), Te, Ti
   REAL(KIND=DP):: rs, rc, rc2
   INTEGER:: K(3), Nfft
   CHARACTER*6:: thermo
   CHARACTER*3:: rest, trace
END TYPE PM
!
! System variable data type
! box = Size of simulation cell
! temp = Temperature of (NVT) simulation
! mv2 = Sum of mv^2(=twice of kinetic energy)
! Epot = Potential energy of system
! Te = In-situ electron temperature
! Ti = In-situ ion temperature
! RNel = Temperature translation parameter for electron
! RNion = Temperature translation parameter for ion
! Del = Diffusion coefficient for electron
! Dion = Diffusion coefficient for ion
! Pel = Electron pressure
! Pion = Ion pressure
TYPE ST
   REAL(KIND=DP):: temp, mv2el, mv2ion, Epot, Te, Ti, Uo, Ur, Um, cel, cion, &
        & Pel, Pion, J0(3), Jei(3)
END TYPE ST
!
! Various number type
! Main objective is transfer specific numbers into print routine
! Other sampling variables
TYPE NM
   INTEGER:: Npt, Nel, Nion, NTpt, NLpt, NTel, NTion, Ng, NCPUS, Npe, Nrefn
   REAL(KIND=DP) :: dx_ee, dx_ei, dx_ii, dv_e, dv_i, &
        RDF_ee(Ndist+1), RDF_ei(Ndist+1), RDF_ii(Ndist+1), &
        & VDF_e(Ndist+1), VDF_i(Ndist+1), dx
   INTEGER, POINTER::NOPT(:,:), CNT(:), DPL(:), NGHB(:,:)
END TYPE NM
!
! MPI variables
TYPE MP
   !
   ! Basic parameter
   INTEGER:: NUMPROC, rc, PTCL, LNCPUS, Nlsize, Ncolmn, iy
   !
   ! Communicator and ID
   INTEGER:: GID, GLOBAL, LID, LOCAL, MYID2, COMM1D
   !
   ! Destiny and source for pipeline
   INTEGER, POINTER:: DEST(:), SRC(:)
   LOGICAL:: TAG, ITAG
END TYPE MP
!
CONTAINS
!
! #############################################################################
! MPI coordinate/data type/communicating pair configuration ###################
SUBROUTINE MPI_initialize(COMM, NS)
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(MP):: COMM
CHARACTER(len=256):: hostname
LOGICAL:: PERIOD, REORDER
INTEGER:: IERR, MYID, offset, Oldtype, blockcnt, COLORR, KEYY, Npe
INTEGER:: NUMPROC, Ncolmn, Nres, ix, iy, i, j
!
! Number of processor for PME routine
Npe = 1 
!
! Basic communicator ##########################################################
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
CALL MPI_GET_PROCESSOR_NAME(hostname, i,COMM%rc)
NS%NCPUS = NUMPROC
!
! ############## Global communicator ##########################################
PERIOD  = .TRUE.
REORDER = .TRUE.
CALL MPI_CART_CREATE(MPI_COMM_WORLD, 1, NUMPROC, PERIOD, REORDER, &
     & COMM%GLOBAL, IERR)
CALL MPI_COMM_RANK(COMM%GLOBAL, COMM%GID, IERR)
!
! ############## Local node and ID decision ###################################
Npe = INT((NUMPROC-1)/Ngmax) + 1
Ncolmn = INT(NUMPROC/Npe)
IF (Ncolmn*Npe /= NUMPROC) Ncolmn = Ncolmn + 1
Nres = MOD(NUMPROC,Ncolmn)
!
iy = INT(COMM%GID/Ncolmn)
ix = MOD(COMM%GID,Ncolmn)
COMM%Ncolmn = Ncolmn
IF (MYID == 0) THEN
   IF (Nres /=0) THEN
      PRINT *, "Distribution of PEs"
      PRINT *, " PME node/pipeline nodes --> "
      DO i=1, Npe-1
         PRINT '(5X, I3, 5X, 31(I3,1X))', (j+(i-1)*Ncolmn-1, j=1,Ncolmn) 
      END DO
      PRINT  '(5X, I3, 5X, 31(I3,1X))', (j+(Npe-1)*Ncolmn-1, j=1,Nres) 
      PRINT *, " PME node/pipeline nodes --> "
   ELSE
      PRINT *, "Distribution of PEs"
      PRINT *, " PME node/pipeline nodes --> "
      DO i=1, Npe
         PRINT '(5X, I3, 5X, 31(I3,1X))', (j+(i-1)*Ncolmn-1, j=1,Ncolmn) 
      END DO
      PRINT *, " PME node/pipeline nodes --> "
   END IF
END IF
!
! ############# COMM1D communicator ###########################################
IF (ix == 0) THEN
   COLORR = 0
   COMM%TAG = .FALSE.
ELSE
   COLORR = 1
   COMM%TAG = .TRUE.
END IF
KEYY = COMM%GID
CALL MPI_COMM_SPLIT(COMM%GLOBAL, COLORR, KEYY, COMM%COMM1D, IERR)
CALL MPI_COMM_RANK(COMM%COMM1D, COMM%MYID2, IERR)
CALL MPI_COMM_SIZE(COMM%COMM1D, COMM%NUMPROC, IERR)
!
! ############## Local communicator ###########################################
COLORR = iy
CALL MPI_COMM_SPLIT(COMM%GLOBAL, COLORR, KEYY, COMM%LOCAL, IERR)
CALL MPI_COMM_RANK(COMM%LOCAL, COMM%LID, IERR)
CALL MPI_COMM_SIZE(COMM%LOCAL, COMM%LNCPUS, IERR)
!
IF (MOD(NUMPROC-Npe,2) == 0 ) THEN
   ! 
   ! Even number of processors for pipeline
   COMM%ITAG = .TRUE.
ELSE
   !
   ! Odd number of processors for pipeline
   COMM%ITAG = .FALSE.
END IF
!
! Destination/Source nodes allocation
! Number of particle sets allocation
IF (Nres > 0 .AND. iy == (Npe-1)) THEN
   Ncolmn = Nres
   IF (Nres < 2) STOP "=========== PE distribution error =============="
END IF
ALLOCATE(COMM%DEST(NUMPROC-1-Npe), COMM%SRC(NUMPROC-1-Npe), &
     & NS%NOPT(NUMPROC-Npe,2), NS%CNT(Ncolmn), NS%DPL(Ncolmn), NS%NGHB(Ncolmn,2))
COMM%Nlsize = Ncolmn
COMM%iy = iy
!
! New data type for moving particles
offset = 0
Oldtype = MPI_REAL8
blockcnt = 3
CALL MPI_TYPE_STRUCT(1, blockcnt, offset, Oldtype, COMM%PTCL, IERR)
CALL MPI_TYPE_COMMIT(COMM%PTCL, IERR)
!
DO i=1, NUMPROC-Npe-1
   COMM%DEST(i) = MOD(COMM%MYID2 + i, COMM%NUMPROC)
   COMM%SRC(i)  = MOD(COMM%MYID2 - i + COMM%NUMPROC, COMM%NUMPROC)
END DO
NS%Npe = Npe
!
RETURN
END SUBROUTINE MPI_INITIALIZE
!
END MODULE DATASTR
