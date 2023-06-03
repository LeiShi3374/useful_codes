!-----------------------------------------------------------------------
!
!     This module defines data structures for cardiac electrophysiology
!     model equation. It also interfaces with individual modules for
!     the cellular activation model.
!
!-----------------------------------------------------------------------

      MODULE CEPMOD
      USE TYPEMOD
      USE ECMOD
      USE APMOD
      USE FNMOD
      USE TTPMOD
      USE BOMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

!     Type of cardiac electrophysiology models: Aliev-Panfilov model,
!     Bueno-Orovio-Cherry-Fenton model, Fitzhugh-Nagumo model,
!     tenTusscher-Panfilov 2006 model
      INTEGER(KIND=IKIND), PARAMETER :: cepModel_NA = 100,
     2   cepModel_DCPLD = 101, cepModel_AP = 102, cepModel_BO = 103,
     3   cepModel_FN = 104, cepModel_TTP = 105

!     Time integration scheme: Forward-Euler, Runge-Kutta 4th order,
!     Crank-Nicholson
      INTEGER(KIND=IKIND), PARAMETER :: tIntType_NA  = 200,
     2   tIntType_FE = 201, tIntType_RK4 = 202, tIntType_CN2 = 203,
     3   tIntType_BE = 204

!!     Type of excitation-contraction coupling for electromechanics
!      INTEGER(KIND=IKIND), PARAMETER :: ecType_NA = 300,
!     2   ecType_ss_cpld = 301, ecType_ss_dcpld = 302,
!     3   ecType_sn_tiso = 303, ecType_sn_ortho = 304,
!     4   ecType_sn_hetero = 305

!     Time integration scheme and related parameters
      TYPE odeType
!        Time integration method type
         INTEGER(KIND=IKIND) :: tIntType = tIntType_NA
!        Max. iterations for Newton-Raphson method
         INTEGER(KIND=IKIND) :: maxItr = 5
!        Absolute tolerance
         REAL(KIND=RKIND) :: absTol = 1.E-8_RKIND
!        Relative tolerance
         REAL(KIND=RKIND) :: relTol = 1.E-4_RKIND
      END TYPE odeType

!     External stimulus type
      TYPE stimType
!        start time
         REAL(KIND=RKIND) :: Ts = 0._RKIND
!        duration of stimulus
         REAL(KIND=RKIND) :: Td = 0._RKIND
!        end time
         REAL(KIND=RKIND) :: Te = 0._RKIND
!        cycle length
         REAL(KIND=RKIND) :: CL = 0._RKIND
!        stimulus amplitude
         REAL(KIND=RKIND) :: A = 0._RKIND
      END TYPE stimType

!     S1S2 protocol data type
      TYPE S1S2type
!        Num S1 repeats before S2 stimulus
         INTEGER(KIND=IKIND) :: nrep = 0
!        Counter to track S1 repeats
         INTEGER(KIND=IKIND) :: cntr = 1
!        Diastolic interval
         REAL(KIND=RKIND) :: DI = 0._RKIND
!        Action potential duration
         REAL(KIND=RKIND) :: APD = 0._RKIND
!        S2 stimulus
         REAL(KIND=RKIND) :: Istim_A
      END TYPE S1S2type

!     Excitation-contraction coupling type
      TYPE ecCpldType
!        Active stress coupling
         LOGICAL :: astress = .FALSE.
!        Active strain coupling
         LOGICAL :: astrain = .FALSE.
!        Time integration options
         TYPE(odeType) :: odeS
!        Excitation-contraction coupling
!         INTEGER(KIND=IKIND) :: ecType = ecType_NA
!        Active stress/strain state variable
         REAL(KIND=RKIND) :: Ta
      END TYPE ecCpldType

!     Cardiac electrophysiology model type
      TYPE cepModelType
!        Type of cardiac electrophysiology model
         INTEGER(KIND=IKIND) :: cepType = cepModel_NA
!        Myocardium zone id: 1-epi; 2-endo; 3-myo
         INTEGER(KIND=IKIND) :: imyo = 1
!        External stimulus
         TYPE(stimType) :: Istim
!        Time integration options
         TYPE(odeType) :: odeS
!        S1S2 type
         TYPE(S1S2type) :: S1S2
!        Excitation-contraction coupling
         TYPE(ecCpldType) :: ec
      END TYPE cepModelType

!     Initialize parameters from file
      LOGICAL :: init_param_file

!     Input parameters file path
      CHARACTER(LEN=stdL) :: fparam_in

!     Number of state variables
      INTEGER(KIND=IKIND) :: nX = 0

!     Number of gating variables
      INTEGER(KIND=IKIND) :: nG = 0

!     State variables
      REAL(KIND=RKIND), ALLOCATABLE :: X(:)

!     Gating variables
      REAL(KIND=RKIND), ALLOCATABLE :: Xg(:)

      END MODULE CEPMOD
!#######################################################################
