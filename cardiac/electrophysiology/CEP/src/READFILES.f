c====================================================================
c
c     This set of subroutines reads inputs from command line and
c     from the input file
c
c====================================================================

      SUBROUTINE READ_INPUTS()
      USE COMMOD
      USE LISTMOD
      IMPLICIT NONE

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: narg
      CHARACTER(LEN=stdL) :: ctmp

      TYPE(listType) :: list
      TYPE(listType), POINTER :: lPtr, lSub

c--------------------------------------------------------------------
c     Set communicator channels and its pointers
      CALL io%o%new(CHNL_O, tag="", OTS=.TRUE.)
      CALL io%e%new(CHNL_E, tag="")
      CALL io%w%new(CHNL_W, tag="")
      CALL io%d%new(CHNL_D, tag="")

      std => io%o
      err => io%e
      wrn => io%w
      dbg => io%d

c--------------------------------------------------------------------
c     Get the command line arguments
      narg = COMMAND_ARGUMENT_COUNT()
      IF (narg .NE. 1) THEN
         err = " ERROR: Required input argument not provided"
      END IF

      CALL GET_COMMAND_ARGUMENT(1, ctmp)
      flag = .FALSE.
      INQUIRE(FILE=TRIM(ctmp), EXIST=flag)
      IF (flag) THEN
         cep_fIn = TRIM(ctmp)
         std = " Reading input file "//TRIM(cep_fIn)
      ELSE
         err = " Input file doesn't exist"
      END IF

c--------------------------------------------------------------------
c     Read input file
      list = listType(cep_fIn, io)

c--------------------------------------------------------------------
c     Load inputs from list

      lSub => list%get(ctmp, "Electrophysiology model")
      CALL TO_LOWER(ctmp)
      SELECT CASE (TRIM(ctmp))
      CASE ("decoupled")
         cep%cepType = cepModel_DCPLD
         nX = 0
         nG = 0

      CASE ("ap", "aliev_panfilov")
         cep%cepType = cepModel_AP
         nX = 2
         nG = 0

      CASE ("bo", "bueno_orovio")
         cep%cepType = cepModel_BO
         nX = 4

      CASE ("fn","fitzhugh_nagumo")
         cep%cepType = cepModel_FN
         nX = 2
         nG = 0

      CASE ("ttp","tentusscher_panfilov")
         cep%cepType = cepModel_TTP
         nX = 7
         nG = 12

      CASE DEFAULT
         err = " Unknown electrophysiology model"
      END SELECT

      fparam_in = ""
      init_param_file = .FALSE.

      lPtr => lSub%get(ctmp, "CEP model parameters file path")
      IF (ASSOCIATED(lPtr)) THEN
         flag = .FALSE.
         INQUIRE(FILE=TRIM(ctmp), EXIST=flag)
         IF (flag) THEN
            fparam_in = TRIM(ctmp)
            init_param_file = .TRUE.
         ELSE
            err = " CEP model parameters file doesn't exist"
         END IF
      END IF

      cep%imyo = 1
      lPtr => list%get(ctmp, "Myocardial zone")
      IF (ASSOCIATED(lPtr)) THEN
         CALL TO_LOWER(ctmp)
         SELECT CASE (TRIM(ctmp))
         CASE ("epi", "epicardium")
            cep%imyo = 1

         CASE ("endo", "endocardium", "pfib", "purkinje")
            cep%imyo = 2

         CASE ("myo", "mid-myo", "myocardium")
            cep%imyo = 3

         CASE DEFAULT
            err = "Undefined myocardium zone"
         END SELECT
      END IF

      cep%odes%tIntType = tIntType_RK4
      lPtr => list%get(ctmp, "Time integration scheme")
      IF (ASSOCIATED(lPtr)) THEN
         CALL TO_LOWER(ctmp)
         SELECT CASE (ctmp)
         CASE ("fe", "euler", "forward_euler", "explicit")
            cep%odes%tIntType = tIntType_FE

         CASE ("rk", "rk4", "runge")
            cep%odes%tIntType = tIntType_RK4

         CASE ("cn", "cn2", "implicit")
            cep%odes%tIntType = tIntType_CN2
            IF (cep%cepType .EQ. cepModel_TTP) THEN
               err = "Implicit time integration for tenTusscher-"//
     2            "Panfilov model can give unexpected results. "//
     3            "Use FE or RK4 instead"
            END IF

         CASE DEFAULT
            err = " Unknown ODE time integration scheme"
         END SELECT
      END IF

      lPtr => list%get(nTS, "Number of time steps")
      lPtr => list%get(dt, "Time step size")

      lSub => list%get(ctmp, "Stimulus")
      IF (ASSOCIATED(lSub)) THEN
         lPtr => lSub%get(cep%Istim%A, "Amplitude")
         IF (.NOT.ISZERO(cep%Istim%A)) THEN
            lPtr => lSub%get(cep%Istim%Ts, "Start time")
            lPtr => lSub%get(cep%Istim%Td, "Duration")
            lPtr => lSub%get(cep%Istim%CL, "Cycle length")
            IF (.NOT.ASSOCIATED(lPtr)) THEN
               cep%Istim%CL = REAL(nTS, KIND=RKIND) * dt
            END IF
            cep%Istim%Te = cep%Istim%Ts + cep%Istim%Td
         END IF
      END IF

      lSub => list%get(ctmp, "Excitation-contraction coupling")
      IF (ASSOCIATED(lSub)) THEN
         CALL TO_LOWER(ctmp)
         SELECT CASE (TRIM(ctmp))
         CASE ("active_stress", "stress")
            cep%ec%astress = .TRUE.
            std = " Excitation-contraction coupling: active stress"

         CASE ("active_strain", "strain")
            cep%ec%astrain = .TRUE.
            std = " Excitation-contraction coupling: active strain"

         CASE DEFAULT
            err = " Unknown excitation-contraction coupling"
         END SELECT

         flag = (cep%cepType.EQ.cepModel_DCPLD) .AND.
     2          (cep%ec%astrain)
         IF (.NOT.flag) THEN
            lPtr => lSub%get(ctmp, "Time integration scheme")
            SELECT CASE (TRIM(ctmp))
            CASE ("fe", "euler", "forward_euler", "explicit")
               cep%ec%odeS%tIntType = tIntType_FE

            CASE ("rk", "rk4", "runge")
               cep%ec%odeS%tIntType = tIntType_RK4

            CASE ("be", "backward_euler", "implicit")
               cep%ec%odeS%tIntType = tIntType_BE

            CASE ("cn", "cn2")
               err = " Crank-Nicholson time integration cannot "//
     2            "be used for excitation-contraction coupling"

            CASE DEFAULT
               err = " Unknown ODE time integration scheme"
            END SELECT
         END IF

         IF (cep%cepType .EQ. cepModel_DCPLD) THEN
            lPtr => lSub%get(ctmp, "Parameters file path")
            IF (ASSOCIATED(lPtr)) THEN
               flag = .FALSE.
               INQUIRE(FILE=TRIM(ctmp), EXIST=flag)
               IF (flag) THEN
                  fparam_in = TRIM(ctmp)
                  init_param_file = .TRUE.
               ELSE
                  err = " Decoupled EC parameters file doesn't exist"
               END IF
            END IF
         END IF
      END IF

      lSub => list%get(ctmp, "S1-S2 protocol")
      IF (ASSOCIATED(lSub)) THEN
         IF (cep%cepType .NE. cepModel_TTP) THEN
            err = " S1S2 protocol applicable for TTP model only"
         END IF
         lPtr => lSub%get(cep%S1S2%nrep,
     2      "Number of S1 repetitions before S2")
         lPtr => lSub%get(cep%S1S2%DI, "Initial diastolic interval")
         lPtr => lSub%get(cep%S1S2%APD, "Basic APD")
         lPtr => lSub%get(cep%S1S2%Istim_A, "S2 amplitude")
         iProg = .FALSE.
      END  IF

c--------------------------------------------------------------------
c     Check inputs for any inconsistencies
      flag = cep%ec%astrain .OR. cep%ec%astress
      IF (flag .AND. cep%cepType.EQ.cepModel_FN) THEN
         err = " Excitation-contraction coupling is not allowed for "//
     2      " Fitzhugh-Nagumo model"
      END IF

      IF (cep%cepType.EQ.cepModel_AP .AND. cep%ec%astrain) THEN
         err = " Active-strain coupling is allowed for TTP and BO"//
     2      " activation models only"
      END IF

      IF (cep%cepType.EQ.cepModel_BO .AND. cep%ec%astrain) THEN
         wrn = " Active-strain coupling with Bueno-Orovio model"//
     2      " can lead to unexpected results"
      END IF

      IF (cep%imyo.LT.1 .OR. cep%imyo.GT.3) THEN
         err = " Invalid myocardial zone specified"
      END IF

      flag = cep%cepType.NE.cepModel_TTP .AND.
     2       cep%cepType.NE.cepModel_BO
      IF (flag .AND. cep%imyo.GT.1) THEN
         err = " Mid-myocardium and endocardium zones are allowed"//
     2      " only for TTP and BO activation models."
      END IF

c--------------------------------------------------------------------
c     Allocate arrays
      ALLOCATE(X(nX), Xg(nG))
      X  = 0._RKIND
      Xg = 0._RKIND

      RETURN
      END SUBROUTINE READ_INPUTS
c====================================================================
