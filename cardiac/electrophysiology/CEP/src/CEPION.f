!-----------------------------------------------------------------------
!
!     This routine embodies formulation for solving electrophysiology
!     cellular activation model at the nodal level.
!
!-----------------------------------------------------------------------

!     Initialize Cardiac Electrophysiology Model
      SUBROUTINE CEPINIT()
      USE COMMOD
      IMPLICIT NONE

      SELECT CASE (cep%cepType)
      CASE (cepModel_DCPLD)
!        Clear log file
         WRITE(oFile,'(A)') "log_dcpld.txt"
         CALL CTXT(oFile)

!        Initialize active stress/strain variable
         cep%ec%Ta = 0._RKIND

!        Overwrite parameters with user-provided input file
         IF (init_param_file) THEN
            std = " Reading parameters from input file"
            CALL EC_READPARFF(fparam_in)
         END IF

      CASE (cepModel_AP)
!        Clear log file
         WRITE(oFile,'(A)') "log_AP.txt"
         CALL CTXT(oFile)

!        Initialize state variables
         std = " Initializing Aliev-Panfilov model"
         CALL AP_INIT(nX, X)

!        Initialize active stress/strain variable
         cep%ec%Ta = 0._RKIND

!        Overwrite parameters with user-provided input file
         IF (init_param_file) THEN
            std = " Reading parameters from input file"
            CALL AP_READPARFF(fparam_in)
         END IF

      CASE (cepModel_BO)
!        Clear log file
         IF (cep%imyo .EQ. 1) THEN
            WRITE(oFile,'(A)') "log_BO_epi.txt"
         ELSE IF (cep%imyo .EQ. 2) THEN
            WRITE(oFile,'(A)') "log_BO_endo.txt"
         ELSE IF (cep%imyo .EQ. 3) THEN
            WRITE(oFile,'(A)') "log_BO_myo.txt"
         END IF
         CALL CTXT(oFile)

!        Initialize state variables
         std = " Initializing Bueno-Orovio model"
         CALL BO_INIT(nX, X)

!        Initialize active stress/strain variable
         cep%ec%Ta = 0._RKIND

!        Overwrite parameters with user-provided input file
         IF (init_param_file) THEN
            std = " Reading parameters from input file"
            CALL BO_READPARFF(fparam_in)
         END IF

      CASE (cepModel_FN)
!        Clear log file
         WRITE(oFile,'(A)') "log_FN.txt"
         CALL CTXT(oFile)

!        Initialize state variables
         std = " Initializing Fitzhugh-Nagumo model"
         CALL FN_INIT(nX, X)

!        Initialize active stress/strain variable
         cep%ec%Ta = 0._RKIND

!        Overwrite parameters with user-provided input file
         IF (init_param_file) THEN
            std = " Reading parameters from input file"
            CALL FN_READPARFF(fparam_in)
         END IF

      CASE (cepModel_TTP)
!        Clear log file
         IF (cep%imyo .EQ. 1) THEN
            WRITE(oFile,'(A)') "log_TTP_epi.txt"
         ELSE IF (cep%imyo .EQ. 2) THEN
            WRITE(oFile,'(A)') "log_TTP_endo.txt"
         ELSE IF (cep%imyo .EQ. 3) THEN
            WRITE(oFile,'(A)') "log_TTP_myo.txt"
         END IF
         CALL CTXT(oFile)

!        Initialize state variables
         std = " Initializing tenTusscher-Panfilov model"
         CALL TTP_INIT(cep%imyo, nX, nG, X, Xg)

!        Initialize active stress/strain variable
         cep%ec%Ta = 0._RKIND

!        Overwrite parameters with user-provided input file
         IF (init_param_file) THEN
            std = " Reading parameters from input file"
            CALL TTP_READPARFF(fparam_in)
         END IF

      END SELECT

      RETURN
      END SUBROUTINE CEPINIT
!#######################################################################
!     State variable integration
      SUBROUTINE CEPINTEG()
      USE COMMOD
      IMPLICIT NONE

      INTEGER, PARAMETER :: iwrite = 1
      INTEGER i, j, k, no
      REAL(KIND=8) :: cTS, Istim, X0, Ksac
      LOGICAL flag, ecCpld

      INTEGER, ALLOCATABLE :: IPAR(:)
      REAL(KIND=8), ALLOCATABLE :: Xo(:), RPAR(:)

!     Output to files
      no = nX
      IF (cep%cepType .EQ. cepModel_BO) THEN
         no = no + 3
         ALLOCATE(RPAR(5))
      ELSE IF (cep%cepType .EQ. cepModel_TTP) THEN
         no = no + 16
         ALLOCATE(RPAR(18))
      ELSE
         ALLOCATE(RPAR(2))
      END IF
      RPAR(:) = 0._RKIND

      ecCpld = cep%ec%astress .OR. cep%ec%astrain
      IF (cep%ec%astress .OR. cep%ec%astrain) no = no + 1

      ALLOCATE(Xo(no))
      Xo = 0._RKIND

      IF (cep%odeS%tIntType .EQ. tIntType_CN2) THEN
         ALLOCATE(IPAR(2))
         IPAR(1) = cep%odeS%maxItr
         IPAR(2) = 0
         RPAR(1) = cep%odeS%absTol
         RPAR(2) = cep%odeS%relTol
      END IF

!     Feedback coefficient for stretch-activated-currents
      Ksac = 0._RKIND

      IF (iProg) THEN
         WRITE(*,'(A)',ADVANCE='NO') " Time integration progress: "
      END IF

      SELECT CASE (cep%cepType)
      CASE (cepModel_DCPLD)
!        Excitation-contraction coupling due to active stress
         IF (cep%ec%astress) THEN
            SELECT CASE (cep%ec%odeS%tIntType)
            CASE (tIntType_FE)
!              Time loop
               j   = 1
               k   = 1
               cTS = 0._RKIND
               DO i=1, nTS
                  CALL EC_ACTVSTRS_FE(cTS, dt, cep%ec%Ta)
                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(1) = cep%ec%Ta
                  cTS   = cTS + dt

!                 Write state variables to a log file
                  IF (MOD(i,iwrite) .EQ. 0)
     2               CALL WTXT(oFile, no, Xo, cTS)

!                 Display progress
                  IF (i.EQ.j .AND. iProg) THEN
                     WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                     k = k + 1
                     j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
                  END IF
               END DO

            CASE (tIntType_RK4)
!              Time loop
               j   = 1
               k   = 1
               cTS = 0._RKIND
               DO i=1, nTS
                  CALL EC_ACTVSTRS_RK(cTS, dt, cep%ec%Ta)
                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(1) = cep%ec%Ta
                  cTS   = cTS + dt

!                 Write state variables to a log file
                  IF (MOD(i,iwrite) .EQ. 0)
     2               CALL WTXT(oFile, no, Xo, cTS)

!                 Display progress
                  IF (i.EQ.j .AND. iProg) THEN
                     WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                     k = k + 1
                     j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
                  END IF
               END DO

            CASE (tIntType_BE)
!              Time loop
               j   = 1
               k   = 1
               cTS = 0._RKIND
               DO i=1, nTS
                  CALL EC_ACTVSTRS_BE(cTS, dt, cep%ec%Ta)
                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(1) = cep%ec%Ta
                  cTS   = cTS + dt

!                 Write state variables to a log file
                  IF (MOD(i,iwrite) .EQ. 0)
     2               CALL WTXT(oFile, no, Xo, cTS)

!                 Display progress
                  IF (i.EQ.j .AND. iProg) THEN
                     WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                     k = k + 1
                     j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
                  END IF
               END DO

            END SELECT

!        Excitation-contraction coupling due to active strain
         ELSE IF (cep%ec%astrain) THEN
            cTS = 0._RKIND
            DO i=1, nTS
               CALL EC_ACTVSTRN(cTS, cep%ec%Ta)
               Xo(1) = cep%ec%Ta

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!              Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF

               cTS = cTS + dt
            END DO
         END IF

      CASE (cepModel_AP)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(1)
               CALL AP_INTEGFE(nX, X, cTS, dt, Istim, Ksac)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(1)
               CALL AP_INTEGRK(nX, X, cTS, dt, Istim, Ksac)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(1)
               CALL AP_INTEGCN2(nX, X, cTS, dt, Istim, Ksac, IPAR, RPAR)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_BO)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (cep%ec%astress) THEN
                  X0 = X(1)
               ELSE IF (cep%ec%astrain) THEN
                  X0 = X(4)
               END IF

!              Integrate local state variables
               CALL BO_INTEGFE(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            RPAR)

!              Output quantities
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (cep%ec%astress) THEN
                  X0 = X(1)
               ELSE IF (cep%ec%astrain) THEN
                  X0 = X(4)
               END IF

!              Integrate local state variables
               CALL BO_INTEGRK(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            RPAR)

!              Output quantities
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (cep%ec%astress) THEN
                  X0 = X(1)
               ELSE IF (cep%ec%astrain) THEN
                  X0 = X(4)
               END IF

!              Integrate local state variables
               CALL BO_INTEGCN2(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            IPAR, RPAR)

!              Output quantities
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_FN)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGFE(nX, X, cTS, dt, Istim)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGRK(nX, X, cTS, dt, Istim)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGCN2(nX, X, cTS, dt, Istim, IPAR, RPAR)

!              Output quantities
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_TTP)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGFE(cep%imyo, nX, nG, X, Xg, dt, Istim,
     2            Ksac, RPAR)

!              Output quantities
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGRK(cep%imyo, nX, nG, X, Xg, dt, Istim,
     2            Ksac, RPAR)

!              Output quantities
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0._RKIND
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGCN2(cep%imyo, nX, nG, X, Xg, dt, Istim,
     2            Ksac, IPAR, RPAR)

!              Output quantities
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) THEN
                  err = " NaN occurence (X). Aborted!"
               END IF

!              Excitation-contraction coupling due to active stress
               IF (cep%ec%astress) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, dt, cep%ec%Ta)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Excitation-contraction coupling due to active strain
               IF (cep%ec%astrain) THEN
                  IF (cep%ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, dt, cep%ec%Ta)

                  ELSE IF (cep%ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, dt, cep%ec%Ta,
     2                  cep%ec%odeS%maxItr, cep%ec%odeS%absTol,
     3                  cep%ec%odeS%relTol)
                  END IF

                  IF (ISNAN(cep%ec%Ta)) THEN
                     err = " NaN occurence (Ta). Aborted!"
                  END IF
                  Xo(no) = cep%ec%Ta
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      END SELECT

      WRITE(*,'(A)')
      IF (cep%odeS%tIntType .EQ. tIntType_CN2) THEN
         WRITE(*,'(4X,A)') "-------------------------------------------"
         WRITE(*,'(4X,A)') " No of Newton-Raphson fails: "//STR(IPAR(2))
         WRITE(*,'(4X,A)') "-------------------------------------------"
      END IF

      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      IF (ALLOCATED(RPAR)) DEALLOCATE(RPAR)
      IF (ALLOCATED(Xo))   DEALLOCATE(Xo)

      RETURN
      END SUBROUTINE CEPINTEG
!-----------------------------------------------------------------------
      SUBROUTINE GET_STIM_AMP(time, cepL, stim_amp, flag)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: time
      TYPE(cepModelType), INTENT(INOUT) :: cepL
      REAL(KIND=8), INTENT(OUT) :: stim_amp
      LOGICAL, INTENT(OUT) :: flag

      flag = .FALSE.
#ifdef S1S2REST
      IF (cepL%S1S2%cntr .LT. cepL%S1S2%nrep) THEN
#endif
         stim_amp = 0D0
         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%Istim%A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%Istim%CL
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
#ifdef S1S2REST
            cepL%S1S2%cntr = cepL%S1S2%cntr + 1
#endif
         END IF

#ifdef S1S2REST
      ELSE IF (cepL%S1S2%cntr .EQ. cepL%S1S2%nrep) THEN

         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%Istim%A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%S1S2%APD + cepL%S1S2%DI
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
            cepL%S1S2%cntr = cepL%S1S2%cntr + 1
            std = " Ts_S2: "//STR(cepL%Istim%Ts)//" DI: "//
     2         STR(cepL%S1S2%DI)
         END IF

      ELSE IF (cep%S1S2%cntr .EQ. cep%S1S2%nrep+1) THEN

         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%S1S2%Istim_A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%Istim%CL
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
            cepL%S1S2%cntr = 0

!           Adjust DI
            IF (cepL%S1S2%DI .GT. 1000.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 1000.0D0
            ELSE IF (cepL%S1S2%DI .GT. 300.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 100.0D0
            ELSE IF (cepL%S1S2%DI .GT. 150.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 20.0D0
            ELSE IF (cepL%S1S2%DI .GT. 50.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 10.0D0
            ELSE IF (cepL%S1S2%DI .GT. 5.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 5.0D0
            ELSE
               flag = .TRUE.
               WRITE(*,'(A)')
               std = " S1S2 restitution protocol complete!"
               RETURN
            END IF
         END IF
      END IF
#endif

      RETURN
      END SUBROUTINE GET_STIM_AMP
!####################################################################
