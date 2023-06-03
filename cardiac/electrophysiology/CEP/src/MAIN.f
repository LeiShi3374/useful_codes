!--------------------------------------------------------------------
!
!
!--------------------------------------------------------------------

      PROGRAM CEPMAIN
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) it1, it2, rate
      REAL(KIND=8) wtime

      WRITE(*,'(A)') REPEAT('=', 64)

      CALL SYSTEM_CLOCK(count_rate=rate)
      CALL SYSTEM_CLOCK(it1)

      CALL READ_INPUTS()

      CALL CEPINIT()

      CALL CEPINTEG()

      CALL FINALIZE()

      CALL SYSTEM_CLOCK(it2)
      wtime = REAL(it2-it1, KIND=RKIND) / REAL(rate, KIND=RKIND)

      std = " Total elapsed time: "//STR(wtime)//" s"
      WRITE(*,'(A)') REPEAT('=', 64)

      RETURN
      END PROGRAM CEPMAIN
!#######################################################################
      SUBROUTINE FINALIZE()
      USE COMMOD
      IMPLICIT NONE

      IF (ALLOCATED(X)) DEALLOCATE(X)
      IF (ALLOCATED(Xg)) DEALLOCATE(Xg)

      CALL std%close()
      CALL wrn%close()
      CALL err%close()
      CALL dbg%close()

      RETURN
      END SUBROUTINE FINALIZE
!#######################################################################
