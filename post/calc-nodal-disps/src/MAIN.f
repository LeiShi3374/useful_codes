!====================================================================
!
!
!
!====================================================================

      PROGRAM MAIN
      USE COMMOD
      IMPLICIT NONE

      INTEGER :: i, j, a, e, Ac, eNoN, fid
      CHARACTER(LEN=stdL) :: fName

      REAL(KIND=RKIND), ALLOCATABLE :: xi(:), dl(:), N(:), Nxi(:,:)

      CALL READFILES()

      ALLOCATE(Dn(nsd,msh(1)%nNo))
      Dn = 0._RKIND

      DO i=1, nprb
         CALL FINDPRBTRACE(msh(1), xprb(:,i), trace(i))
      END DO

      eNoN = msh(1)%eNoN
      ALLOCATE(xi(nsd), dl(nsd), N(eNoN), Nxi(nsd,eNoN))

      fid = 1125
      DO i=1, nprb
         WRITE(fName,'(A)') "dispf_p"//STR(i)//".dat"
         OPEN(fid, FILE=TRIM(fName))
         CLOSE(fid, STATUS='DELETE')

         OPEN(fid, FILE=TRIM(fName))
         WRITE(fid,'(A)') "Variables = t, u_x, u_y, u_z"
         CLOSE(fid)
      END DO

      DO cTS=startTS, endTS, incrTS
         time = REAL(cTS,KIND=RKIND) * dt

         Dn   = 0._RKIND
         IF (cTS .LT. 100) THEN
            WRITE(fName,'(A,I3.3,A)') TRIM(saveName)//"_", cTS, ".vtu"
         ELSE
            WRITE(fName,'(A)') TRIM(saveName)//"_"//STR(cTS)//".vtu"
         END IF
         CALL READVTUDISP(fName, msh(1)%nNo, Dn)

         DO i=1, nprb
            e  = trace(i)%gE
            xi = trace(i)%xi
            CALL GETGNN(nsd, msh(1)%eType, eNoN, xi, N, Nxi)

            dl = 0._RKIND
            DO a=1, eNoN
               Ac = msh(1)%IEN(a,e)
               dl = dl + N(a)*Dn(:,Ac)
            END DO

            WRITE(fName,'(A)') "dispf_p"//STR(i)//".dat"
            OPEN(fid, FILE=TRIM(fName), POSITION='APPEND')
            WRITE(fid,'(F9.4)',ADVANCE='NO') time
            DO j=1, nsd
               WRITE(fid,'(3X,1pE15.6)',ADVANCE='NO') dl(j)
            END DO
            WRITE(fid,'(A)')
            CLOSE(fid)
         END DO
      END DO

      DEALLOCATE(xi, dl, N, Nxi)

      CALL FINALIZE()

      END PROGRAM MAIN

!====================================================================
      SUBROUTINE STOPSIM()
      IMPLICIT NONE

      CALL FINALIZE()
      STOP

      END SUBROUTINE STOPSIM
!====================================================================
      SUBROUTINE FINALIZE()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) i, iM

!     Deallocate mesh
      IF (ALLOCATED(msh)) THEN
         DO iM=1, nMsh
            CALL DESTROY(msh(iM))
         END DO
         DEALLOCATE(msh)
      END IF

      IF (ALLOCATED(trace)) THEN
         DO i=1, nprb
            IF (ALLOCATED(trace(i)%xi)) DEALLOCATE(trace(i)%xi)
         END DO
         DEALLOCATE(trace)
      END IF

      IF (ALLOCATED(x))    DEALLOCATE(x)
      IF (ALLOCATED(Dn))   DEALLOCATE(Dn)
      IF (ALLOCATED(xprb)) DEALLOCATE(xprb)

!     Closing the output channels
      CALL std%close()
      CALL wrn%close()
      CALL err%close()
      CALL dbg%close()

      RETURN
      END SUBROUTINE FINALIZE
!====================================================================
