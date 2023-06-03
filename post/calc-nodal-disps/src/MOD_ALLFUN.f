!--------------------------------------------------------------------
!
!     All the routines that require an interface are included here.
!     This mainly involves small routines with a very well-defined
!     functionality. To use these routines, just add "USE ALLFUN" to
!     your routine.
!
!--------------------------------------------------------------------

      MODULE ALLFUN
      USE MATFUN
      IMPLICIT NONE

      INTERFACE DESTROY
         MODULE PROCEDURE DESTROYFACE, DESTROYMSH, DESTROYADJ,
     2      DESTROYSTACK, DESTROYQUEUE
      END INTERFACE DESTROY

      INTERFACE GETNADJCNCY
         MODULE PROCEDURE GETNADJ_MSH, GETNADJ_FACE
      END INTERFACE

      INTERFACE GETEADJCNCY
         MODULE PROCEDURE GETEADJ_MSH, GETEADJ_FACE
      END INTERFACE

      CONTAINS
!====================================================================
!     This routine (get to block) searches for "kwd" in file "fid"
!     and returns the position at the next line. If "m" and/or "n"
!     are present, the size of matrix is checked to be compatible
      SUBROUTINE GTBLK(fid, kwd, n)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fid
      CHARACTER(LEN=*), INTENT(IN) :: kwd
      INTEGER(KIND=IKIND), INTENT(OUT) :: n

      INTEGER(KIND=IKIND) l
      CHARACTER(LEN=stdL) rLine

      l = LEN(kwd)
      DO
         READ(fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(2:l+1) .EQ. kwd) EXIT
      END DO
      rLine = rLine(l+2:stdL)
      CALL GET(rLine,n)
      RETURN

 001  err = "Block with keyword <"//kwd//"> was not found"

      END SUBROUTINE GTBLK
!====================================================================
!     These set of routines destroy an object.
      PURE SUBROUTINE DESTROYFACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(OUT) :: lFa

      IF (ALLOCATED(lFa%gE))     DEALLOCATE(lFa%gE)
      IF (ALLOCATED(lFa%gN))     DEALLOCATE(lFa%gN)
      IF (ALLOCATED(lFa%lN))     DEALLOCATE(lFa%lN)
      IF (ALLOCATED(lFa%IEN))    DEALLOCATE(lFa%IEN)
      IF (ALLOCATED(lFa%w))      DEALLOCATE(lFa%w)
      IF (ALLOCATED(lFa%x))      DEALLOCATE(lFa%x)
      IF (ALLOCATED(lFa%xi))     DEALLOCATE(lFa%xi)
      IF (ALLOCATED(lFa%N))      DEALLOCATE(lFa%N)
      IF (ALLOCATED(lFa%nV))     DEALLOCATE(lFa%nV)
      IF (ALLOCATED(lFa%Nx))     DEALLOCATE(lFa%Nx)
      IF (ALLOCATED(lFa%Nxx))    DEALLOCATE(lFa%Nxx)

      CALL DESTROYADJ(lFa%nAdj)
      CALL DESTROYADJ(lFa%eAdj)

      lFa%eType = eType_NA
      lFa%nEl   = 0
      lFa%nNo   = 0

      RETURN
      END SUBROUTINE DESTROYFACE
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYMSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(OUT) :: lM

      INTEGER(KIND=IKIND) iFa

      IF (ALLOCATED(lM%IEN))     DEALLOCATE(lM%IEN)
      IF (ALLOCATED(lM%w))       DEALLOCATE(lM%w)
      IF (ALLOCATED(lM%xib))     DEALLOCATE(lM%xib)
      IF (ALLOCATED(lM%xi))      DEALLOCATE(lM%xi)
      IF (ALLOCATED(lM%N))       DEALLOCATE(lM%N)
      IF (ALLOCATED(lM%Nb))      DEALLOCATE(lM%Nb)
      IF (ALLOCATED(lM%Nx))      DEALLOCATE(lM%Nx)
      IF (ALLOCATED(lM%Nxx))     DEALLOCATE(lM%Nxx)

      IF (ALLOCATED(lM%fa)) THEN
         DO iFa=1, lM%nFa
            CALL DESTROYFACE(lM%fa(iFa))
         END DO
         DEALLOCATE(lM%fa)
      END IF

      CALL DESTROYADJ(lM%nAdj)
      CALL DESTROYADJ(lM%eAdj)

      lM%eType = eType_NA
      lM%nEl   = 0
      lM%nFa   = 0
      lM%nNo   = 0

      RETURN
      END SUBROUTINE DESTROYMSH
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYADJ(adj)
      USE COMMOD
      IMPLICIT NONE
      TYPE(adjType), INTENT(OUT) :: adj

      IF (ALLOCATED(adj%prow)) DEALLOCATE(adj%prow)
      IF (ALLOCATED(adj%pcol)) DEALLOCATE(adj%pcol)
      adj%nnz = 0

      RETURN
      END SUBROUTINE DESTROYADJ
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYSTACK(stk)
      USE COMMOD
      IMPLICIT NONE
      TYPE(stackType), INTENT(OUT) :: stk

      IF (ALLOCATED(stk%v)) DEALLOCATE(stk%v)

      stk%n    = 0
      stk%maxN = 0

      RETURN
      END SUBROUTINE DESTROYSTACK
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYQUEUE(que)
      USE COMMOD
      IMPLICIT NONE
      TYPE(queueType), INTENT(OUT) :: que

      IF (ALLOCATED(que%v)) DEALLOCATE(que%v)

      que%n    = 0
      que%maxN = 0

      RETURN
      END SUBROUTINE DESTROYQUEUE
!====================================================================
!     Finding the mesh ID based on the mesh name
      SUBROUTINE FINDMSH(mshName, iM)
      USE COMMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: mshName
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM

      DO iM=1, nMsh
         IF (msh(iM)%name .EQ. mshName) EXIT
      END DO
      IF (iM .GT. nMsh) err="Unable to find msh <"//TRIM(mshName)//">"

      RETURN
      END SUBROUTINE FINDMSH
!--------------------------------------------------------------------
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE FINDFACE(faceName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: faceName
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM, iFa

      iFa = 0
      MY_LOOP : DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            IF (msh(iM)%fa(iFa)%name .EQ. faceName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP
      IF (iM .GT. nMsh) err="Unable to find face <"//TRIM(faceName)//">"

      RETURN
      END SUBROUTINE FINDFACE
!====================================================================
!     Find nodal adjacency of a given mesh. Computes list of all nodes
!     around a given node of a mesh.
      SUBROUTINE GETNADJ_MSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, i, j, k, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjL(:,:)

      ALLOCATE(incNd(lM%nNo))
      incNd = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            DO b=1, lM%eNoN
               IF (b .EQ. a) CYCLE
               incNd(Ac) = incNd(Ac) + 1
            END DO
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(adjL(maxAdj, lM%nNo))
      incNd = 0
      adjL  = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            DO b=1, lM%eNoN
               IF (b .EQ. a) CYCLE
               Bc = lM%IEN(b,e)
               flag = .TRUE.
               j = 1
               DO i=1, incNd(Ac)
                  IF (Bc .EQ. adjL(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  ELSE IF (Bc .GT. adjL(i,Ac)) THEN
                     j = i + 1
                  END IF
               END DO
               IF (flag) THEN
                  IF (incNd(Ac) .EQ. 0) THEN
                     incNd(Ac)  = 1
                     adjL(1,Ac) = Bc
                  ELSE
                     DO k=incNd(Ac), j, -1
                        adjL(k+1,Ac) = adjL(k,Ac)
                     END DO
                     adjL(j,Ac) = Bc
                     incNd(Ac) = incNd(Ac) + 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE(incNd)

      lM%nAdj%nnz = 0
      DO a=1, lM%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            lM%nAdj%nnz = lM%nAdj%nnz + 1
         END DO
      END DO

      ALLOCATE(lM%nAdj%prow(lM%nNo+1), lM%nAdj%pcol(lM%nAdj%nnz))
      j = 0
      lM%nAdj%prow(1) = j + 1
      DO a=1, lM%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            j = j + 1
            lM%nAdj%pcol(j) = adjL(i,a)
         END DO
         lM%nAdj%prow(a+1) = j + 1
      END DO
      DEALLOCATE(adjL)

      RETURN
      END SUBROUTINE GETNADJ_MSH
!--------------------------------------------------------------------
!     Find nodal adjacency of a given face. Computes list of all nodes
!     around a given node of a face.
      SUBROUTINE GETNADJ_FACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType),  INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, i, j, k, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjL(:,:)

      ALLOCATE(incNd(lFa%nNo))
      incNd = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO b=1, lFa%eNoN
               IF (b .EQ. a) CYCLE
               incNd(Ac) = incNd(Ac) + 1
            END DO
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(adjL(maxAdj, lFa%nNo))
      incNd = 0
      adjL  = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO b=1, lFa%eNoN
               IF (b .EQ. a) CYCLE
               Bc = lFa%IEN(b,e)
               Bc = lfa%lN(Bc)
               flag = .TRUE.
               j = 1
               DO i=1, incNd(Ac)
                  IF (Bc .EQ. adjL(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  ELSE IF (Bc .GT. adjL(i,Ac)) THEN
                     j = i + 1
                  END IF
               END DO
               IF (flag) THEN
                  IF (incNd(Ac) .EQ. 0) THEN
                     incNd(Ac)  = 1
                     adjL(1,Ac) = Bc
                  ELSE
                     DO k=incNd(Ac), j, -1
                        adjL(k+1,Ac) = adjL(k,Ac)
                     END DO
                     adjL(j,Ac) = Bc
                     incNd(Ac)  = incNd(Ac) + 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE(incNd)

      lFa%nAdj%nnz = 0
      DO a=1, lFa%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            lFa%nAdj%nnz = lFa%nAdj%nnz + 1
         END DO
      END DO

      ALLOCATE(lFa%nAdj%prow(lFa%nNo+1), lFa%nAdj%pcol(lFa%nAdj%nnz))
      j = 0
      lFa%nAdj%prow(1) = j + 1
      DO a=1, lFa%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            j = j + 1
            lFa%nAdj%pcol(j) = adjL(i,a)
         END DO
         lFa%nAdj%prow(a+1) = j + 1
      END DO
      DEALLOCATE(adjL)

      RETURN
      END SUBROUTINE GETNADJ_FACE
!====================================================================
!     Find element adjacency of a given mesh. Computes list of all
!     elements around a given element of a mesh.
      SUBROUTINE GETEADJ_MSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjList(:,:),
     2   tmpI(:,:)

      ALLOCATE(incNd(lM%nNo))
      incNd = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            incNd(Ac) = incNd(Ac) + 1
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(tmpI(maxAdj, lM%nNo))
      incNd = 0
      tmpI  = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            incNd(Ac) = incNd(Ac) + 1
            tmpI(incNd(Ac), Ac) = e
         END DO
      END DO
      b = 2*maxAdj

 001  b = b + maxAdj
      DEALLOCATE(incNd)
      ALLOCATE(incNd(lM%nEl))
      IF (ALLOCATED(adjList)) DEALLOCATE(adjList)
      ALLOCATE(adjList(b, lM%nEl))
      adjList = 0
      incNd   = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            DO i=1, maxAdj
               IF (tmpI(i,Ac) .EQ. 0) EXIT
               flag = .TRUE.
               DO j=1, incNd(e)
                  IF (adjList(j,e) .EQ. tmpI(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(e) = incNd(e) + 1
                  IF (incNd(e) .GE. b) GOTO 001
                  adjList(incNd(e),e) = tmpI(i,Ac)
               END IF
            END DO
         END DO
      END DO
      maxAdj = MAXVAL(incNd)
      DEALLOCATE(tmpI, incNd)

      lM%eAdj%nnz = 0
      DO e=1, lM%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               lM%eAdj%nnz = lM%eAdj%nnz + 1
            ELSE
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE(lM%eAdj%prow(lM%nEl+1), lM%eAdj%pcol(lM%eAdj%nnz))
      j = 0
      lM%eAdj%prow(1) = j + 1
      DO e=1, lM%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               j = j + 1
               lM%eAdj%pcol(j) = adjList(i,e)
            ELSE
               EXIT
            END IF
         END DO
         lM%eAdj%prow(e+1) = j + 1
      END DO
      DEALLOCATE(adjList)

      RETURN
      END SUBROUTINE GETEADJ_MSH
!--------------------------------------------------------------------
!     Find element adjacency of a given face. Computes list of all
!     elements around a given element of a face.
      SUBROUTINE GETEADJ_FACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType),  INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjList(:,:),
     2   tmpI(:,:)

      ALLOCATE(incNd(lFa%nNo))
      incNd = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(tmpI(maxAdj, lFa%nNo))
      incNd = 0
      tmpI  = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
            tmpI(incNd(Ac), Ac) = e
         END DO
      END DO
      b = 2*maxAdj

 001  b = b + maxAdj
      DEALLOCATE(incNd)
      ALLOCATE(incNd(lFa%nEl))
      IF (ALLOCATED(adjList)) DEALLOCATE(adjList)
      ALLOCATE(adjList(b, lFa%nEl))
      adjList = 0
      incNd   = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO i=1, maxAdj
               IF (tmpI(i,Ac) .EQ. 0) EXIT
               flag = .TRUE.
               DO j=1, incNd(e)
                  IF (adjList(j,e) .EQ. tmpI(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(e) = incNd(e) + 1
                  IF (incNd(e) .GE. b) GOTO 001
                  adjList(incNd(e),e) = tmpI(i,Ac)
               END IF
            END DO
         END DO
      END DO
      maxAdj = MAXVAL(incNd)
      DEALLOCATE(tmpI, incNd)

      lFa%eAdj%nnz = 0
      DO e=1, lFa%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               lFa%eAdj%nnz = lFa%eAdj%nnz + 1
            ELSE
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE(lFa%eAdj%prow(lFa%nEl+1), lFa%eAdj%pcol(lFa%eAdj%nnz))
      j = 0
      lFa%eAdj%prow(1) = j + 1
      DO e=1, lFa%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               j = j + 1
               lFa%eAdj%pcol(j) = adjList(i,e)
            ELSE
               EXIT
            END IF
         END DO
         lFa%eAdj%prow(e+1) = j + 1
      END DO
      DEALLOCATE(adjList)

      RETURN
      END SUBROUTINE GETEADJ_FACE
!====================================================================
!     Finds an element of a mesh for any probe. Returns 0 if not found.
!     Uses Newton method to compute the parametric coordinate (xi) of
!     the probe with respect to an element for the given physical
!     coordinate (xp). Elements are searched prescribed by eList.
      SUBROUTINE FINDE(xp, lM, nNo, xg, ne, eList, Ec, xi, lDebug)
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, ne, eList(ne)
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd), xg(nsd,nNo)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
      REAL(KIND=RKIND), INTENT(OUT) :: xi(nsd)
      LOGICAL, INTENT(IN), OPTIONAL :: lDebug

      LOGICAL :: ldbg, l1, l2, l3, l4
      INTEGER(KIND=IKIND) :: a, e, i, Ac, eNoN
      REAL(KIND=RKIND) :: rt, xi0(nsd)

      LOGICAL, ALLOCATABLE :: eChck(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:)

      ldbg = .FALSE.
      IF (PRESENT(lDebug)) ldbg = lDebug
      Ec   = 0
      eNoN = lM%eNoN
      ALLOCATE(eChck(lM%nEl), xl(nsd,eNoN), N(eNoN), Nxi(nsd,eNoN))

!     Initialize parameteric coordinate for Newton's iterations
      xi0 = 0._RKIND
      DO i=1, lM%nG
         xi0 = xi0 + lM%xi(:,i)
      END DO
      xi0 = xi0 / REAL(lM%nG, KIND=RKIND)

      IF (ldbg) THEN
         WRITE(1000,'(A)') "=================================="
         WRITE(1000,'(A)',ADVANCE='NO') "Finding trace for xp: "
         DO i=1,  nsd
            WRITE(1000,'(A)',ADVANCE='NO') " "//STR(xp(i))
         END DO
         WRITE(1000,'(A)')
         WRITE(1000,'(A)')
         WRITE(1000,'(A)') "List of elems: <"//STR(ne)//">"
         DO e=1, ne
            WRITE(1000,'(A)',ADVANCE='NO') " "//STR(eList(e))
         END DO
         WRITE(1000,'(A)')
         WRITE(1000,'(A)') "=================================="
      END IF

      eChck = .FALSE.
      DO e=1, ne
         Ec = eList(e)
         IF (Ec .EQ. 0) CYCLE
         IF (eChck(Ec)) CYCLE
         eChck(Ec) = .TRUE.

         DO a=1, eNoN
            Ac = lM%IEN(a,Ec)
            xi(:) = xi(:) + lM%xi(:,a)
            xl(:,a) = xg(:,Ac)
         END DO

         IF (ldbg) THEN
            WRITE(1000,'(A)') "----------------------------"
            WRITE(1000,'(A)') "Probe el: "//STR(Ec)
            DO a=1, eNoN
               Ac = lM%IEN(a,Ec)
               WRITE(1000,'(4X,A)',ADVANCE='NO') STR(Ac)
               DO i=1, nsd
                  WRITE(1000,'(A)',ADVANCE='NO') " "//STR(xl(i,a))
               END DO
               WRITE(1000,'(A)')
            END DO
         END IF

         xi = xi0
         CALL GETXI(lM%eType, eNoN, xl, xp, xi, l1)

         IF (ldbg) THEN
            WRITE(1000,'(4X,A)',ADVANCE='NO') "xi: "
            DO i=1, nsd
               WRITE(1000,'(A)',ADVANCE='NO') " "//STR(xi(i))
            END DO
            WRITE(1000,'(A)')
         END IF

!        Check if parameteric coordinate is within bounds
         a = 0
         DO i=1, nsd
            IF (xi(i).GE.lM%xib(1,i) .AND. xi(i).LE.lM%xib(2,i))
     2         a = a + 1
         END DO
         l2 = a .EQ. nsd

!        Check for shape function even if the Newton's method returns OK
         CALL GETGNN(nsd, lM%eType, eNoN, xi, N, Nxi)

         IF (ldbg) THEN
            WRITE(1000,'(4X,A)',ADVANCE='NO') "N: "
            DO a=1, eNoN
               WRITE(1000,'(A)',ADVANCE='NO') " "//STR(N(a))
            END DO
            WRITE(1000,'(A)')
            CALL FLUSH(1000)
         END IF

!        Check if shape functions are within bounds and sum to unity
         i  = 0
         rt = 0._RKIND
         DO a=1, eNoN
            rt = rt + N(a)
            IF (N(a).GT.lM%Nb(1,a) .AND. N(a).LT.lM%Nb(2,a)) i = i + 1
         END DO
         l3 = i .EQ. eNoN
         l4 = rt.GE.0.9999_RKIND .AND. rt.LE.1.0001_RKIND

         l1 = ALL((/l1, l2, l3, l4/))
         IF (l1) RETURN
      END DO

      Ec = 0
      xi = 0._RKIND

      RETURN
      END SUBROUTINE FINDE
!====================================================================
      END MODULE ALLFUN
!====================================================================
