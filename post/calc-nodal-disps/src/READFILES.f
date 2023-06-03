!====================================================================
!
!
!
!====================================================================
      SUBROUTINE READFILES()
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE

      LOGICAL flag
      INTEGER i
      CHARACTER(LEN=stdL) :: fname, ctmp

      TYPE(fileType) :: ftmp
      TYPE(listType) :: list
      TYPE(listType), POINTER :: lPtr

!     Set the communicator channels
      CALL io%o%new(CHNL_O, tag="", oTS=.TRUE., oTF=.FALSE.)
      CALL io%e%new(CHNL_E, tag="", oTF=.FALSE.)
      CALL io%w%new(CHNL_W, tag="", oTF=.FALSE.)
      CALL io%d%new(CHNL_D, tag="", oTF=.FALSE.)

      std => io%o
      err => io%e
      wrn => io%w
      dbg => io%d

!     Set the default values
      saveName = ""
      startTS  = 1
      endTS    = 0
      incrTS   = 1

      i = IARGC()
      IF (i .NE. 0) THEN
         IF (i .GT. 1) err = " Too many arguments"
         CALL GETARG(1,ctmp)
         fname = ctmp
      END IF

!     Load all the LPN parameters from the file into list structure
      std = " Reading input file "//CLR(fname)
      list = listType(fname, io)

!     Copy read values into local variables
      lPtr => list%get(nsd, "Number of spatial dimensions", 1,ll=2,ul=3)
      lPtr => list%get(dt, "Time step size", 1, lb=0._RKIND)

      lPtr => list%get(saveName, "Saved results folder path", 1)
      saveName = TRIM(saveName)//delimiter

      ctmp = "result"
      lPtr => list%get(ctmp, "Name prefix of saved VTU files")
      saveName = TRIM(saveName)//ctmp

      lPtr => list%get(startTS, "VTU files start time step", 1)
      lPtr => list%get(endTS,   "VTU files end time step", 1)
      lPtr => list%get(incrTS,  "VTU files increment")

!     Read mesh
      nMsh = 1
      ALLOCATE(msh(nMsh))
      lPtr => list%get(ftmp, "Mesh file path", 1)
      CALL READVTU(msh(1), ftmp%fname)
      CALL SELECTELE(msh(1))

!     Read probe coordinates
      nprb = list%srch("Probe coordinates", ll=1)
      ALLOCATE(xprb(nsd,nprb), trace(nprb))
      xprb = 0._RKIND
      DO i=1, nprb
         lPtr => list%get(xprb(:,i), "Probe coordinates", i)
         ALLOCATE(trace(i)%xi(nsd))
      END DO

!     Destroy list once all the data is copied
      CALL DESTROYLIST(list)

      RETURN
      END SUBROUTINE READFILES
!====================================================================
      SUBROUTINE READVTU(lM, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND) :: iStat
      TYPE(vtkXMLType) :: vtu

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (init)"

      CALL getVTK_numPoints(vtu, lM%nNo, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num points)"

      CALL getVTK_numElems(vtu, lM%nEl, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num cells)"

      CALL getVTK_nodesPerElem(vtu, lM%eNoN, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (nodes per cell)"

      ALLOCATE(x(nsd,lM%nNo), tmpX(maxNSD,lM%nNo))
      CALL getVTK_pointCoords(vtu, tmpX, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (coords)"
      x(:,:) = tmpX(1:nsd,:)
      DEALLOCATE(tmpX)

      ALLOCATE(lM%IEN(lM%eNoN,lM%nEl))
      CALL getVTK_elemIEN(vtu, lM%IEN, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (ien)"
      lM%IEN = lM%IEN + 1

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READVTU
!====================================================================
      SUBROUTINE READVTUDISP(fName, nNo, lD)
      USE COMMOD
      USE vtkXMLMod
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: fName
      INTEGER, INTENT(IN) :: nNo
      REAL(KIND=RKIND), INTENT(INOUT) :: lD(nsd,nNo)

      INTEGER(KIND=IKIND) :: iStat, i
      TYPE(vtkXMLType) :: vtu

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (init)"

      CALL getVTK_numPoints(vtu, i, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num points)"
      IF (i .NE. nNo) THEN
         err = "Incompatible mesh and results file "//TRIM(fName)
      END IF

      ALLOCATE(tmpX(maxNSD,nNo))
      CALL getVTK_pointData(vtu, "Displacement", tmpX, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (Displacement)"
      lD(:,:) = tmpX(1:nsd,:)
      DEALLOCATE(tmpX)

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READVTUDISP
!====================================================================
      SUBROUTINE FINDPRBTRACE(lM, xp, ltrc)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd)
      TYPE(traceType), INTENT(INOUT) :: ltrc

      LOGICAL ldbg
      INTEGER e, Ec
      REAL(KIND=RKIND) :: xi(nsd)

      INTEGER, ALLOCATABLE :: eList(:)

      ldbg  = .FALSE.

      ALLOCATE(eList(msh(1)%nEl))
      DO e=1, lM%nEl
         eList(e) = e
      END DO

      ltrc%gE = 0
      CALL FINDE(xp, lM, lM%nNo, x, lM%nEl, eList, Ec, ltrc%xi, ldbg)
      IF (Ec .EQ. 0) err = " Probe could be outside the domain"
      ltrc%gE = Ec

      RETURN
      END SUBROUTINE FINDPRBTRACE
!====================================================================
