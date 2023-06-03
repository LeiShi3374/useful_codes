!**************************************************

      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      use vtkLegacyMod
      implicit none

      type meshType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: V(:)
      end type meshType

      integer :: nsd, nstart, nend, nfreq, nfiles, ntime
      real(kind=8) :: dt, time
      character(len=strL) :: fhdr

      integer, allocatable :: incNd(:)
      real(kind=8), allocatable :: acTime(:)

      type(meshType) :: msh

      contains
!-------------------------------------------------
         subroutine destroyMesh(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         if (allocated(lM%IEN)) deallocate(lM%IEN)
         if (allocated(lM%x))   deallocate(lM%x)
         if (allocated(lM%V))   deallocate(lM%V)

         return
         end subroutine destroyMesh
!-------------------------------------------------
      end module variables

!**************************************************

      program actvnTime
      use variables
      implicit none
      type(meshType) :: lM
      integer i, a
      character(len=strL) :: fname

      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

      call readInputs(fname)
      nfiles = (nend-nstart)/nfreq + 1
      do i=1, nfiles
         write(*,'(A)') repeat('*',48)
         ntime = nstart + (i-1)*nfreq
         time  = real(ntime,kind=8)*dt
         if (ntime .lt. 100) then
            write(fname,'(A,I3.3,A)') TRIM(fhdr), ntime, ".vtu"
         else
            write(fname,'(A)') TRIM(fhdr)//TRIM(STR(ntime))//".vtu"
         end if

         write(*,ftab1) "Reading numerical solution <--- "//trim(fname)
         call readVTU(msh, fname)
         if (i .eq. 1) then
            lM%nNo = msh%nNo
            lM%nEl = msh%nEl
            lM%eNoN = msh%eNoN
            allocate(lM%x(nsd,lM%nNo), lM%IEN(lM%eNoN,lM%nEl))
            lM%x   = msh%x
            lM%IEN = msh%IEN

            allocate(incNd(msh%nNo), acTime(msh%nNo))
            incNd  = 0
            acTime = 0D0
         end if

         do a=1, msh%nNo
            if (incNd(a) .eq. 1) cycle
            if (msh%V(a) .gt. -eps) then
               incNd(a)  = 1
               acTime(a) = time
            end if
         end do
         call destroyMesh(msh)

         if (SUM(incNd) .eq. msh%nNo) then
            if (i .lt. nfiles) then
               write(stdout,ftab1) "All nodes are activated"
               exit
            else
               write(stdout,ftab1) "Warning: all nodes are not yet "// &
     &            "activated"
            end if
         end if
      end do

      allocate(lM%V(lM%nNo))
      lM%V = acTime

      write(*,'(A)') repeat('*',48)
      write(fname,'(A)') "activationTime.vtu"
      write(stdout,ftab1) "Writing output file ---> "//trim(fName)
      call writeVTU(lM, fname)
      write(*,'(A)') repeat('*',48)

      call destroyMesh(lM)
      deallocate(incNd, acTime)

      return
      end program actvnTime

!**************************************************

      subroutine readInputs(fname)
      use variables
      implicit none
      character(len=strL), intent(in) :: fname
      integer fid

      fid = 1265
      open(fid,file=trim(fname))
      read(fid,*) ! #Num spatial dimensions
      read(fid,*) nsd
      read(fid,*) ! empty line

      read(fid,*) ! #nStart, nEnd, nFreq
      read(fid,*) nstart, nend, nfreq
      read(fid,*) ! empty line

      read(fid,*) ! #time stepping
      read(fid,*) dt
      read(fid,*) ! empty line

      read(fid,*) ! #vtu files path
      read(fid,'(A)') fhdr
      read(fid,*) ! empty line
      close(fid)

      end subroutine readInputs

!**************************************************

      subroutine readVTU(lM, fname)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtu
      integer :: istat, i
      real(kind=8), allocatable :: tmpX(:,:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtu, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtu, lM%nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtu, lM%nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtu, lM%eNoN, istat)
         if (istat .lt. 0) exit

         allocate(lM%ien(lM%eNoN,lM%nEl), lM%x(nsd,lM%nNo))
         allocate(tmpX(maxNSD,lM%nNo))
         call getVTK_pointCoords(vtu, tmpX, istat)
         if (istat .lt. 0) exit
         lM%x(:,:) = tmpX(1:nsd,:)

         call getVTK_elemIEN(vtu, lM%ien, istat)
         if (istat .lt. 0) exit
         lM%IEN(:,:) = lM%IEN(:,:) + 1

         allocate(lM%V(lM%nNo))
         call getVTK_pointData(vtu, "EP_Action_potential", lM%V, istat)
         if (istat .lt. 0) exit

         deallocate(tmpX)
         call flushVTK(vtu)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk file read error"
         STOP
      end if

      return
      end subroutine readVTU

!**************************************************

      subroutine writeVTU(lM, fname)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtu
      integer :: istat

      if (nsd .EQ. 3) then
         select case (lM%eNoN)
         case(4)
            lM%vtkType = 10  ! Tet !
         case(6)
            lM%vtkType = 13  ! Wedge !
         case(8)
            lM%vtkType = 12  ! Brick/Hex !
         case default
            write(stdout,ftab2) "Unknown vtk type"
         end select
      else
         select case (lM%eNoN)
         case(3)
            lM%vtkType = 5   ! Tri !
         case(4)
            lM%vtkType = 9   ! Bilinear/Quad !
         case(9)
            lM%vtkType = 28  ! Biquad/Quadr-Quad !
         case default
            write(stdout,ftab2) "Unknown vtk type"
         end select
      end if

      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtu, fname, istat)
         if (istat .lt. 0) exit

         call putVTK_pointCoords(vtu, lM%x, istat)
         if (istat .lt. 0) exit

         lM%IEN = lM%IEN - 1
         call putVTK_elemIEN(vtu, lM%IEN, lM%vtkType, istat)
         if (istat .lt. 0) exit
         lM%IEN = lM%IEN + 1

         call putVTK_pointData(vtu,"Activation_Time", lM%V, istat)
         if (istat .lt. 0) exit

         call vtkWriteToFile(vtu, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtu)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk write error"
         STOP
      end if

      return
      end subroutine writeVTU

!**************************************************
