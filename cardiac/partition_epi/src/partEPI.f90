!**************************************************

      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      implicit none

      integer, parameter :: nsd = 3

      type faceType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: vtkType
         integer, allocatable :: gN(:)
         integer, allocatable :: gE(:)
         integer, allocatable :: eID(:)
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: x(:,:)
         character(len=strL) :: fname
      end type faceType

      type(faceType), allocatable :: sfa(:), tfa(:)

      integer :: nfa
      character(len=strL) :: fout

      end module variables

!**************************************************

      program partEPI
      use variables
      implicit none
      character(len=strL) :: fname
      type(faceType) :: lfa

      if (iargc() .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1, fname)

      write(stdout,ftab1)
      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)

      call readInputs(fname, lfa)

      call partFace(lFa)

      write(stdout,ftab1)
      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)

      end program partEPI

!**************************************************

      subroutine readInputs(fin, lfa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lfa
      character(len=strL), intent(in) :: fin

      integer ifa, fid
      character(len=strL) :: fname

      fid = 10
      write(stdout,ftab1) "Reading inputs   <---   "//trim(fin)
      open(fid,file=trim(fin))
      read(fid,*) ! num. partitions
      read(fid,*) nfa

      allocate(sfa(nfa), tfa(nfa))
      do ifa=1, nfa
         read(fid,*)
         read(fid,*) ! source partition file
         read(fid,'(A)') fname
         call readVTP(sfa(ifa), fname)
         write(stdout,ftab2) "Loading face   <---   "//trim(fname)
      end do
      read(fid,*)
      read(fid,*) ! target file to be partitioned
      read(fid,'(A)') fname
      call readVTP(lfa, fname)
      write(stdout,ftab2) "Loading face   <---   "//trim(fname)
      close(fid)

      fid  = scan(trim(fname), '/', back=.true.)
      write(fout,'(A)') fname(fid+1:len(trim(fname))-4)//"_part.vtp"

      return
      end subroutine readInputs

!**************************************************

      subroutine readVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: vtkType, istat, iopt
      integer :: e, nNo, eNoN, nEl
      real(kind=8), allocatable :: tmpX(:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtp, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtp, nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtp, nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtp, eNoN, istat)
         if (istat .lt. 0) exit

         lFa%nNo = nNo
         lFa%nEl = nEl
         lFa%eNon = eNoN
         allocate(lFa%ien(eNoN,nEl), lFa%x(nsd,nNo), lFa%gN(nNo), &
               lFa%gE(nEl))

         call getVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call getVTK_elemIEN(vtp, lFa%ien, istat)
         if (istat .lt. 0) exit

         call getVTK_pointData(vtp, 'GlobalNodeID', lFa%gN, istat)
         if (istat .lt. 0) exit

         call getVTK_elemData(vtp, 'GlobalElementID', lFa%gE, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "VTP file read error. Missing data"
         STOP
      end if

      return
      end subroutine readVTP

!**************************************************

      subroutine partFace(lfa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lfa

      write(stdout,ftab1) "Partitioning face.. "
      call setElemID(lfa)

      write(stdout,ftab2) "Writing face   --->   "//trim(fout)
      call writeVTP(lFa, fout)

      end subroutine partFace

!**************************************************

      subroutine setElemID(lfa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lfa

      integer a, b, e, f, Ac, Bc, ifa
      real(kind=8) :: smin, ds, x1(nsd), x2(nsd)

      write(stdout,ftab2) "Setting element IDs on face.."
      allocate(lfa%eID(lfa%nel))
      lfa%eID = 0
      do e=1, lfa%nel
         x1 = 0D0
         do a=1, lfa%eNoN
            Ac = lfa%IEN(a,e) + 1
            x1 = x1 + lfa%x(:,Ac)
         end do
         x1 = x1 / real(lfa%eNoN,kind=8)

         smin = huge(smin)
         do ifa=1, nfa
            do f=1, sfa(ifa)%nel
               x2 = 0D0
               do b=1, sfa(ifa)%eNoN
                  Bc = sfa(ifa)%IEN(b,f) + 1
                  x2 = x2 + sfa(ifa)%x(:,Bc)
               end do
               x2 = x2/real(sfa(ifa)%eNoN,kind=8)

               ds = SQRT(SUM( (x1(:)-x2(:))**2 ))
               if (smin .gt. ds) then
                  smin = ds
                  lfa%eID(e) = ifa
               end if
            end do
         end do
      end do

      return
      end subroutine setElemID

!**************************************************

      subroutine writeVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat, vtkType

      if (nsd .EQ. 3) then
         select case (lFa%eNoN)
         case(3)
            vtkType = 5   ! Tri !
         case(4)
            vtkType = 9   ! Quad !
         end select
      else
         select case (lFa%eNoN)
         case(2)
            vtkType = 3   ! Line !
         case(3)
            vtkType = 21   ! Quadr-Edge !
         end select
      end if

      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtp, fname, istat)
         if (istat .lt. 0) exit

         call putVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call putVTK_elemIEN(vtp, lFa%IEN, vtkType, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"GlobalNodeID", lFa%gN, istat)
         if (istat .lt. 0) exit

         call putVTK_elemData(vtp,"GlobalElementID", lFa%gE, istat)
         if (istat .lt. 0) exit

         if (allocated(lfa%eID)) then
            call putVTK_elemData(vtp,"PART", lFa%eID, istat)
            if (istat .lt. 0) exit
         end if

         call vtkWriteToFile(vtp, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk write error"
         STOP
      end if

      return
      end subroutine writeVTP

!**************************************************
