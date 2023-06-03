!**************************************************

      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      implicit none

      integer, parameter :: nd=3, nEd=2

      type faceType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer, allocatable :: gN(:)
         integer, allocatable :: gE(:)
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: d(:,:)
         character(len=strl) :: name
         character(len=strl), allocatable :: edge_name(:)
      end type faceType

      type(faceType), allocatable :: fa(:), ed(:)

      integer :: nFa
      integer, allocatable :: lmrks(:,:)
      character(len=strL) :: fdir

      integer :: writeSingleFile

      end module variables

!**************************************************

      program main
      use variables
      implicit none

      integer :: a, i, iFa, iEd
      character(len=strL) :: fname, cmd

      if (iargc() .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         stop
      end if
      call getarg(1,fname)

      write(*,ftab1) "------------------------------------------"
      write(*,ftab1) "Reading inputs.."
      call readInputs(fname)

      write(*,ftab1) "Reading face data.."
      do iFa=1, nFa
         write(*,ftab2) "<"//trim(fa(iFa)%name)//">"
         call readVTU(fa(iFa))
         fa(iFa)%IEN = fa(iFa)%IEN + 1
      end do

      write(cmd,'(A)') "mkdir -p "//trim(fdir)
      call system(trim(cmd))

      if (writeSingleFile .ne. 0) then
         do iEd=1, nEd
            write(fname,'(A)') trim(fdir)//"/edge_"//trim(str(iEd))//&
               ".txt"
            open(100+iEd,file=trim(fname))
            write(100+iEd,'(A)') "VARIABLES = X, Y, Z, Dx, DY, DZ"//&
               ", INDX"
         end do
      end if

      do iFa=1, nFa
         call getEdges(fa(iFa), ed, lmrks(:,iFa))

         do iEd=1, nEd       
            i = len(trim(fa(iFa)%name))
            write(fname,'(A)') trim(fdir)//"/"//trim(fa(iFa)%edge_name(iEd))//".vtp"
            ed(iEd)%IEN = ed(iEd)%IEN - 1
            call writeVTP(ed(iEd), fname)
            ed(iEd)%IEN = ed(iEd)%IEN + 1

            if (writeSingleFile .ne. 0) then
               do a=1, ed(iEd)%nNo
                  do i=1, nd
                     write(100+iEd,'(F12.6,X)',advance='no') &
                        ed(iEd)%x(i,a)
                  end do
                  do i=1, nd
                     write(100+iEd,'(F12.6,X)',advance='no') &
                        ed(iEd)%d(i,a)
                  end do
                  write(100+iEd,'(2X,I2)') iFa
               end do
            end if
         end do
      end do

      if (writeSingleFile .ne. 0) then
         do iEd=1, nEd
            close(100+iEd)
         end do
      end if

      return
      end program main

!**************************************************

      subroutine readInputs(fname)
      use variables
      implicit none

      character(len=strl), intent(in) :: fname

      character(len=strl) :: aaa
      integer, parameter :: fid = 10
      integer :: iFa, iEd

      writeSingleFile = 0
      open(fid,file=trim(fname))
      read(fid,*) ! num leaflets !
      read(fid,*) nFa

      allocate(fa(nFa), ed(nEd), lmrks(nEd,nFa))

      read(fid,*)
      do iFa=1, nFa
         read(fid,'(A)') fa(iFa)%name
      end do

      read(fid,*)
      do iFa=1, nFa
         read(fid,*) (lmrks(iEd,iFa),iEd=1, 2)
      end do

      read(fid,*)
      do iFa=1, nFa
         allocate(fa(iFa)%edge_name(nEd))
         do iEd = 1, nEd
            read(fid,'(A)') fa(iFa)%edge_name(iEd)
         end do
      end do
      
      read(fid,*)
      read(fid,'(A)') fdir
      read(fid,*)
      read(fid,*) writeSingleFile
      close(fid)

      return
      end subroutine readInputs

 !**************************************************

      subroutine readVTU(lFa)
      use variables
      implicit none

      type(faceType), intent(inout) :: lFa

      type(vtkXMLType) :: vtu
      integer :: vtkType, istat
      integer :: e, nNo, eNoN, nEl
      real(kind=8), allocatable :: tmpX(:,:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtu, trim(lFa%name), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtu, lFa%nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtu, lFa%nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtu, lFa%eNoN, istat)
         if (istat .lt. 0) exit

         allocate(lFa%ien(lFa%eNoN,lFa%nEl), lFa%x(nd,lFa%nNo))

         call getVTK_pointCoords(vtu, lFa%x, istat)
         if (istat .lt. 0) exit

         if (allocated(tmpX)) deallocate(tmpX)
         allocate(tmpX(maxNSD,lFa%nNo))
         call getVTK_pointData(vtu, "Displacement", tmpX, istat)
         if (istat .ge. 0) then
            write(*,ftab3) "Loading displacement field.."
            allocate(lFa%d(nd, lFa%nNo))
            lFa%d(:,:) = tmpX(1:nd,:)
            deallocate(tmpX)
         else
            write(*,ftab4) "Ignoring.."
            istat = 0
         end if
 
         call getVTK_elemIEN(vtu, lFa%ien, istat)
         if (istat .lt. 0) exit

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

      subroutine getEdges(lFa, lEd, refNd)
      use variables
      implicit none

      type(faceType), intent(in) :: lFa
      type(faceType), intent(out) :: lEd(nEd)
      integer :: refNd(nEd)

      logical :: flag, dflag
      integer :: a, b, e, Ac, Bc
      integer :: a1, b1, e1, Ac1, Bc1
      integer :: i, j, nNb, nEb, maxNEb
      integer, allocatable :: IENb(:,:), bNds(:), tmpI(:,:), ePtr(:)
  
!     Collect all the boundary edges of the face
      maxNEb = 250
      nEb = 0
 100  maxNEb = maxNEb*4
      if (allocated(IENb)) deallocate(IENb, ePtr)
      allocate(IENb(3,maxNEb))
      IENb = 0
      do e=1, lFa%nEl
         do a=1, lFa%eNoN
            b = a + 1
            if (a .eq. lFa%eNoN) b = 1
            Ac = lFa%IEN(a,e)
            Bc = lFa%IEN(b,e)
            flag = .true.
            do e1=1, lFa%nEl
               if (e1 .eq. e) cycle
               do a1=1, lFa%eNoN
                  b1 = a1 + 1
                  if (a1 .eq. lFa%eNoN) b1 = 1
                  Ac1 = lFa%IEN(a1,e1)
                  Bc1 = lFa%IEN(b1,e1)
                  i = 0
                  if (Ac.eq.Ac1 .or. Ac.eq.Bc1) i = i + 1
                  if (Bc.eq.Ac1 .or. Bc.eq.Bc1) i = i + 1
                  if (i .eq. 2) then
                     flag = .false.
                     exit
                  end if
               end do
               if (.not.flag) exit
            end do
            if (flag) then
               nEb = nEb + 1
               if (nEb .gt. maxNEb) goto 100
               IENb(1,nEb) = Ac
               IENb(2,nEb) = Bc
               IENb(3,nEb) = e
            end if
         end do
      end do

      allocate(tmpI(3,nEb))
      tmpI = IENb(:,1:nEb)
      deallocate(IENb)
      allocate(IENb(3,nEb))
      IENb = tmpI
      deallocate(tmpI)   

      do i=1, nEd
   !     Choose a seed as one of the landmarks for sorting the border nodes
         flag = .false.
         do e=1, nEb
            Ac = IENb(1,e)
            Bc = IENb(2,e)
            if (refNd(i) .eq. Ac) then
               flag = .true.
               exit
            end if
         end do
         if (.not.flag) write(*,ftab4) &
            "Error: choosen landmark not on the boundary"

         allocate(bNds(nEb), ePtr(nEb))
         bNds = 0
         ePtr = 0
         nNb  = 1
         bNds(nNb) = Ac
         ePtr(nNb) = IENb(3,e)
         do
            do e=1, nEb
               Ac = IENb(1,e)
               if (Bc .eq. Ac) then
                  Bc  = IENb(2,e)
                  nNb = nNb + 1
                  bNds(nNb) = Ac
                  ePtr(nNb) = IENb(3,e)
                  exit
               end if
            end do
            if (Bc .eq. bNds(1)) exit
         end do

         lEd(i)%nNo = nNb
         lEd(i)%nEl = lEd(i)%nNo

         allocate(lEd(i)%x(nd,lEd(i)%nNo), lEd(i)%gN(lEd(i)%nNo), &
            lEd(i)%IEN(2,lEd(i)%nEl), lEd(i)%gE(lEd(i)%nEl))
         if (dflag) allocate(lEd(i)%d(nd,lEd(i)%nNo))

         do j=1, lEd(i)%nNo
            Ac = bNds(j)
            lEd(i)%x(:,j) = lFa%x(:,Ac)
            lEd(i)%gN(j)  = Ac
            if (dflag) lEd(i)%d(:,j) = lFa%d(:,Ac)

            lEd(i)%gE(j) = ePtr(j)

         end do

!        IEN array
         do e=1, lEd(i)%nEl - 1
            lEd(i)%IEN(1,e) = e
            lEd(i)%IEN(2,e) = e+1
            ! print *, lEd(i)%IEN(:,e)
         end do
         e = lEd(i)%nEl
         lEd(i)%IEN(1,e) = e
         lEd(i)%IEN(2,e) = 1

         ! print *, lEd(i)%gE(:)

         deallocate(bNds, ePtr)
      end do   

      return
      end subroutine getEdges

!**************************************************

      subroutine writeVTP(lFa, fname)
      use variables
      implicit none

      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat, vtkType

      select case (lFa%eNoN)
      case(2)
         vtkType = 3   ! Line !
      case(3)
         vtkType = 5   ! Tri !
      case(4)
         vtkType = 9   ! Quad !
      end select

      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtp, fname, istat)
         if (istat .lt. 0) exit

         call putVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, istat)
         if (istat .lt. 0) exit

         if (allocated(lFa%d)) then
            call putVTK_pointData(vtp, "Displacement", lFa%d, istat)
            if(istat .lt. 0) exit
         end if

         call putVTK_elemIEN(vtp, lFa%IEN, vtkType, istat)
         if (istat .lt. 0) exit

         call putVTK_elemData(vtp, "GlobalElementID", lFa%gE, istat)
         if (istat .lt. 0) exit

         call vtkWriteToFile(vtp, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk write error"
         stop
      end if

      return
      end subroutine writeVTP

!**************************************************
