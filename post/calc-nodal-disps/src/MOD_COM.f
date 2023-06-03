!--------------------------------------------------------------------
!
!     All the data structures are defined in this module.
!
!--------------------------------------------------------------------

      MODULE COMMOD
      USE CHNLMOD

!     Include parameters
      INCLUDE "CONSTS.f"

!--------------------------------------------------------------------
!     Define all the derived data types

!     Mesh adjacency (neighboring element for each element)
      TYPE adjType
!        No of non-zeros
         INTEGER(KIND=IKIND) :: nnz = 0
!        Column pointer
         INTEGER(KIND=IKIND), ALLOCATABLE :: pcol(:)
!        Row pointer
         INTEGER(KIND=IKIND), ALLOCATABLE :: prow(:)
      END TYPE adjType

!     The face type containing mesh at boundary
      TYPE faceType
!        Parametric direction normal to this face (NURBS)
         INTEGER(KIND=IKIND) d
!        Number of nodes (control points) in a single element
         INTEGER(KIND=IKIND) eNoN
!        Element type
         INTEGER(KIND=IKIND) :: eType = eType_NA
!        The mesh index that this face belongs to
         INTEGER(KIND=IKIND) :: iM
!        Number of elements
         INTEGER(KIND=IKIND) :: nEl = 0
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        Number of nodes
         INTEGER(KIND=IKIND) :: nNo = 0
!        Global element Ids
         INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:)
!        Global node Ids
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
!        Global to local maping tnNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: lN(:)
!        Connectivity array
         INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:)
!        Surface area
         REAL(KIND=RKIND) area
!        Gauss point weights
         REAL(KIND=RKIND), ALLOCATABLE :: w(:)
!        Position coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!        Gauss points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Shape functions at Gauss points
         REAL(KIND=RKIND), ALLOCATABLE :: N(:,:)
!        Normal vector to each nodal point
         REAL(KIND=RKIND), ALLOCATABLE :: nV(:,:)
!        Shape functions derivative at Gauss points
         REAL(KIND=RKIND), ALLOCATABLE :: Nx(:,:,:)
!        Second derivatives of shape functions - for shells & IGA
         REAL(KIND=RKIND), ALLOCATABLE :: Nxx(:,:,:)
!        Face name for flux files
         CHARACTER(LEN=stdL) name
!        Face nodal adjacency
         TYPE(adjType) :: nAdj
!        Face element adjacency
         TYPE(adjType) :: eAdj
      END TYPE faceType

!     This is the container for a mesh
      TYPE mshType
!        Whether the shape function is linear
         LOGICAL lShpF
!        Element type
         INTEGER(KIND=IKIND) :: eType = eType_NA
!        Number of nodes in a single element
         INTEGER(KIND=IKIND) eNoN
!        Number of elements
         INTEGER(KIND=IKIND) :: nEl = 0
!        Number of faces
         INTEGER(KIND=IKIND) :: nFa = 0
!        Number of nodes
         INTEGER(KIND=IKIND) :: nNo = 0
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        The element type recognized by VTK format
         INTEGER(KIND=IKIND) vtkType
!        The connectivity array mapping eNoN,nEl
         INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:)
!        Gauss weights
         REAL(KIND=RKIND), ALLOCATABLE :: w(:)
!        Gauss integration points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Bounds on parameteric coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: xib(:,:)
!        Parent shape function
         REAL(KIND=RKIND), ALLOCATABLE :: N(:,:)
!        Shape function bounds
         REAL(KIND=RKIND), ALLOCATABLE :: Nb(:,:)
!        Parent shape functions gradient
         REAL(KIND=RKIND), ALLOCATABLE :: Nx(:,:,:)
!        Second derivatives of shape functions - used for shells & IGA
         REAL(KIND=RKIND), ALLOCATABLE :: Nxx(:,:,:)
!        Mesh Name
         CHARACTER(LEN=stdL) :: name
!        Mesh nodal adjacency
         TYPE(adjType) :: nAdj
!        Mesh element adjacency
         TYPE(adjType) :: eAdj
!        Faces are stored in this variable
         TYPE(faceType), ALLOCATABLE :: fa(:)
      END TYPE mshType

      TYPE traceType
         INTEGER :: gE
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:)
      END TYPE traceType
!--------------------------------------------------------------------
!     Declare variables to be read from input file

!     Number of spatial dimensions
      INTEGER(KIND=IKIND) nsd
!     Starting time step
      INTEGER(KIND=IKIND) startTS
!     Ending time step
      INTEGER(KIND=IKIND) endTS
!     Increment time step
      INTEGER(KIND=IKIND) incrTS
!     Number of probes
      INTEGER(KIND=IKIND) nprb

!     Time step size
      REAL(KIND=RKIND) dt

!     Saved output file name
      CHARACTER(LEN=stdL) saveName

!--------------------------------------------------------------------
!     Declare variables evaluated during runtime
!     Current time step
      INTEGER(KIND=IKIND) cTS
!     Number of meshes
      INTEGER(KIND=IKIND) nMsh

!     Time
      REAL(KIND=RKIND) time

!     Probe coordinates
      REAL(KIND=RKIND), ALLOCATABLE :: xprb(:,:)
!     Position vector
      REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!     Displacement vector
      REAL(KIND=RKIND), ALLOCATABLE :: Dn(:,:)

!     Mesh
      TYPE(mshType), ALLOCATABLE :: msh(:)

!     Traces of probes
      TYPE(traceType), ALLOCATABLE :: trace(:)

!     Input/output to the screen is handled by this structure
      TYPE(chnlType), POINTER :: std, err, wrn, dbg

!     To group above channels
      TYPE(ioType), TARGET :: io

      END MODULE COMMOD
!--------------------------------------------------------------------
