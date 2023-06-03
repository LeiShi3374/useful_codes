!--------------------------------------------------------------------
!
!     Parameter definitions
!
!--------------------------------------------------------------------

!     Maximum possible nsd
      INTEGER(KIND=IKIND), PARAMETER :: maxnsd = 3
!--------------------------------------------------------------------
!     Gauss point coordinates, upto 5 points
      REAL(KIND=RKIND), PARAMETER :: gW(5,5) =
     2   RESHAPE( (/2._RKIND, 0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     3   1._RKIND, 1._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     4   0.5555555555555556_RKIND, 0.8888888888888889_RKIND,
     4   0.5555555555555556_RKIND, 0._RKIND, 0._RKIND,
     5   0.3478548451374538_RKIND, 0.6521451548625462_RKIND,
     5   0.6521451548625462_RKIND, 0.3478548451374538_RKIND, 0._RKIND,
     6   0.2369268850561890_RKIND, 0.4786286704993665_RKIND,
     6   0.5688888888888889_RKIND, 0.4786286704993665_RKIND,
     6   0.2369268850561890_RKIND/) , (/5, 5/) )
!--------------------------------------------------------------------
!     Gauss point weights, upto 5 points
      REAL(KIND=RKIND), PARAMETER :: gXi(5,5) =
     2   RESHAPE( (/0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     3   -0.57735026918962584_RKIND, 0.57735026918962584_RKIND,
     3   0._RKIND, 0._RKIND, 0._RKIND,
     4   -0.7745966692414834_RKIND, 0._RKIND, 0.7745966692414834_RKIND,
     4   0._RKIND, 0._RKIND,
     5   -0.86113631159405257_RKIND, -0.33998104358485631_RKIND,
     5   0.33998104358485631_RKIND, 0.86113631159405257_RKIND, 0._RKIND,
     6   -0.90617984593866396_RKIND, -0.53846931010568311_RKIND,
     6   0._RKIND, 0.53846931010568311_RKIND,
     6   0.90617984593866396_RKIND/) , (/5, 5/) )
!--------------------------------------------------------------------
      CHARACTER, PARAMETER :: delimiter = "/"
!--------------------------------------------------------------------
!     Types of accepted elements
!     Point, Line (linear), Line (quadratic), Triangle (linear),
!     Triangle (quadratic), Quads (bilinear), Quads (serendipity),
!     Quads (biquadratic), Tetrahedron (linear), Tets (quadratic),
!     Hexgonal bricks (trilinear), Hex (quadratic/serendipity),
!     Hex (triquadratic), Wedge, NURBS
      INTEGER(KIND=IKIND), PARAMETER :: eType_NA = 100, eType_PNT = 101,
     2   eType_LIN1 = 102, eType_LIN2 = 103, eType_TRI3 = 104,
     3   eType_TRI6 = 105, eType_QUD4 = 106, eType_QUD8 = 107,
     4   eType_QUD9 = 108, eType_TET4 = 109, eType_TET10 = 110,
     5   eType_HEX8 = 111, eType_HEX20 = 112, eType_HEX27 = 113,
     6   eType_WDG = 114
!--------------------------------------------------------------------
