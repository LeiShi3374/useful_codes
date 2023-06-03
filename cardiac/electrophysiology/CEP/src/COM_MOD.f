!-----------------------------------------------------------------------
!
!     This module defines common and other auxiliary data structures for
!     cardiac electrophysiology model equation.
!
!-----------------------------------------------------------------------

      MODULE COMMOD
      USE CHNLMOD
      USE CEPMOD
      IMPLICIT NONE

!     Input file
      CHARACTER(LEN=stdL) :: cep_fIn

!     Display progress
      LOGICAL :: iProg = .TRUE.

!     No. of time steps
      INTEGER(KIND=IKIND) nTS

!     Time step increment
      REAL(KIND=RKIND) :: dt = 0._RKIND

!     Output log file
      CHARACTER(LEN=stdL) :: oFile

!     Cardiac electrophysiology type
      TYPE(cepModelType) :: cep

!     Communicator channel and its pointers
      TYPE(ioType), TARGET :: io
      TYPE(chnlType), POINTER :: std, wrn, err, dbg

      END MODULE COMMOD
!#######################################################################
