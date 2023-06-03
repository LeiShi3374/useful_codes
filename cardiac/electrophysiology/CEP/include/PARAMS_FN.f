!--------------------------------------------------------------------
!
!     Parameters used for Fitzhugh-Nagumo Myocyte Activation Model.
!
!     Reference for Aliev-Panfilov electrophysiology model:
!        Goktepe, S., & Kuhl, E. (2009). Computational modeling of
!        cardiac electrophysiology: A novel finite element approach.
!        Int. J. Numer. Meth. Engng, 79, 156–178.
!        https://doi.org/10.1002/nme
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=RKIND) :: alpha = -0.5_RKIND
      REAL(KIND=RKIND) :: a = 0._RKIND
      REAL(KIND=RKIND) :: b = -0.6_RKIND
      REAL(KIND=RKIND) :: c = 50._RKIND
!#######################################################################

