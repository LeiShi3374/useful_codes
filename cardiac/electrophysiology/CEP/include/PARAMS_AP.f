!--------------------------------------------------------------------
!
!     Parameters for Aliev-Panfilov cellular activation model.
!     Parameters are chosen based on below references.
!
!     Reference for Aliev-Panfilov electrophysiology model:
!        Goktepe, S., & Kuhl, E. (2009). Computational modeling of
!        cardiac electrophysiology: A novel finite element approach.
!        Int. J. Numer. Meth. Engng, 79, 156–178.
!        https://doi.org/10.1002/nme
!
!     Reference for active stress model:
!        Goktepe, S., & Kuhl, E. (2010). Electromechanics of the heart:
!        A unified approach to the strongly coupled excitation-
!        contraction problem. Computational Mechanics, 45(2–3), 227–243.
!        https://doi.org/10.1007/s00466-009-0434-z
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling (mV)
      REAL(KIND=RKIND) :: Vscale  = 100._RKIND
!     Time scaling (ms)
      REAL(KIND=RKIND) :: Tscale  = 12.9_RKIND
!     Voltage offset parameter (mV)
      REAL(KIND=RKIND) :: Voffset = -80._RKIND
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=RKIND) :: alpha = 1.E-2_RKIND
      REAL(KIND=RKIND) :: a     = 2.E-3_RKIND
      REAL(KIND=RKIND) :: b     = 0.15_RKIND
      REAL(KIND=RKIND) :: c     = 8._RKIND
      REAL(KIND=RKIND) :: mu1   = 0.2_RKIND
      REAL(KIND=RKIND) :: mu2   = 0.3_RKIND
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Resting voltage (mV)
      REAL(KIND=RKIND) :: Vrest = -80._RKIND
!     Critical voltage (mV)
      REAL(KIND=RKIND) :: Vcrit = -30._RKIND
!     Saturation potential (Pa/mV)
      REAL(KIND=RKIND) :: K_T = 5.0E3_RKIND
!     Minimum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_0 = 0.1_RKIND
!     Maximum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_i = 1._RKIND
!     Transition rate (mV^{-1})
      REAL(KIND=RKIND) :: xi_T  = 1._RKIND

!     Cm: Cell capacitance per unit surface area
      REAL(KIND=RKIND) :: Cm  = 1._RKIND
!     sV: Surface to volume ratio
      REAL(KIND=RKIND) :: sV  = 1._RKIND
!     rho: Cellular resistivity
      REAL(KIND=RKIND) :: rho = 1._RKIND
!#######################################################################
