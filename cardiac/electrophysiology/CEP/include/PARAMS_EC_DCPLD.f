!--------------------------------------------------------------------
!
!     Parameters for decoupled and uniformly activated excitation model
!     for excitation-contraction coupling. Parameters are chosen based
!     on below references.
!
!     Reference for active stress model:
!        Pfaller, M. R., et al.(2019). The importance of the pericardium
!        for cardiac biomechanics: from physiology to computational
!        modeling. Biomechanics and Modeling in Mechanobiology,
!        18(2), 503–529. https://doi.org/10.1007/s10237-018-1098-4
!
!     Reference for active strain model:
!        Barbarotta, L., et al.(2018). A transmurally heterogeneous
!        orthotropic activation model for ventricular contraction and
!        its numerical validation. International Journal for Numerical
!        Methods in Biomedical Engineering, 34(12), 1–24.
!        https://doi.org/10.1002/cnm.3137
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     Cardiac cycle length (ms)
      REAL(KIND=RKIND) :: t_CL = 1000._RKIND
!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Contractility (Pa)
      REAL(KIND=RKIND) :: sigm0 = 9.0E4_RKIND

!     Min activation rate (ms^-1)
      REAL(KIND=RKIND) :: min_alpha = -0.03_RKIND

!     Max activation rate (ms^-1)
      REAL(KIND=RKIND) :: max_alpha = 0.005_RKIND

!     Onset of systole (ms)
      REAL(KIND=RKIND) :: t_sys = 170._RKIND

!     Onset of diastole (ms)
      REAL(KIND=RKIND) :: t_dia = 484._RKIND

!     Gamma (ms)
      REAL(KIND=RKIND) :: gamma = 5._RKIND
!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active strain model
!     Min fiber shortening
      REAL(KIND=RKIND) :: gf_min = -0.13_RKIND

!     Initial time for activation (ms)
      REAL(KIND=RKIND) :: ta_s = 170._RKIND

!     Duration of activation (ms)
      REAL(KIND=RKIND) :: ta_e = 480._RKIND
!--------------------------------------------------------------------

