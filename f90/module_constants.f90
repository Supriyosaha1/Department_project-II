module module_constants
  
  public

  ! from the NIST Constant Index
  real(kind=8),parameter :: clight         = 2.99792458d+10          ![cm/s] light speed
  real(kind=8),parameter :: planck         = 6.62607015d-27          ![erg s] Planck's constant
  real(kind=8),parameter :: kb             = 1.380 649d-16           ![erg/K] Boltzman constant
  real(kind=8),parameter :: mp             = 1.67262192369d-24       ![g] proton mass
  real(kind=8),parameter :: me             = 9.1093837015d-28        ![g] electron 
  real(kind=8),parameter :: amu            = 1.66053906660d-24       ![g] atomic mass unit
  real(kind=8),parameter :: e_ch           = 4.80320471d-10          ![esu] electron charge
  real(kind=8),parameter :: cmtoA          = 1.0d8                   ! from cm to A
  real(kind=8),parameter :: evtoerg        = 1.602176634d-12         ! from eV to erg
  real(kind=8),parameter :: grtoev         = 1.782661d-33*clight**2  ! from gr to eV
  real(kind=8),parameter :: pi             = 3.1415926535898d0       ! pi
  real(kind=8),parameter :: twopi          = 2.0d0 * pi              ! 2 x pi
  real(kind=8),parameter :: fourpi         = 4.0d0 * pi              ! 4 x pi 
  real(kind=8),parameter :: sqrtpi         = 1.77245387556703d0      ! sqrt(pi)
  real(kind=8),parameter :: XH             = 0.76d0
  real(kind=8),parameter :: mpc            = 3.08d24                 ![cm] from Mpc to cm
  real(kind=8),parameter :: msun           = 1.989d33                ![g] solar mass 
  real(kind=8),parameter :: lambda_LyA_Ang = 1215.67d0               ![A] wavelength of LyA line

end module module_constants

