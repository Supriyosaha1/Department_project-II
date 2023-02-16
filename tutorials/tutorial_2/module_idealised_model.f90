module module_idealised_model

  ! template that implements the model of Burchett+2021

  ! This model has a H number density defined by n_H(r) = n_H,0 (r_inner/r)^2  [Eq. 2 of Burchett+21]  
  ! The velocity profile is linear in r (increasing or decreasing) from v_inner at r_inner to v_outer at r_outer:
  !            v(r) = v_in + (r - r_inner)/(r_outer - r_inner) * (v_outer - v_inner)
  ! 
  ! -> 5 independent parameters which we chose to be : n_H,0, r_inner, r_outer, v_inner, v_outer
  
  implicit none

  private

  ! custom parameters read from configuration file 
  real(kind=8) :: nH_norm_cgs = 1.0d0       ! normalisation of the H density profile, in cgs
  real(kind=8) :: r_inner_kpc = 1.0d0       ! inner radius of the CGM, in kpc
  real(kind=8) :: r_outer_kpc = 30.0d0      ! outer radius of the CGM, in kpc
  real(kind=8) :: v_inner_kms = 50.0d0      ! velocity at r_inner_kpc, in km/s
  real(kind=8) :: v_outer_kms = 500.0d0     ! velocity at r_outer_kpc, in km/s
  real(kind=8) :: Temperature = 1.0d4       ! gas temperature [K]
  real(kind=8) :: TurbulentVel_kms = 20.0d0 ! velocity dispersion (sigma) due to turbulent velocity [km/s]

  ! useful constants for conversions
  real(kind=8),parameter :: kpc2cm = 3.0856775814913673d21
  real(kind=8),parameter :: kms2cms = 1d5

  ! mandatory parameters 
  real(kind=8),public :: idealised_model_box_size_cm  ! size of computational box [cm]

  ! useful derived parameters
  real(kind=8) :: inner_radius_squared_box_units ! square of inner radius of sphere, in units of box size.
  real(kind=8) :: outer_radius_squared_box_units ! square of outer radius of sphere, in units of box size.
  real(kind=8) :: r_inner_cgs,r_outer_cgs,v_inner_cgs,v_outer_cgs,TurbulentVel_cgs ! user-defined parameters in cgs units.
  real(kind=8) :: r_inner_box_units,r_outer_box_units ! user-defined parameters in box units
  
  ! public functions:
  public :: idealised_model_get_velocity,idealised_model_get_scatterer_density,idealised_model_get_turbulent_velocity
  public :: idealised_model_get_temperature,idealised_model_get_dust_density
  public :: idealised_model_read_params,idealised_model_print_params
  
contains
  
  function idealised_model_get_velocity(x)
    ! return bulk velocity of gas at position x
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_velocity is in cm/s 
    implicit none
    real(kind=8),dimension(3),intent(in)  :: x
    real(kind=8),dimension(3) :: idealised_model_get_velocity
    real(kind=8) :: d, v, r(3), r_norm
    r(:) = x(:) - 0.5d0 
    d = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
    r_norm = sqrt(d)
    if (d <= outer_radius_squared_box_units .and. d >= inner_radius_squared_box_units) then
       ! we are within the CGM -> use linear velocity profile
       idealised_model_get_velocity(:) = r(:)/r_norm * (v_inner_cgs + (r_norm - r_inner_box_units)/(r_outer_box_units - r_inner_box_units)*(v_outer_cgs - v_inner_cgs))
    else
       idealised_model_get_velocity(:) = 0.0d0   ! [cm/s]
    end if
    return
  end function idealised_model_get_velocity


  function idealised_model_get_scatterer_density(x)
    ! return number density of first scatterer in list at position x
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_scatterer_density is in [cm^-3]
    implicit none
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8)  :: idealised_model_get_scatterer_density
    real(kind=8) :: d

    d = (x(1)-0.5d0)*(x(1)-0.5d0) + (x(2)-0.5d0)*(x(2)-0.5d0) + (x(3)-0.5d0)*(x(3)-0.5d0)
    if (d <= outer_radius_squared_box_units .and. d >= inner_radius_squared_box_units) then 
       idealised_model_get_scatterer_density = nH_norm_cgs * r_inner_box_units*r_inner_box_units / d ! [cm^-3]
    else
       idealised_model_get_scatterer_density = 0.0d0 ! [cm^-3]
    end if
    return
  end function idealised_model_get_scatterer_density


  function idealised_model_get_turbulent_velocity(x)
    ! return the turbulent velocity dispersion at position x
    ! NB: does not need to go to zero where beyond the sphere (where density is zero). 
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_turbulent_velocity is in [cm/s]
    implicit none
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8) :: idealised_model_get_turbulent_velocity
    idealised_model_get_turbulent_velocity =  TurbulentVel_cgs ! [cm/s]
    return
  end function idealised_model_get_turbulent_velocity


  function idealised_model_get_temperature(x)
    ! return the temperature at position x
    ! NB: does not need to go to zero where beyond the sphere (where density is zero). 
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_temperature is in [K]
    implicit none
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8) :: idealised_model_get_temperature
    idealised_model_get_temperature = Temperature
    return
  end function idealised_model_get_temperature


  function idealised_model_get_dust_density(x)
    ! NB: we should use the module_dust_model here to invert a Tau into a density (at a given wavelength) ? 
    ! return number density of dust at position x
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_dust_density is in [cm^-3]
    implicit none
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8) :: idealised_model_get_dust_density
    idealised_model_get_dust_density = 0.0d0 ! [cm^-3]
    return
  end function idealised_model_get_dust_density


  subroutine idealised_model_read_params(pfile)   
    implicit none
    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[IdealisedModel]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('nH_norm','nH_norm_cgs')
             read(value,*) nH_norm_cgs
          case ('r_inner_kpc')
             read(value,*) r_inner_kpc
          case ('r_outer_kpc')
             read(value,*) r_outer_kpc
          case ('v_inner_kms')
             read(value,*) v_inner_kms
          case ('v_outer_kms')
             read(value,*) v_outer_kms
          case ('Temperature')
             read(value,*) Temperature
          case ('TurbulentVel_kms')
             read(value,*) TurbulentVel_kms
          end select
       end do
    end if
    close(10)

    ! convert parameters to cgs
    r_inner_cgs      = r_inner_kpc * kpc2cm
    r_outer_cgs      = r_outer_kpc * kpc2cm
    v_inner_cgs      = v_inner_kms * kms2cms
    v_outer_cgs      = v_outer_kms * kms2cms
    TurbulentVel_cgs = TurbulentVel_kms * kms2cms

    
    ! make sure to define mandatory parameter idealised_model_box_size_cm 
    idealised_model_box_size_cm = 3d0 * r_outer_cgs ! JB - this is not ideal ... 

    ! convert distances to box units too
    r_inner_box_units = r_inner_cgs / idealised_model_box_size_cm
    r_outer_box_units = r_outer_cgs / idealised_model_box_size_cm
    ! compute useful parameters 
    inner_radius_squared_box_units = r_inner_box_units**2
    outer_radius_squared_box_units = r_outer_box_units**2
    
    return
  end subroutine idealised_model_read_params

  
  subroutine idealised_model_print_params(unit)
    implicit none
    integer(kind=4),optional,intent(in)  :: unit
    if (present(unit)) then
       write(unit,'(a,a,a)')     '[IdealisedModel]'
       write(unit,'(a,ES13.6)')  '  nH_norm_cgs           = ',nH_norm_cgs
       write(unit,'(a,ES13.6)')  '  r_inner_kpc           = ',r_inner_kpc
       write(unit,'(a,ES13.6)')  '  r_outer_kpc           = ',r_outer_kpc
       write(unit,'(a,ES13.6)')  '  v_inner_kms           = ',v_inner_kms
       write(unit,'(a,ES13.6)')  '  v_outer_kms           = ',v_outer_kms
       write(unit,'(a,ES13.6)')  '  Temperature           = ',Temperature
       write(unit,'(a,ES13.6)')  '  TurbulentVel_kms      = ',TurbulentVel_kms
       write(unit,'(a,ES13.6)')  '  box_size_cm           = ',idealised_model_box_size_cm
    else
       write(*,'(a,a,a)')        '[IdealisedModel]'
       write(*,'(a,ES13.6)')  '  nH_norm_cgs           = ',nH_norm_cgs
       write(*,'(a,ES13.6)')  '  r_inner_kpc           = ',r_inner_kpc
       write(*,'(a,ES13.6)')  '  r_outer_kpc           = ',r_outer_kpc
       write(*,'(a,ES13.6)')  '  v_inner_kms           = ',v_inner_kms
       write(*,'(a,ES13.6)')  '  v_outer_kms           = ',v_outer_kms
       write(*,'(a,ES13.6)')  '  Temperature           = ',Temperature
       write(*,'(a,ES13.6)')  '  TurbulentVel_kms      = ',TurbulentVel_kms
       write(*,'(a,ES13.6)')  '  box_size_cm           = ',idealised_model_box_size_cm
    end if
  end subroutine idealised_model_print_params

  
end module module_idealised_model
