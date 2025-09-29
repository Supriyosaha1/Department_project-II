module module_idealised_model

  ! template implements a uniform and static sphere with 1 scatterer and no dust 
  
  !use module_constants, only : kb
  
  implicit none

  private

  ! custom parameters read from configuration file 
  real(kind=8) :: ColumnDensity_cgs = 1d15           ! column density from center to edge of sphere [cm^-2]
  real(kind=8) :: idealised_model_box_size_cm = 1e24 ! size of the box in cm
  real(kind=8) :: Radius_boxUnits = 0.48             ! radius of sphere in box-size units
  real(kind=8) :: Temperature = 1d4                  ! gas temperature [K]
  real(kind=8) :: TurbulentVelocity_kms = 10.0d0     ! velocity dispersion (sigma) due to turbulent velocity [km/s]
  
  ! useful derived parameters
  real(kind=8) :: Radius_cgs ! Radius of sphere in cm
  real(kind=8) :: Radius_squared_box_units ! square of radius of sphere, in units of box size.
  
  ! public functions:
  public :: idealised_model_get_velocity,idealised_model_get_scatterer_density,idealised_model_get_turbulent_velocity
  public :: idealised_model_get_temperature,idealised_model_get_dust_density
  public :: idealised_model_read_params,idealised_model_print_params
  ! public parameters:
  public :: idealised_model_box_size_cm
  
contains
  
  function idealised_model_get_velocity(x)
    ! return bulk velocity of gas at position x
    ! -- x normalised to box size (i.e. values between 0 and 1)
    ! -- idealised_model_get_velocity is in cm/s 
    implicit none
    real(kind=8),dimension(3),intent(in)  :: x
    real(kind=8),dimension(3) :: idealised_model_get_velocity
    idealised_model_get_velocity(:) = 0.0d0   ! [cm/s]
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
    if (d < Radius_squared_box_units) then 
       idealised_model_get_scatterer_density = ColumnDensity_cgs / Radius_cgs ! [cm^-3]
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
    idealised_model_get_turbulent_velocity =  TurbulentVelocity_kms * 1d5 ! [cm/s]
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
          case ('ColumnDensity_cgs')
             read(value,*) ColumnDensity_cgs
          case ('boxsize_cm')
             read(value,*) idealised_model_box_size_cm
          case ('Radius_boxUnits')
             read(value,*) Radius_boxUnits
          case ('Temperature')
             read(value,*) Temperature
          case ('TurbulentVelocity_kms')
             read(value,*) TurbulentVelocity_kms
          end select
       end do
    end if
    close(10)
    
    ! Convert radius to cgs  
    Radius_cgs = Radius_boxUnits * idealised_model_box_size_cm
    ! compute useful parameters 
    Radius_squared_box_units = Radius_boxUnits**2
    
    return
  end subroutine idealised_model_read_params

  
  subroutine idealised_model_print_params(unit)
    implicit none
    integer(kind=4),optional,intent(in)  :: unit
    if (present(unit)) then
       write(unit,'(a,a,a)')     '[IdealisedModel]'
       write(unit,'(a,ES13.6)')  '  ColumnDensity_cgs     = ',ColumnDensity_cgs
       write(unit,'(a,ES13.6)')  '  boxsize_cm            = ',idealised_model_box_size_cm
       write(unit,'(a,ES13.6)')  '  Radius_boxUnits       = ',Radius_boxUnits
       write(unit,'(a,ES13.6)')  '  Temperature           = ',Temperature
       write(unit,'(a,ES13.6)')  '  TurbulentVelocity_kms = ',TurbulentVelocity_kms
    else
       write(*,'(a,a,a)')        '[IdealisedModel]'
       write(*,'(a,ES13.6)')     '  ColumnDensity_cgs     = ',ColumnDensity_cgs
       write(*,'(a,ES13.6)')     '  boxsize_cm            = ',idealised_model_box_size_cm
       write(*,'(a,ES13.6)')     '  Radius_boxUnits       = ',Radius_boxUnits
       write(*,'(a,ES13.6)')     '  Temperature           = ',Temperature
       write(*,'(a,ES13.6)')     '  TurbulentVelocity_kms = ',TurbulentVelocity_kms
    end if
  end subroutine idealised_model_print_params

  
end module module_idealised_model
