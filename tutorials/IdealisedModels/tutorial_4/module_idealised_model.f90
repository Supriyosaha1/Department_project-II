module module_idealised_model

  ! template that implements the model of Burchett+2021

  ! This model has a H number density defined by n_H(r) = n_H,0 (r_inner/r)^2  [Eq. 2 of Burchett+21]  
  ! The velocity profile is linear in r (increasing or decreasing) from v_inner at r_inner to v_outer at r_outer:
  !            v(r) = [v_in + (r - r_inner)/(r_outer - r_inner) * (v_outer - v_inner)] * u_r
  !  (with u_r the normalised radial vector). 
  ! 
  ! -> 5 independent parameters which we chose to be : n_H,0, r_inner, r_outer, v_inner, v_outer

  ! We add to this a few cold streams.
  ! The geometry of the cold streams is such that each one is a circular cone, originating at the center of the box,
  !      and limitted by r_inner and r_outer from above. 
  ! The density in the cone is assumed to be a scaled-up version of the outflow (defined with cold_stream_to_outflow_nH_ratio).
  ! The velocity of the cone scales linearly with radius, from v_out_cone_kms to v_in_cone_kms
  
  implicit none

  private

  ! custom parameters read from configuration file
  ! ------------------------------------------------
  ! Properties of the volume-filling outflow 
  real(kind=8) :: nH_norm_cgs = 1.0d0       ! normalisation of the H density profile, in cgs
  real(kind=8) :: r_inner_kpc = 1.0d0       ! inner radius of the CGM, in kpc
  real(kind=8) :: r_outer_kpc = 30.0d0      ! outer radius of the CGM, in kpc
  real(kind=8) :: v_inner_kms = 50.0d0      ! velocity at r_inner_kpc, in km/s
  real(kind=8) :: v_outer_kms = 500.0d0     ! velocity at r_outer_kpc, in km/s
  real(kind=8) :: Temperature = 1.0d4       ! gas temperature [K]
  real(kind=8) :: TurbulentVel_kms = 20.0d0 ! velocity dispersion (sigma) due to turbulent velocity [km/s]
  ! ------------------------------------------------
  ! Properties of the volume-filling outflow 
  real(kind=8)    :: cold_stream_covering_fraction = 0.1  ! each cold stream will have a solid angle of 4 pi * this. 
  real(kind=8),allocatable :: cold_stream_direction(:,:)  ! direction of each cold stream
  integer(kind=4) :: nb_of_cold_streams = 3 ! nb of cold streams
  real(kind=8)    :: cold_stream_to_outflow_nH_ratio = 10.0d0 
  real(kind=8)    :: v_out_cone_kms   ! velocity of the gas in the stream at r_outer
  real(kind=8)    :: v_in_cone_kms    ! velocity of the gas in the stream at r_inner

  
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

  ! useful parameters for cold streams
  real(kind=8) :: cold_stream_solid_angle  ! the solid angle of each stream (defined using 4pi * cold_stream_covering_fraction)
  real(kind=8) :: cold_stream_min_dotprod  ! any (normalised) vector r within the stream has x.kstream = costheta > cold_stream_min_dotprod
  real(kind=8)    :: v_out_cone_cgs, v_in_cone_cgs    ! velocity of the gas in the stream at r_inner

  
  ! public functions:
  public :: idealised_model_get_velocity,idealised_model_get_scatterer_density,idealised_model_get_turbulent_velocity
  public :: idealised_model_get_temperature,idealised_model_get_dust_density
  public :: idealised_model_read_params,idealised_model_print_params

contains


  ! --------------------------------------------------------------------------------  
  ! private functions for convenience ... 
  ! --------------------------------------------------------------------------------
  function is_in_a_stream(x)
    implicit none 
    logical                 :: is_in_a_stream
    real(kind=8),intent(in) :: x(3)
    real(kind=8)            :: r(3)
    real(kind=8)            :: dotprod
    integer(kind=4)         :: istream
    is_in_a_stream = .False.
    r = x - 0.5d0  ! coords relative to center of model.
    r = r / sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))  ! normalised r vector. 
    do istream = 1,nb_of_cold_streams
       dotprod = r(1)*cold_stream_direction(1,istream) + r(2)*cold_stream_direction(2,istream) + r(3)*cold_stream_direction(3,istream)
       if (dotprod >= cold_stream_min_dotprod) then
          is_in_a_stream = .True.
          return
       end if
    end do
    return
  end function is_in_a_stream
  
  ! --------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------
  
  
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
       ! Check if the point is in a stream or in the outflow
       if (is_in_a_stream(x)) then 
          idealised_model_get_velocity(:) = r(:)/r_norm * (v_in_cone_cgs + (r_norm - r_inner_box_units)/(r_outer_box_units - r_inner_box_units)*(v_out_cone_cgs - v_in_cone_cgs))
       else
          idealised_model_get_velocity(:) = r(:)/r_norm * (v_inner_cgs + (r_norm - r_inner_box_units)/(r_outer_box_units - r_inner_box_units)*(v_outer_cgs - v_inner_cgs))
       end if
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
       if (is_in_a_stream(x)) then 
          idealised_model_get_scatterer_density = nH_norm_cgs * r_inner_box_units*r_inner_box_units / d * cold_stream_to_outflow_nH_ratio ! [cm^-3]
       else
          idealised_model_get_scatterer_density = nH_norm_cgs * r_inner_box_units*r_inner_box_units / d ! [cm^-3]
       end if
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
    integer(kind=4)         :: err,i,istream
    logical                 :: section_present
    real(kind=8) :: norm
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
          case ('cold_stream_covering_fraction')
             read(value,*) cold_stream_covering_fraction
          case ('cold_stream_to_outflow_nH_ratio')
             read(value,*) cold_stream_to_outflow_nH_ratio
          case ('v_out_cone_kms')
             read(value,*) v_out_cone_kms
          case ('v_in_cone_kms')
             read(value,*) v_in_cone_kms
          case ('nb_of_cold_streams')
             read(value,*) nb_of_cold_streams
          case ('cold_stream_directions')
             if (nb_of_cold_streams > 0) then
                allocate(cold_stream_direction(3,nb_of_cold_streams))
                read(value,*) (cold_stream_direction(1:3,istream),istream = 1, nb_of_cold_streams)
             end if
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

    ! define derived parameters for cold streams
    cold_stream_solid_angle = 12.56637 * cold_stream_covering_fraction   ! 12.566370 is 4*pi
    cold_stream_min_dotprod = 1.0d0 - cold_stream_solid_angle / 6.28318  ! 6.28318 is 2*pi
    
    ! make sure directions are normalised vectors
    do istream = 1,nb_of_cold_streams
       norm = cold_stream_direction(1,istream)**2 + cold_stream_direction(2,istream)**2 + cold_stream_direction(3,istream)**2 
       cold_stream_direction(:,istream) = cold_stream_direction(:,istream)/sqrt(norm)
    end do
    v_in_cone_cgs      = v_in_cone_kms * kms2cms
    v_out_cone_cgs     = v_out_cone_kms * kms2cms

    
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
