module module_scatterer_model

  use module_constants
  use module_utils, only : isotropic_direction, anisotropic_direction_HIcore, anisotropic_direction_Rayleigh
  !--PEEL--
  use module_utils, only : anisotropic_probability_HIcore,anisotropic_probability_Rayleigh
  !--LEEP--
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Defining a type corresponding to a spectral transition
  
  type, public :: scatterer
     ! definition of atomic values
     real(kind=8)                  :: m_ion             ! Atomic mass of the ion
     character(10)                 :: name_ion          ! Name of the ion, corresponding to the nomenclature of the Krome code ('Ions')
     real(kind=8)                  :: lambda_cm         ! Wavelength of the transition
     real(kind=8)                  :: f                 ! Oscillator strength, low - up (absorption from the lower level to the upper level)
     real(kind=8)                  :: A                 ! Einstein coefficient of the transition, low - up
     integer(kind=4)               :: n_fluo            ! Number of fluorescent channels
     real(kind=8),allocatable      :: lambda_fluo_cm(:) ! Wavelengths of the fluorescent channels
     real(kind=8),allocatable      :: A_fluo(:)         ! Einstein coefficients of the fluorescent channels, up - low
     
     ! useful pre-computed quantities
     real(kind=8)                  :: A_over_fourpi     ! Einstein coefficient of the low-up transition divided by four * pi
     real(kind=8)                  :: A_tot             ! Sum of the (one) low-up and all the up - low Einstein coefficients
     real(kind=8)                  :: nu                ! Frequency of the transition
     real(kind=8),allocatable      :: nu_fluo(:)        ! Frequencies of the fluorescent channels
     real(kind=8)                  :: sigma_factor      ! Cross-section factor-> multiply by Voigt(x,a)/delta_nu_doppler to get sigma.

     logical                       :: recoil       = .false.     ! if set to true, recoil effect is computed
     logical                       :: isotropic    = .true.      ! if set to true, scattering events will be isotropic
     logical                       :: core_skip    = .false.     ! if true, skip scatterings in the core of the line (as in Smith+15).
     real(kind=8)                  :: xcritmax     = -1d10       ! core-skipping will truncate at min(xcrit, xcritmax) -> set to a low (but positive) value to activate. 

  end type scatterer
  

  public :: get_tau, scatter, read_scatterer_params, print_scatterer_params 
  public :: peeloff_weight
  
  
contains

  
  function get_tau(nscat, vth_sq_times_m, vturb, distance_to_border_cm, nu_cell, s)

    ! --------------------------------------------------------------------------
    ! compute optical depth of Hydrogen over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nscat                 : number density of scatterers                    [ cm^-3 ]
    ! - vth_sq_times_m        : 2 * kb * T / amu 
    ! - distance_to_border_cm : distance over which we compute tau              [ cm ]
    ! - nu_cell               : photon's frequency in the frame of the cell     [ Hz ]
    ! OUTPUT :
    ! - get_tau               : optical depth of scatterer over distance_to_border_cm
    ! --------------------------------------------------------------------------
    type(scatterer),intent(in) :: s
    real(kind=8),intent(in)    :: nscat,vth_sq_times_m,vturb,distance_to_border_cm,nu_cell
    real(kind=8)               :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau

    ! compute Doppler width and a-parameter

    delta_nu_doppler = sqrt(vth_sq_times_m / s%m_ion + vturb**2*1d10) / s%lambda_cm
    a = s%A_over_fourpi / delta_nu_doppler
 
    ! Cross section of scatterer
    x_cell = (nu_cell - s%nu)/delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma = s%sigma_factor / delta_nu_doppler * h_cell

    get_tau = sigma * nscat * distance_to_border_cm

    return

  end function get_tau


  
  subroutine scatter(vcell,vth_sq_times_m,vturb,nu_cell,k,nu_ext,iran,xcrit,s)
    ! ---------------------------------------------------------------------------------
    ! perform scattering event on an ion
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth_sq_times_m : thermal velocity square times atomic mass
    ! - vturb    : turbulent velocity                                  [km/s]
    ! - nu_cell  : frequency of incoming photon in cell's rest-frame   [ Hz ] 
    ! - k        : propagaction vector (normalized) 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - iran     : random number generator seed
    ! - s        : Variable de type scatterer, 
    ! OUTPUTS :
    ! - nu_cell  : updated frequency in cell's frame   [ Hz ]
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! - k        : updated propagation direction
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------
    !
    ! Notes on the phase function, for HI :
    ! -----------------------------
    ! - for core photons (|x| < 0.2) we use P(mu) = 11/24 + 3/24 * mu**2
    ! - for wing photons (|x| > 0.2) we use P(mu) = 3/8 * (1 + mu**2) [this is Rayleigh]
    ! where mu = cos(theta), (and theta in [0,pi]).
    ! ---------------------------------------------------------------------------------

    type(scatterer),intent(in)              :: s
    real(kind=8),intent(inout)              :: nu_cell, nu_ext
    real(kind=8),dimension(3),intent(inout) :: k
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth_sq_times_m, vturb
    !--CORESKIP--
    real(kind=8),intent(in)                 :: xcrit
    real(kind=8)                            :: xc
    !--PIKSEROC--
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba
    real(kind=8)                            :: x_atom, dopwidth
    real(kind=8),dimension(3)               :: knew
    integer(kind=4)                         :: i

    ! !--CORESKIP--  sanity check ... 
    ! if (.not. s%core_skip .and. xcrit .ne. 0.0d0) then
    !    print*,'ERROR: core skipping is off but xcrit is not zero ... '
    !    stop
    ! end if
    ! if (s%core_skip)  then
    !    xc = min(xcrit,s%xcritmax)
    ! else
    !    xc=0.0d0
    ! endif
    ! !--PIKSEROC--

    dopwidth = sqrt(vth_sq_times_m / s%m_ion + vturb**2*1d10)
    
    ! define x_cell & a
    delta_nu_doppler = dopwidth / s%lambda_cm 
    a = s%A_over_fourpi / delta_nu_doppler
    x_cell = (nu_cell - s%nu) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * dopwidth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    !--CORESKIP--
    uper   = sqrt(xc**2-log(ruper))*cos(twopi*r2)
    !--PIKSEROC--
    uper   = uper * dopwidth  ! from x to velocity
    

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    ! a) There are no fluorescent channels, the decay is resonant
    if(s%n_fluo == 0) then
       nu_atom = nu_cell - nu_ext * upar/clight
    ! b) There are fluorescent channels, a random number decides the decay channel
    else
       r2 = ran3(iran)
       ! i) The decay is resonant
       if(r2 < s%A / s%A_tot) then
          nu_atom = nu_cell - nu_ext * upar/clight
       ! ii) The decay is through a fluorescent channel
       else
          proba = s%A / s%A_tot
          do i=1,s%n_fluo
             if(r2 < s%A_fluo(i) / s%A_tot + proba) then
                ! For fluroescent channels, the 'scattering' event is non-coherent,
                ! we lose the information of frequency of the incoming photon.
                nu_atom = s%nu_fluo(i)
                exit
             end if
             proba = proba +  s%A_fluo(i) / s%A_tot
          end do
       end if
    end if


    ! 4/ determine direction of scattered photon
    if (s%isotropic) then
       call isotropic_direction(knew,iran)
       mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
       bu = sqrt(1.0d0 - mu*mu)
    else ! Works only for HI
       x_atom  = (nu_atom - s%nu) / delta_nu_doppler
       if (abs(x_atom) < 0.2d0) then ! core scattering 
          call anisotropic_direction_HIcore(k,knew,mu,bu,iran)
       else ! wing scattering 
          call anisotropic_direction_Rayleigh(k,knew,mu,bu,iran)
       end if
    end if

    ! 5/ recoil effect 
    if (s%recoil) then ! Works only for HI
       nu_atom = nu_atom / (1.0d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.0d0-mu))
    end if
    
    ! 6/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.0d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter



  !--PEEL--
  function peeloff_weight(vcell,vth_sq_times_m,vturb,nu_ext,kin,kout,iran,s)
    
    ! ---------------------------------------------------------------------------------
    ! Compute probability that a photon coming along kin scatters off in direction kout.
    ! Also update nu_ext to external-frame frequency along kout
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth_sq_times_m : thermal velocity square times atomic mass
    ! - vturb    : turbulent velocity                                  [km/s]
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - kin      : propagation vector (normalized)
    ! - kout     : direction after interaction (fixed)
    ! - iran     : random number generator seed
    ! OUTPUTS :
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------
    !
    ! Notes on the phase function, for HI :
    ! -----------------------------
    ! - for core photons (|x| < 0.2) we use P(mu) = 11/24 + 3/24 * mu**2
    ! - for wing photons (|x| > 0.2) we use P(mu) = 3/8 * (1 + mu**2) [this is Rayleigh]
    ! where mu = cos(theta), (and theta in [0,pi]).
    ! ---------------------------------------------------------------------------------

    type(scatterer),intent(in)              :: s
    real(kind=8),intent(inout)              :: nu_ext
    real(kind=8),dimension(3),intent(in)    :: kin, kout
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth_sq_times_m, vturb
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba
    real(kind=8)                            :: x_atom,nu_cell,dopwidth
    integer(kind=4)                         :: i

    ! compute frequency in cell's frame 
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext

    dopwidth = sqrt(vth_sq_times_m / s%m_ion + vturb**2*1d10)

    ! define x_cell & a
    delta_nu_doppler = dopwidth / s%lambda_cm 
    a = s%A_over_fourpi / delta_nu_doppler
    x_cell = (nu_cell - s%nu) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * dopwidth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * dopwidth  ! from x to velocity

    
    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    ! a) There are no fluorescent channels, the decay is resonant
    if(s%n_fluo == 0) then
       nu_atom = nu_cell - nu_ext * upar/clight
    ! b) There are fluorescent channels, a random number decides the decay channel
    else
       r2 = ran3(iran)
       ! i) The decay is resonant
       if(r2 < s%A / s%A_tot) then
          nu_atom = nu_cell - nu_ext * upar/clight
       ! ii) The decay is through a fluorescent channel
       else
          proba = s%A / s%A_tot
          do i=1,s%n_fluo
             if(r2 < s%A_fluo(i) / s%A_tot + proba) then
                ! For fluroescent channels, the 'scattering' event is non-coherent,
                ! we lose the information of frequency of the incoming photon.
                nu_atom = s%nu_fluo(i)
                exit
             end if
             proba = proba +  s%A_fluo(i) / s%A_tot
          end do
       end if
    end if

    ! 4/ determine direction of scattered photon
    if (s%isotropic) then
       peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
       mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
       bu = sqrt(1.0d0 - mu*mu)
    else ! Works only for HI
       x_atom  = (nu_atom - s%nu) / delta_nu_doppler
       if (abs(x_atom) < 0.2) then ! core scattering 
          peeloff_weight = anisotropic_probability_HIcore(kin,kout,mu,bu)
       else ! wing scattering 
          peeloff_weight = anisotropic_probability_Rayleigh(kin,kout,mu,bu)
       end if
    end if

    ! 5/ recoil effect 
    if (s%recoil) then  ! Works only for HI
       nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
    end if
    
    ! 6/ compute freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    
  end function peeloff_weight
!--LEEP--



  subroutine read_scatterer_params(sfile,pfile,s,index_call_read)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters scatterer in the atomic parameter file pfile
    ! ---------------------------------------------------------------------------------

    character(*),intent(in)       :: sfile,pfile
    type(scatterer),intent(inout) :: s
    integer(kind=4),intent(in)    :: index_call_read ! If 1, call read voigt and uparallel params
    character(1000)               :: line,name,value
    integer(kind=4)               :: err,i
    logical                       :: file_exists

    INQUIRE(FILE=sfile, EXIST=file_exists)
    if(.not. file_exists) then
       print*, 'File '//trim(sfile)//' does not exist, stopping the program.'
       stop
    else
       open(unit=10,file=trim(sfile),status='old',form='formatted')
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('m_ion')
             read(value,*) s%m_ion
          case ('name_ion')
             read(value,*) s%name_ion
          case ('lambda_cm')
             read(value,*) s%lambda_cm
             s%nu = clight / s%lambda_cm
          case ('A')
             read(value,*) s%A
             s%A_over_fourpi = s%A / fourpi
          case ('f')
             read(value,*) s%f
             s%sigma_factor = sqrtpi*e_ch**2*s%f/me/clight
          case ('n_fluo')
             read(value,*) s%n_fluo
             allocate(s%lambda_fluo_cm(s%n_fluo), s%A_fluo(s%n_fluo), s%nu_fluo(s%n_fluo))
          case ('lambda_fluo_cm')
             read(value,*) s%lambda_fluo_cm(:)
             s%nu_fluo(:) = clight / s%lambda_fluo_cm(:)
          case ('A_fluo')
             read(value,*) s%A_fluo(:)
             s%A_tot = s%A + sum(s%A_fluo(:))

             !    !--CORESKIP--
             ! case ('core_skip') 
             !    read(value,*) core_skip
             ! case ('xcritmax')
             !    read(value,*) xcritmax
             !    !--PIKSEROC--

          end select
       end do
       close(10)
    end if

    ! !--CORESKIP--  sanity check ... 
    ! if (core_skip .and. xcritmax <= 0.0d0) then
    !    print*,'ERROR: core skipping is on but xcritmax is not set... '
    !    stop
    ! end if
    ! !--PIKSEROC-- 
    
    if(index_call_read == 1) then
       call read_uparallel_params(pfile)
       call read_voigt_params(pfile)
    end if

    return

  end subroutine read_scatterer_params


  
  subroutine print_scatterer_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    write(unit,*) 'print_scatterer_params: to redo'

    ! if (present(unit)) then 
    !    write(unit,'(a,a,a)') '[HI]'
    !    write(unit,'(a,L1)') '  recoil    = ',recoil
    !    write(unit,'(a,L1)') '  isotropic = ',isotropic
    !    !--CORESKIP--
    !    write(unit,'(a,L1)')     '  core_skip = ',core_skip
    !    write(unit,'(a,ES10.3)') '  xcritmax     = ',xcritmax
    !    !--PIKSEROC--
    !    write(unit,'(a)') ''
    !    call print_uparallel_params(unit)
    !    call print_voigt_params(unit)
    ! else
    !    write(*,'(a,a,a)') '[HI]'
    !    write(*,'(a,L1)') '  recoil    = ',recoil
    !    write(*,'(a,L1)') '  isotropic = ',isotropic
    !    !--CORESKIP--
    !    write(*,'(a,L1)')     '  core_skip = ',core_skip
    !    write(*,'(a,ES10.3)') '  xcritmax     = ',xcritmax
    !    !--PIKSEROC--
    !    write(*,'(a)') ''
    !    call print_uparallel_params()
    !    call print_voigt_params()
    ! end if
    
    return
    
  end subroutine print_scatterer_params


end module module_scatterer_model
