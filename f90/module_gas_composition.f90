module module_gas_composition

  ! arbitrary mix of scatterers and dust ... 

  use module_scatterer_model
  use module_dust_model
  use module_random
  use module_constants, only: kb, amu
  use module_idealised_model

  implicit none

  private

  type, public :: gas
     ! fluid
     real(kind=8)                :: v(3)                   ! gas velocity [cm/s]
     real(kind=8),allocatable    :: number_density(:)      ! number density of each scatterer [cm^-3]
     real(kind=8)                :: vth_sq_times_m         ! 2 kb T / amu  , dopwidth = sqrt(vth_sq_times_m / m_ion + vturb**2)
     real(kind=8)                :: vturb                  ! Turbulent velocity in every cell [cm/s]
     ! DUST -> model of Laursen, Sommer-Larsen and Andersen 2009.
     ! ->  ndust = (nHI + f_ion nHII)*Z/Zref
     ! f_ion and Zref are two free parameters . 
     real(kind=8)                :: ndust                  ! pseudo-numerical density of dust particles [#/cm3]
  end type gas
  real(kind=8),public         :: box_size_cm          ! size of simulation box in cm.
  
  type(scatterer),allocatable :: scatterer_list(:)    ! List of scatterers, defined from scatterer_list (below)
  logical                     :: DI_present
  integer(kind=4)             :: metal_number         ! Number of metallic ions, equal to nscatterer if no HI and scatterer-1 if there is HI
  integer(kind=4)             :: element_number
  integer(kind=4),allocatable :: scatterer_ion_index(:)
  character(10),allocatable   :: name_ions_no_repetition(:)

  integer(kind=4) :: HI_1216_flag = -1  ! will be set to the flag value for interaction through HI-Lya (i.e. the value of flag from get_scatter_flag) if HI-Lya is present. Stays -1 otherwise.
  integer(kind=4) :: HI_1216_index = -1 ! the index of HI-1216 in the scatterer list. 
  
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! mixture parameters
  integer(kind=4)             :: nscatterer                              ! Number of scatterers in the run
  character(100),allocatable  :: scatterer_names(:)                      ! List of names of scatterer (i.e. lines)  ! Follow nomenclature of the tuto
  character(1000)             :: atomic_data_dir = '../ions_parameters/' ! directory where the atomic data files are (usually in .../rascas/ions_parameters/)
  character(2000)             :: krome_data_dir  = './'                  ! directory where the density of metallic ions are, from the Krome runs
  real(kind=8)                :: f_ion           = 0.01                  ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)                :: Zref            = 0.005                 ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC.
  real(kind=8)                :: vturb_kms       = 2d1                   ! Constant turbulent velocity accross the simulation,  in km/s
  logical                     :: ignoreDust      = .false.               ! if true, no dust in the mix
                                                 ! (NB: in this case, dust is ignored during RT, but dust properties are still defined)
    
  ! ! possibility to overwrite ramses values with an ad-hoc model 
  ! logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  ! real(kind=8),allocatable :: fix_nscat(:)                  ! ad-hoc HI density (H/cm3)
  ! real(kind=8)             :: fix_vth_sq_times_m            ! ad-hoc thermal velocity (cm/s)
  ! real(kind=8)             :: fix_ndust           = 0.0d0   ! ad-hoc dust number density (/cm3)
  ! real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static...
  ! real(kind=8)             :: fix_vturb           = 0.0d0   ! ad-hoc turbulent velocity
  ! real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  ! --------------------------------------------------------------------------


  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  public :: get_metal_ion_names

  public :: gas_from_idealised_model, gas_get_dopvel

  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau
  !--LEEP--
  public :: nscatterer, metal_number, element_number, krome_data_dir
  
contains

  
  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var,g) 

    ! define gas contents from ramses raw data

    use module_ramses
    
    character(2000),intent(in)                                  :: repository 
    integer(kind=4),intent(in)                                  :: snapnum
    integer(kind=4),intent(in)                                  :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar+metal_number,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out)              :: g
    integer(kind=4)                                             :: ileaf, i, j
    real(kind=8),dimension(:),allocatable                       :: T, nhi, metallicity, nhii, ndi
    real(kind=8),dimension(:,:),allocatable                     :: v

    ! allocate gas-element array
    allocate(g(nleaf))

    ! Allocate the density variable
    do ileaf=1,nleaf
       allocate(g(ileaf)%number_density(element_number))
    end do


    box_size_cm = ramses_get_box_size_cm(repository,snapnum)

    ! compute velocities in cm / s
    write(*,*) '-- module_gas_composition : extracting velocities from ramses '
    allocate(v(3,nleaf))
    call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar+metal_number,ramses_var,v)

    do ileaf = 1,nleaf
       g(ileaf)%v = v(:,ileaf)
    end do
    deallocate(v)

    ! get nHI and temperature from ramses
    write(*,*) '-- module_gas_composition : extracting nHI and T from ramses '
    allocate(T(nleaf),nhi(nleaf))
    call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar+metal_number,ramses_var,T,nhi)

    ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
    write(*,*) '-- module_gas_composition : extracting ndust from ramses '
    allocate(metallicity(nleaf),nhii(nleaf))
    call ramses_get_metallicity(nleaf,nvar+metal_number,ramses_var,metallicity)

    call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar+metal_number,ramses_var,nhii)
    nhii = nhii - nhi
    do ileaf = 1,nleaf
       g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
    end do
    deallocate(metallicity,nhii)

    ! Get deuterium if needed
    if(DI_present) then
       allocate(ndi(nleaf))
       call ramses_get_deuterium(repository,snapnum,nleaf,nvar+metal_number,ramses_var,ndi)
    end if

    ! For each element, determine if it is hydrogen (or deuterium), or a metallic ion. For metallic ions, the data is stored in ramses_var, from nvar+1 to nvar+metal_number
    j = 1
    do i=1,element_number
       if(name_ions_no_repetition(i) == 'HI') then
          do ileaf = 1,nleaf
             g(ileaf)%number_density(i) = nhi(ileaf)
          end do
       else if(name_ions_no_repetition(i) == 'DI') then
          do ileaf = 1,nleaf
             g(ileaf)%number_density(i) = ndi(ileaf)
          end do
       else
          do ileaf = 1,nleaf
             g(ileaf)%number_density(i) = ramses_var(nvar+j,ileaf)
          end do
          j = j+1
       end if
    end do
    deallocate(nhi)
    if(DI_present) deallocate(ndi)
    
    ! compute thermal velocity
    g(:)%vth_sq_times_m = 2d0*kb*T / amu
    deallocate(T)

    ! Set turbulent velocity
    g(:)%vturb          = vturb_kms * 1d5  ! [cm/s]


    return

  end subroutine gas_from_ramses_leaves



  subroutine gas_from_idealised_model(nleaf, g, x_leaf, leaf_level)
    ! work only for one ion/element
    
    integer(kind=4),intent(in)                       :: nleaf
    type(gas),dimension(:),allocatable,intent(out)   :: g
    real(kind=8),intent(in),dimension(nleaf,3)       :: x_leaf
    integer(kind=4),intent(in),dimension(nleaf)      :: leaf_level ! not useful anymore (was useful for MC sampling of cell properties), kept in for now...
    integer(kind=4)                                  :: ileaf
    real(kind=8)                                     :: T
    real(kind=8),dimension(3)                        :: xtmp
    
    ! allocate gas-element array
    allocate(g(nleaf))

    ! Allocate the density variable
    ! work only for one ion/element
    if(element_number /= 1)then
       print*,'element_number =',element_number
       print*,'idealised model works only for one element/ion, check the config file...'
       stop
    endif
    do ileaf=1,nleaf
       allocate(g(ileaf)%number_density(element_number))
    end do
    
    box_size_cm = idealised_model_box_size_cm

    do ileaf = 1, nleaf
       xtmp(1:3) = x_leaf(ileaf,1:3)
       ! get gas velocity
       g(ileaf)%v(:) = idealised_model_get_velocity(xtmp)
       ! get density of scatterer
       g(ileaf)%number_density(1) = idealised_model_get_scatterer_density(xtmp)
       ! get temperature
       T = idealised_model_get_temperature(xtmp)
       g(ileaf)%vth_sq_times_m = 2 * kb * T / amu
       ! get vturb
       g(ileaf)%vturb = idealised_model_get_turbulent_velocity(xtmp)
       ! get ndust
       g(ileaf)%ndust = idealised_model_get_dust_density(xtmp)
    end do

    return
  end subroutine gas_from_idealised_model


  !function gas_get_min_dopvel(x)
  !  implicit none
  !  real(kind=8),dimension(3),intent(in) :: x
  !  real(kind=8) :: T,vth,vth_sq_times_m,vturb,vmin
  !  real(kind=8) ::gas_get_min_dopvel
  !  integer(kind=4) :: i
  !  T = idealised_model_get_temperature(x)
  !  vth_sq_times_m = 2 * kb * T / amu
  !  vturb = idealised_model_get_turbulent_velocity(x)
  !  vmin = 1d20
  !  do i = 1, nscatterer
  !     vmin = min(vmin,vth_sq_times_m/scatterer_list(i)%m_ion)
  !  end do
  !  gas_get_min_dopvel = sqrt(vmin + vturb*vturb)
  !  return
  !end function gas_get_min_dopvel


  function gas_get_dopvel(x)
    ! only for idealised_model with only one element/ion!
    implicit none
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8)                         :: T,vth_sq_times_m,vturb
    real(kind=8)                         ::gas_get_dopvel
    T = idealised_model_get_temperature(x)
    vth_sq_times_m = 2 * kb * T / amu
    vturb = idealised_model_get_turbulent_velocity(x)
    gas_get_dopvel = sqrt(vth_sq_times_m/scatterer_list(1)%m_ion + vturb*vturb)
    return
  end function gas_get_dopvel
  
  
  ! subroutine overwrite_gas(g)
  !   ! overwrite ramses values with an ad-hoc model

  !   type(gas),dimension(:),intent(inout) :: g
  !   integer(kind=4)                      :: ileaf, nleaf

  !   box_size_cm   = fix_box_size_cm

  !   g(:)%v(1)     = fix_vel
  !   g(:)%v(2)     = fix_vel
  !   g(:)%v(3)     = fix_vel
  !   g(:)%vth_sq_times_m = fix_vth_sq_times_m
  !   g(:)%vturb    = fix_vturb
  !   g(:)%ndust    = fix_ndust

  !   nleaf = size(g)
  !   do ileaf=1,nleaf
  !      g(ileaf)%number_density(:) = fix_nscat(:)
  !   end do
    

  ! end subroutine overwrite_gas


  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity



  !--CORESKIP--
  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran, CS_dist_cm, CS_xcrit)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of ions and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon
    ! - CS_dist_cm: true distance to cell border, not along k. For Core-Skipping only.
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [dust], 2<->nscatterer+1 [scatterer]
    ! - CS_xcrit: critical frequency above which core-skipping applies. For Core-Skipping only.
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer(kind=4),intent(inout)         :: iran
    integer(kind=4)                       :: gas_get_scatter_flag, i
    real(kind=8)                          :: tau_ions(nscatterer), tau_dust, tau_cell, tau, tirage
    !--CORESKIP--
    real(kind=8),intent(in)               :: CS_dist_cm
    real(kind=8),intent(inout)            :: CS_xcrit
    real(kind=8)                          :: CS_tau_cell,delta_nu_doppler,a,xcw,x
    !--PIKSEROC--
    
    ! compute optical depths for different components of the gas.
    tau_cell = 0.0d0 
    do i = 1, nscatterer
       tau_ions(i) = get_tau(cell_gas%number_density(scatterer_ion_index(i)),cell_gas%vth_sq_times_m, cell_gas%vturb, distance_to_border_cm, nu_cell, scatterer_list(i))
       tau_cell = tau_cell + tau_ions(i) 
    end do
    if(ignoreDust)then
       tau_dust = 0.0d0
    else
       tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    endif
    tau_cell = tau_cell + tau_dust
    

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 

       gas_get_scatter_flag = 0   ! no scatter
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    
    else  ! the scattering happens inside the cell -> either on dust or on a line scatterer.

       tirage = ran3(iran)
       if(tirage < tau_dust/tau_cell) then

          ! Absorption by dust
          gas_get_scatter_flag = 1            

       else

          ! Absorption by a line scatterer
          tau = tau_dust
          do i=1,nscatterer
             tau = tau + tau_ions(i)
             if(tirage < tau/tau_cell) then
                gas_get_scatter_flag = i+1
                exit
             end if
          end do

          if (gas_get_scatter_flag == HI_1216_flag) then ! if HI Lya scattering setup core-skipping acceleration 
             CS_xcrit = 0.0d0
             if (scatterer_list(HI_1216_index)%core_skip) then
                delta_nu_doppler = sqrt(cell_gas%vth_sq_times_m + cell_gas%vturb**2) / scatterer_list(HI_1216_index)%lambda_cm 
                a = scatterer_list(HI_1216_index)%A_over_fourpi/delta_nu_doppler
                xcw = 6.9184721d0 + 81.766279d0 / (log10(a)-14.651253d0)  ! Smith+15, Eq. 21
                x = (nu_cell - scatterer_list(HI_1216_index)%nu)/delta_nu_doppler
                if (abs(x) < xcw) then ! apply core-skipping
                   CS_tau_cell = tau_ions(HI_1216_index) * (CS_dist_cm / distance_to_border_cm)
                   if (CS_tau_cell > 1.0d0) CS_xcrit = (CS_tau_cell * a)**(1./3.)/5.  ! Smith+15, Eq. 35
                end if
             end if
          end if
          
       end if

       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
       
    end if

    return
  end function gas_get_scatter_flag



  !--CORESKIP-- 
  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran,xcrit)
  !--PIKSEROC--
    
    integer(kind=4),intent(inout)            :: flag
    type(gas),intent(in)                     :: cell_gas
    real(kind=8),intent(inout)               :: nu_cell, nu_ext
    real(kind=8),dimension(3), intent(inout) :: k
    integer(kind=4),intent(inout)            :: iran
    !--CORESKIP--
    real(kind=8),intent(in)                  :: xcrit
    !--PIKSEROC--
    integer(kind=4)                          :: ilost

    if(flag==1) then
       call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran, ilost)
       if(ilost==1)flag=-1
    else
       call scatter(cell_gas%v, cell_gas%vth_sq_times_m, cell_gas%vturb, nu_cell, k, nu_ext, iran, xcrit, scatterer_list(flag-1))
    end if


  end subroutine gas_scatter



  !--PEEL--
  function  gas_get_tau(cell_gas, distance_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm at freq. nu_cell
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of scatterers and dust
    ! - distance_cm : the distance along which to compute tau [cm]
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! OUTPUTS:
    ! - gas_get_tau : the total optical depth
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)    :: cell_gas
    real(kind=8),intent(in) :: distance_cm
    real(kind=8),intent(in) :: nu_cell
    real(kind=8)            :: gas_get_tau
    integer(kind=4)         :: i

    ! compute optical depths for different components of the gas.
    gas_get_tau = 0d0
    do i = 1, nscatterer
       gas_get_tau = gas_get_tau + get_tau(cell_gas%number_density(scatterer_ion_index(i)) ,cell_gas%vth_sq_times_m, cell_gas%vturb, distance_cm, nu_cell, scatterer_list(i))
    end do

    if(.not.(ignoreDust)) gas_get_tau = gas_get_tau + get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)
    
    return
    
  end function gas_get_tau



  function gas_peeloff_weight(flag,cell_gas,nu_ext,kin,kout,iran)

    integer(kind=4),intent(in)            :: flag
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: nu_ext
    real(kind=8),dimension(3), intent(in) :: kin, kout
    integer(kind=4),intent(inout)         :: iran
    real(kind=8)                          :: gas_peeloff_weight

    if(flag==1) then
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    else
       gas_peeloff_weight = peeloff_weight(cell_gas%v, cell_gas%vth_sq_times_m, cell_gas%vturb, nu_ext, kin, kout, iran, scatterer_list(flag-1))
    end if

  end function gas_peeloff_weight
  !--LEEP--



  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,j,nleaf
    
    nleaf = size(g)
    
    write(unit) (g(i)%v(:), i=1,nleaf)
    do i=1,element_number
       write(unit) (g(j)%number_density(i), j=1,nleaf)
    end do
    write(unit) (g(i)%vth_sq_times_m, i=1,nleaf)
    write(unit) (g(i)%vturb, i=1,nleaf)
    write(unit) (g(i)%ndust, i=1,nleaf)
    write(unit) box_size_cm
  end subroutine dump_gas


  subroutine read_gas(unit,n,g)
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i,j
    
    allocate(g(1:n))
    do i=1,n
       allocate(g(i)%number_density(element_number))
    end do
    
    read(unit) (g(i)%v(:),i=1,n)
    do j=1,element_number
       read(unit) (g(i)%number_density(j), i=1,n)
    end do
    read(unit) (g(i)%vth_sq_times_m, i=1,n)
    read(unit) (g(i)%vturb, i=1,n)
    read(unit) (g(i)%ndust,i=1,n)
    read(unit) box_size_cm

  end subroutine read_gas


  subroutine gas_destructor(g)
    type(gas),dimension(:),allocatable,intent(inout) :: g
    deallocate(g)
  end subroutine gas_destructor


  subroutine set_elements_index

    ! From the list of lines given by the user, this subroutine computes:
    ! - element_number, the number of distinct ions (elements)
    ! - scatterer_ion_index(), which is a list of length nscatterer, assigning an ion to each scatterer(line). Exemple, if the list if SiII1260 - CII1334 - SiII1526, then scatterer_ion_index = (1,2,1)
    ! - name_ions_no_repetition, which is a list of strings of length element_number, giving the names of each ion present, in order of apparition, without repetition
    ! - metal_number, the number of elements (ions) minus hydrogen and deuterium if they are present

    implicit none

    integer(kind=4)              :: i,j, increment,first_index_element(100)
    logical                      :: ion_repeat, HI_present

    if(nscatterer == 0) then
       
       element_number = 0
       metal_number   = 0

    else

       HI_present = .false.
       DI_present = .false.

       ! The first item of the list of scatterers defines the first element.
       element_number = 1
       allocate(scatterer_ion_index(nscatterer))
       scatterer_ion_index(1) = 1
       if(scatterer_list(1)%name_ion == 'HI') HI_present = .true.
       if(scatterer_list(1)%name_ion == 'DI') DI_present = .true.
       first_index_element(1) = 1

       increment = 1

       ! For each line (scatterer) after the first one, checking if the corresponding ion already appeared before.
       do i=2,nscatterer
          if(scatterer_list(i)%name_ion == 'HI') HI_present = .true.
          if(scatterer_list(i)%name_ion == 'DI') DI_present = .true.
          ion_repeat = .false.
          do j=1,i-1
             if(scatterer_list(i)%name_ion == scatterer_list(j)%name_ion) then   ! if true, the ion of scatterer i already appeared before, for example SiII1260 - SiII1526
                scatterer_ion_index(i) = j                                       ! assigns the ion j to the scatterer i.
                ion_repeat = .true.
             end if
          end do
          ! If the ion of scatterer i never appeared, increment the number of different ions.
          if(.not. ion_repeat) then
             element_number = element_number + 1
             first_index_element(1+increment) = i   ! Saving where the new element appeared in the list of scatterers
             scatterer_ion_index(i) = 1 + increment 
             increment = increment + 1
          end if
       end do

       if(DI_present .and. .not. HI_present) then
          print*,'Problem, deuterium cannot be used without HI, stopping the program.'
          stop
       end if


       allocate(name_ions_no_repetition(element_number))
       do i=1,element_number
          name_ions_no_repetition(i) = scatterer_list(first_index_element(i))%name_ion
       end do

       ! Defining metal_number
       if(HI_present) then
          metal_number = element_number-1
          if(DI_present) metal_number = metal_number - 1
       else
          metal_number = element_number
       end if

    end if


    return

  end subroutine set_elements_index


  ! Routine which returns a list of the names of all metallic elements, in order of appearance in the parameter file.
  subroutine get_metal_ion_names(metal_ion_names)

    implicit none

    character(10),allocatable, intent(inout) :: metal_ion_names(:)
    integer(kind=4)                          :: i, j

    allocate(metal_ion_names(metal_number))

    j=1
    do i=1,element_number
       if(name_ions_no_repetition(i) /= 'HI' .and. name_ions_no_repetition(i) /= 'DI') then
          metal_ion_names(j) = name_ions_no_repetition(i)
          j=j+1
       end if
    end do

  end subroutine get_metal_ion_names



  subroutine read_gas_composition_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (HI, D, dust models)
    ! ---------------------------------------------------------------------------------

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
       if (line(1:17) == '[gas_composition]') then
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
          case ('nscatterer')
             read(value,*) nscatterer
             allocate(scatterer_names(nscatterer), scatterer_list(nscatterer))
          case ('scatterer_names')
             read(value,*) scatterer_names(:)
          case ('atomic_data_dir')
             write(atomic_data_dir,'(a)') trim(value)
          case ('krome_data_dir')
             write(krome_data_dir,'(a)') trim(value)
          case ('f_ion')
             read(value,*) f_ion
          case ('Zref')
             read(value,*) Zref
          case ('vturb_kms')
             read(value,*) vturb_kms
          case ('ignoreDust')
             read(value,*) ignoreDust
          ! case ('gas_overwrite')
          !    read(value,*) gas_overwrite
          ! case ('fix_nscat')
          !    read(value,*) fix_nscat
          ! case ('fix_vth_sq_times_m')
          !    read(value,*) fix_vth_sq_times_m
          ! case ('fix_ndust')
          !    read(value,*) fix_ndust
          ! case ('fix_vel')
          !    read(value,*) fix_vel
          ! case ('fix_box_size_cm')
          !    read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    ! define flag for HI-Lyman-alpha if present
    HI_1216_flag = -1
    do i = 1, nscatterer
       if (trim(scatterer_names(i)) == 'HI-1216') then
          HI_1216_index = i
          HI_1216_flag  = i + 1  ! +1 is for dust -> scatterer flags go from 2 to nscatterer + 1.
       end if
    end do
    
    ! Reading scatterer properties
    do i=1,nscatterer
       call read_scatterer_params(trim(atomic_data_dir)//'/'//trim(scatterer_names(i))//'.dat', pfile, scatterer_list(i))
    end do

    ! Calling a routine to compute the number of metallic ions, and others helpful variable defined at the beginning of this routine:
    call set_elements_index
    
    call read_dust_params(pfile)

    return

  end subroutine read_gas_composition_params


  subroutine print_gas_composition_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    integer(kind=4)                     :: i
             
    if (present(unit)) then 
       write(unit,'(a,a,a)') '[gas_composition]'
       write(unit,'(a)')        '# mixture parameters'
       write(unit,'(a,i2)')     '  nscatterer      = ',nscatterer
       do i=1,nscatterer
          write(unit,'(a,i2,a,a)')'  scatterer ',i,'   = ',scatterer_names(i)
       end do
       write(unit,'(a,a)')      '  atomic_data_dir = ',trim(atomic_data_dir)
       write(unit,'(a,a)')      '  krome_data_dir  = ',trim(krome_data_dir)
       write(unit,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref            = ',Zref
       write(unit,'(a,ES10.3)') '  vturb_kms       = ',vturb_kms
       write(unit,'(a,L1)')     '  ignoreDust      = ',ignoreDust
       ! write(unit,'(a)')        '# overwrite parameters'
       ! write(unit,'(a,L1)')     '  gas_overwrite   = ',gas_overwrite
       ! if(gas_overwrite)then
       !    do i=1,nscatterer
       !       write(unit,'(a,i2,a,ES10.3)') '  fix_nscat_',i,' = ',fix_nscat(i)
       !    end do
       !    write(unit,'(a,ES10.3)') '  fix_vth_sq_times_m = ',fix_vth_sq_times_m
       !    write(unit,'(a,ES10.3)') '  fix_ndust       = ',fix_ndust
       !    write(unit,'(a,ES10.3)') '  fix_vel         = ',fix_vel
       !    write(unit,'(a,ES10.3)') '  fix_vturb_kms   = ',fix_vturb_kms
       !    write(unit,'(a,ES10.3)') '  fix_box_size_cm = ',fix_box_size_cm
       ! endif
       ! scatterer params
       write(unit,'(a)')             ' '
       do i=1,nscatterer
          write(unit,'(a,a,a)') '[',trim(scatterer_names(i)),']'
          call print_scatterer_params(scatterer_list(i),unit)
       end do
       call print_dust_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,i2)')     '  nscatterer      = ',nscatterer
       do i=1,nscatterer
          write(*,'(a,i2,a,a)')'  scatterer ',i,' = ',scatterer_names(i)
       end do
       write(*,'(a,a)')      '  atomic_data_dir = ',trim(atomic_data_dir)
       write(*,'(a,a)')      '  krome_data_dir  = ',trim(krome_data_dir)
       write(*,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(*,'(a,ES10.3)') '  Zref            = ',Zref
       write(*,'(a,ES10.3)') '  vturb_kms       = ',vturb_kms
       write(*,'(a,L1)')     '  ignoreDust      = ',ignoreDust
       ! write(*,'(a)')        '# overwrite parameters'
       ! write(*,'(a,L1)')     '  gas_overwrite   = ',gas_overwrite
       ! if(gas_overwrite)then
       !    do i=1,nscatterer
       !       write(*,'(a,i2,a,ES10.3)') '  fix_nscat_',i,' = ',fix_nscat(i)
       !    end do
       !    write(*,'(a,ES10.3)') '  fix_vth_sq_times_m = ',fix_vth_sq_times_m
       !    write(*,'(a,ES10.3)') '  fix_ndust       = ',fix_ndust
       !    write(*,'(a,ES10.3)') '  fix_vel         = ',fix_vel
       !    write(*,'(a,ES10.3)') '  fix_vturb       = ',fix_vturb
       !    write(*,'(a,ES10.3)') '  fix_box_size_cm = ',fix_box_size_cm
       ! endif

       ! scatterer params
       write(*,'(a)')             ' '
       do i=1,nscatterer
          write(*,'(a,a,a)') '[',trim(scatterer_names(i)),']'
          call print_scatterer_params(scatterer_list(i))
       enddo
       call print_dust_params
    end if

    return

  end subroutine print_gas_composition_params


end module module_gas_composition
