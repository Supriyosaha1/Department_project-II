module module_ssp_lib
  
  ! Module which contains declarations and routines/functions relative to the SSPs library

  use module_utils, only: locatedb
  use module_constants
  
  private
  ! SSP/SED parameters (read from files)
  integer(kind=4)           :: nagebins
  real(kind=4),allocatable  :: agebins(:)        ! values of these ages (in Gyr in the code, yr in the file ...
  real(kind=4),allocatable  :: mid_agebins(:)    ! age boundaries for SED use
  integer(kind=4)           :: num_mets          ! number of metallicities produced by the SSP model
  real(kind=4),allocatable  :: metallicities(:)  ! array of model metallicities (in absolute unit), in log
  integer(kind=4)           :: nlambda           ! number of wavelengths
  real(kind=8),allocatable  :: lambda(:)         ! wavelength
  real(kind=8),allocatable  :: sed_list(:,:,:)   ! seds of SSPs (lambda,age,met)

  type, public :: SSPgrid
     real(kind=8),allocatable :: grid(:,:,:)
     real(kind=8),allocatable :: age(:),met(:),lambda(:)
     integer(kind=4)          :: nage,nmet,nlambda
  end type SSPgrid

  type, public :: PWLgrid
     real(kind=8),allocatable :: F0(:,:), beta(:,:), Ndot(:,:)
     real(kind=8),allocatable :: age(:),met(:)
     real(kind=8),allocatable :: lmin, lmax, lambda0
     integer(kind=4)          :: nage,nmet     
  end type PWLgrid
  
  type(SSPgrid) :: fullLib
  logical       :: isLoaded = .false.
  real(kind=8)  :: Lsun_in_cgs = 3.826d33          ! erg/s using the same convention as in BC03

  
  public :: init_ssp_lib, ssp_lib_extract_subset, ssp_lib_interpolate, ssp_lib_integrate
  public :: ssp_lib_fit_powerlaw, ssp_lib_interpolate_powerlaw


contains

!*******************************************************************************************************

  subroutine init_ssp_lib(SSPdir,lambdamin,lambdamax)

    implicit none

    character(1000),intent(in)       :: SSPdir
    real(kind=8),intent(in),optional :: lambdamin, lambdamax
    real(kind=8),allocatable         :: nu(:)
    integer(kind=4)                  :: i,j
    
    ! read data
    call read_met_bins(SSPdir)
    call read_age_bins(SSPdir)
    call read_all_seds(SSPdir,lambdamin,lambdamax)

    print*,'reading SSPs...'
    print*,minval(agebins),maxval(agebins)
    print*,minval(metallicities),maxval(metallicities)
    print*,minval(lambda), maxval(lambda)
    print*, nagebins,num_mets,nlambda

    fullLib%nage = nagebins
    fullLib%nmet = num_mets
    fullLib%nlambda = nlambda    
    allocate(fullLib%age(nagebins)) ; fullLib%age(:) = agebins        ! convention x = ages
    allocate(fullLib%met(num_mets)) ; fullLib%met(:) = metallicities  ! convention y = metallicities
    allocate(fullLib%lambda(nlambda)) ; fullLib%lambda(:) = lambda    ! convention z = SEDs
    allocate(fullLib%grid(nagebins,num_mets,nlambda))
    fullLib%grid(:,:,:) =  reshape(source=sed_list,shape=(/nagebins,num_mets,nlambda/),order=(/3,1,2/)) ! convention z = SEDs

    isLoaded=.true.
    
    
    ! deallocate
    deallocate(metallicities, agebins, lambda, sed_list)
    
    ! spectra are in [Lsun / A / Msun]
    ! with Lsun = 3.826e33 erg/s and Msun = 2e33g
    ! lambda are in Angstrom
    ! age are in yr

    ! Change units 
    fullLib%grid = fullLib%grid * Lsun_in_cgs   ! in erg/s/A/Msun
    
  end subroutine init_ssp_lib
    
!*******************************************************************************************************
  subroutine ssp_lib_extract_subset(lmin,lmax,subLib)

    implicit none

    real(kind=8),intent(inout) :: lmin,lmax
    type(SSPgrid), intent(out) :: subLib
    integer(kind=4)            :: i0,i1,nl,i,j
    real(kind=8)               :: lambda0
    real(kind=8),allocatable   :: nu(:)

    if (lmin==lmax) then
       call locatedb(fullLib%lambda, fullLib%nlambda, lmin, i0)
       if( (lmin-fullLib%lambda(i0)) > (fullLib%lambda(i0+1)-lmin))then
          i0 = i0+1
       endif
       lambda0 = fullLib%lambda(i0)
       nl = 1
       i1 = i0
    else
       call locatedb(fullLib%lambda, fullLib%nlambda, lmin, i0)
       if( (lmin-fullLib%lambda(i0)) > (fullLib%lambda(i0+1)-lmin))then
          i0 = i0+1
       endif
       lmin = fullLib%lambda(i0)
       call locatedb(fullLib%lambda, fullLib%nlambda, lmax, i1)
       if( (lmax-fullLib%lambda(i1)) > (fullLib%lambda(i1+1)-lmax))then
          i1 = i1+1
       endif
       lmax = fullLib%lambda(i1)
       nl = i1-i0+1
    endif
    
    subLib%nage = fullLib%nage
    subLib%nmet = fullLib%nmet
    subLib%nlambda = nl
    allocate(subLib%age(subLib%nage))
    allocate(subLib%met(subLib%nmet))
    allocate(subLib%lambda(subLib%nlambda))
    allocate(subLib%grid(subLib%nage,subLib%nmet,subLib%nlambda))
    
    subLib%age = fullLib%age
    subLib%met = fullLib%met
    subLib%lambda = fullLib%lambda(i0:i1)

    !JB- convert ergs to nb of photons
    allocate(nu(nl))
    nu(:) = clight / (subLib%lambda(:)*1e-8) ! [Hz]
    do i = 1, subLib%nlambda
       subLib%grid(:,:,i) = fullLib%grid(:,:,i0+i-1) / planck / nu(i)  ! nb of photons / s / A / Msun 
    enddo
    deallocate(nu)
    !-JB
    
    return
  end subroutine ssp_lib_extract_subset

!*******************************************************************************************************
  subroutine ssp_lib_interpolate(lib, x, y, res)

    implicit none
    type(SSPgrid),intent(in)              :: lib
    real(kind=8), intent(in)              :: x,y
    real(kind=8),dimension(:),intent(out) :: res
    integer(kind=4)                       :: j1,j2,ny,i1,i2,nx
    real(kind=8)                          :: wy1,wy2,norm,wx1,wx2
    
    ny = size(lib%met)
    ! locate bins and compute weights
    if (y <= lib%met(1)) then 
       j1 = 1
       wy1  = 1.0d0
       j2 = 2
       wy2  = 0.0d0
    else if (y >= lib%met(ny)) then 
       j1  = ny-1
       wy1 = 0.0d0
       j2  = ny
       wy2 = 1.0d0
    else
       do j2 = 1,ny
          if (y <= lib%met(j2)) exit
       end do
       j1 = j2 - 1
       wy1  = lib%met(j2) - y 
       wy2  = y - lib%met(j1)
       norm = wy1 + wy2
       wy1  = wy1 / norm
       wy2  = wy2 / norm
    end if
    

    nx = size(lib%age)
    ! locate bins and compute weights
    if (x <= lib%age(1)) then 
       i1 = 1
       wx1 = 1.0d0
       i2 = 2
       wx2 = 0.0d0
    else if (x >= lib%age(nx)) then 
       i1  = nx-1
       wx1 = 0.0d0
       i2  = nx
       wx2 = 1.0d0
    else
       do i2 = 1,nx
          if (x <= lib%age(i2)) exit
       end do
       i1 = i2 - 1
       wx1  = lib%age(i2) - x
       wx2  = x - lib%age(i1)
       norm = wx1 + wx2
       wx1  = wx1 / norm
       wx2  = wx2 / norm
    end if
    
    res(:) =  wx2 * wy2 * lib%grid(i2, j2, :)  &
            + wx1 * wy2 * lib%grid(i1, j2, :)  &
            + wx2 * wy1 * lib%grid(i2, j1, :)  &
            + wx1 * wy1 * lib%grid(i1, j1, :)
    return
  end subroutine ssp_lib_interpolate

  !*******************************************************************************************************

  subroutine ssp_lib_fit_powerlaw(lminfit,lmaxfit,lambdamin,lambdamax,clip,fulllib_powerlaws)

    ! fit a powerlaw over the spectral range [lminfit,lmaxfit] of spectra at each metallicity and age,
    ! and integrate the powerlaw from lambdamin to lambdamax to get photon rate. 
    ! The powerlaw fit has the for F_lambda = F_0 * (lambda / lambda_0) ** beta
    ! the results are stored in the module variable fulllib_powerlaws

    implicit none

    real(kind=8),intent(in)   :: lminfit,lmaxfit,lambdamin,lambdamax
    logical,intent(in)        :: clip
    type(PWLgrid),intent(out) :: fulllib_powerlaws
    integer(kind=4)           :: iage,imet,ilam,i0,i1,nl
    real(kind=8),allocatable  :: f(:),ffit(:),logl(:)
    logical,allocatable       :: mask(:)
    real(kind=8)              :: a,b,lambda_0_cm

    ! define basic elements
    fulllib_powerlaws%lmin    = lambdamin
    fulllib_powerlaws%lmax    = lambdamax
    fulllib_powerlaws%lambda0 = 0.5d0 * (lambdamin + lambdamax)  ! Angstrom
    fulllib_powerlaws%nmet    = fulllib%nmet
    fulllib_powerlaws%nage    = fulllib%nage
    allocate(fulllib_powerlaws%met(fulllib_powerlaws%nmet),fulllib_powerlaws%age(fulllib_powerlaws%nage))
    fulllib_powerlaws%met     = fulllib%met
    fulllib_powerlaws%age     = fulllib%age
    ! allocate fit parameters 
    allocate(fulllib_powerlaws%F0(fulllib_powerlaws%nage,fulllib_powerlaws%nmet))
    allocate(fulllib_powerlaws%beta(fulllib_powerlaws%nage,fulllib_powerlaws%nmet))
    ! allocate integrated photon rate
    allocate(fulllib_powerlaws%Ndot(fulllib_powerlaws%nage,fulllib_powerlaws%nmet))
    
    ! find wavelength indexes corresponding to lambdamin and lambdamax in fulllib.
    call locatedb(fullLib%lambda, fullLib%nlambda, lminfit, i0)
    if( (lminfit-fullLib%lambda(i0)) > (fullLib%lambda(i0+1)-lminfit))then
       i0 = i0+1
    endif
    call locatedb(fullLib%lambda, fullLib%nlambda, lmaxfit, i1)
    if( (lmaxfit-fullLib%lambda(i1)) > (fullLib%lambda(i1+1)-lmaxfit))then
       i1 = i1+1
    endif
    nl = i1-i0+1
    if (nl < 20) then  ! require more than 20 points for a fit... 
       print*,'ERROR: the wavelength interval should be larger in ssp_lib_fit_powerlaw'
       stop
    end if

    print*,i0,i1
    print*,lminfit,fullLib%lambda(i0)
    print*,lmaxfit,fullLib%lambda(i1)
    
    allocate(f(nl),mask(nl),ffit(nl),logl(nl))
    logl = log10(fulllib%lambda(i0:i1)/fulllib_powerlaws%lambda0)
    ! loop on age
    do iage = 1, fulllib_powerlaws%nage
       ! loop on metallicity 
       do imet = 1, fulllib_powerlaws%nmet
          f = log10(fulllib%grid(iage,imet,i0:i1))  ! local copy of spectrum, in erg/s/A
          ! compute fit coefficients -> Log10(F) = a + b*Log10(lambda)
          mask = .True.
          call get_fit_coeffs(nl,logl,f,mask,a,b)
          if (clip) then 
             ! remove big absorption lines (points much lower than fit) and recompute fit.
             ffit = a + b * logl
             mask = (f-ffit > -0.097)  ! ie flux / fit > 0.8
             call get_fit_coeffs(nl,logl,f,mask,a,b)
             ! second pass with stronger constraint
             ffit = a + b * logl
             mask = (f-ffit > -0.0223)  ! ie flux / fit > 0.95
             call get_fit_coeffs(nl,logl,f,mask,a,b)
          end if
          ! Extract F_0 and beta from a,b
          fulllib_powerlaws%beta(iage,imet) = b
          fulllib_powerlaws%F0(iage,imet)  = 10.0d0**(a) ! + b * log10(fulllib_powerlaws%lambda0))
       end do
    end do
    deallocate(f,mask,ffit,logl)

    ! integrate power-laws to get nb of photons per sec on wavelength range
    lambda_0_cm = fulllib_powerlaws%lambda0 * 1e-8
    a = fulllib_powerlaws%lambda0
    do iage = 1, fulllib_powerlaws%nage
       do imet = 1, fulllib_powerlaws%nmet
          ! Given that F_lbda = F_0 (lbda / lbda_0)**beta (in erg/s/A),
          ! the number of photons (in /s/A) is N_lbda = F_0*lbda_0/hc * (lbda/lbda_0)**(1+beta).
          ! (NB: the first lbda_0 here has to be in cm)
          ! This integrates to (in #/s) :
          ! (F_0 lbda_0 / hc) * lbda_0/(beta+2)  * [ (lbda_max/lbda_0)**(2+beta) - (lbda_min/lbda_0)**(2+beta)]
          ! (NB: the first lbda_0 here is in cm, the second in A). 
          ! OR, if beta == -2, the integral is
          ! (F_0*lbda_0/hc) * lbda_0 * ln(lbda_max/lbda_min)     [again, first lbda_0 in cm, second in A]
          fulllib_powerlaws%Ndot(iage,imet) = fulllib_powerlaws%F0(iage,imet) * lambda_0_cm / planck / clight  ! phot/s/A
          b = fulllib_powerlaws%beta(iage,imet)
          if (b == -2.0d0) then
             fulllib_powerlaws%Ndot(iage,imet) = fulllib_powerlaws%Ndot(iage,imet) * a * log(lambdamax/lambdamin)
          else
             fulllib_powerlaws%Ndot(iage,imet) = fulllib_powerlaws%Ndot(iage,imet) * a/(b+2)*( (lambdamax/a)**(2+b) - (lambdamin/a)**(2+b) )
          end if
       end do
    end do

    return

  contains
    
    subroutine get_fit_coeffs(n,l,f,mask,a,b)
      ! smart-ass algebra ... 
      implicit none
      integer(kind=4),intent(in) :: n
      real(kind=8),intent(in)    :: l(n),f(n)
      logical,intent(in)         :: mask(n)
      real(kind=8),intent(out)   :: a,b
      real(kind=8) :: Xbar,Dbar,XDbar,X2bar,neff
      real(kind=8) :: ones(n)

      ones = 1.0d0 
      neff = sum(ones,mask=mask)
      Xbar = sum(l,mask=mask)/neff
      Dbar = sum(f,mask=mask)/neff
      XDbar = sum(f*l,mask=mask)/neff
      X2bar = sum(l*l,mask=mask)/neff
      a = (Xbar*XDbar - Dbar*X2bar) / (Xbar*Xbar - X2bar)
      b = (Dbar - a)/Xbar
      return
    end subroutine get_fit_coeffs

  end subroutine ssp_lib_fit_powerlaw

  !*******************************************************************************************************
!!$
!!$  subroutine ssp_lib_integrate_powerlaw(lmin,lmax,lambda_0,F_0,beta,Ndot)
!!$
!!$    implicit none
!!$
!!$    real(kind=8),intent(in)              :: lmin,lmax ! min and max wavelength over which to integrate. [A]
!!$    real(kind=8),intent(out)             :: lambda_0
!!$    real(kind=8),allocatable,intent(out) :: F_0(:,:), beta(:,:)
!!$    integer(kind=4) :: iage,imet
!!$    real(kind=8)    ::lambda_0_cm
!!$
!!$    lambda_0_cm = lambda_0 * 1e-8
!!$    
!!$    ! loop on age
!!$    do iage = 1, lib%nage
!!$       ! loop on metallicity 
!!$       do imet = 1, lib%nmet
!!$          ! Given that F_lbda = F_0 (lbda / lbda_0)**beta (in erg/s/A),
!!$          ! the number of photons (in /s/A) is N_lbda = F_0*lbda_0/hc * (lbda/lbda_0)**(1+beta).
!!$          ! (NB: the first lbda_0 here has to be in cm)
!!$          ! This integrates to (in #/s) :
!!$          ! (F_0 lbda_0 / hc) * lbda_0/(beta+2)  * [ (lbda_max/lbda_0)**(2+beta) - (lbda_min/lbda_0)**(2+beta)]
!!$          ! (NB: the first lbda_0 here is in cm, the second in A). 
!!$          ! OR, if beta == -2, the integral is
!!$          ! (F_0*lbda_0/hc) * lbda_0 * ln(lbda_max/lbda_min)     [again, first lbda_0 in cm, second in A]
!!$          Ndot(iage,imet) = F_0(iage,imet) * lambda_0_cm / planck / clight  ! phot/s/A
!!$          if (beta == -2.0d0) then
!!$             Ndot(iage,imet) = Ndot(iage,imet) * lambda_0 * log(lmax/lmin)
!!$          else
!!$             Ndot(iage,imet) = Ndot(iage,imet) * lambda_0 / (beta(iage,imet)+2) * ((lmax/lambda_0)**(2+beta) - (lmin/lambda_0)**(2+beta))
!!$          end if
!!$       end do
!!$    end do
!!$    
!!$    return
!!$  end subroutine ssp_lib_integrate_powerlaw
!!$  
  !*******************************************************************************************************
  
  subroutine ssp_lib_interpolate_powerlaw(lib,x,y,res,i1,j1,wy1,wx1)
    implicit none
    type(PWLgrid),intent(in)    :: lib
    real(kind=8),intent(in)     :: x,y
    real(kind=8),intent(out)    :: res  ! nb of phots per sec.
    integer(kind=4),intent(out) :: i1,j1
    real(kind=8),intent(out)    :: wy1,wx1
    integer(kind=4)             :: ny,nx,j2,i2
    real(kind=8)                :: norm,wx2,wy2
    
    ny = size(lib%met)
    ! locate bins and compute weights
    if (y <= lib%met(1)) then 
       j1 = 1
       wy1  = 1.0d0
       j2 = 2
       wy2  = 0.0d0
    else if (y >= lib%met(ny)) then 
       j1  = ny-1
       wy1 = 0.0d0
       j2  = ny
       wy2 = 1.0d0
    else
       do j2 = 1,ny
          if (y <= lib%met(j2)) exit
       end do
       j1 = j2 - 1
       wy1  = lib%met(j2) - y 
       wy2  = y - lib%met(j1)
       norm = wy1 + wy2
       wy1  = wy1 / norm
       wy2  = wy2 / norm
    end if
    
    nx = size(lib%age)
    ! locate bins and compute weights
    if (x <= lib%age(1)) then 
       i1 = 1
       wx1 = 1.0d0
       i2 = 2
       wx2 = 0.0d0
    else if (x >= lib%age(nx)) then 
       i1  = nx-1
       wx1 = 0.0d0
       i2  = nx
       wx2 = 1.0d0
    else
       do i2 = 1,nx
          if (x <= lib%age(i2)) exit
       end do
       i1 = i2 - 1
       wx1  = lib%age(i2) - x
       wx2  = x - lib%age(i1)
       norm = wx1 + wx2
       wx1  = wx1 / norm
       wx2  = wx2 / norm
    end if
    
    res =  wx2 * wy2 * lib%Ndot(i2, j2)  &
         + wx1 * wy2 * lib%Ndot(i1, j2)  &
         + wx2 * wy1 * lib%Ndot(i2, j1)  &
         + wx1 * wy1 * lib%Ndot(i1, j1)
 
    return
    
  end subroutine ssp_lib_interpolate_powerlaw

  !*******************************************************************************************************

  subroutine ssp_lib_integrate(x,y,n,integral)

    implicit none 

    real(kind=8),intent(in)  :: x(n), y(n)
    integer,intent(in)       :: n
    real(kind=8),intent(out) :: integral
    integral = trapz1(x,y,n)

    return

  end subroutine ssp_lib_integrate
  
  !*************************************************************************

  FUNCTION trapz1(X,Y,N,cum)
    
    ! Integrates function Y(X) along the whole interval 1..N, using a very 
    ! simple staircase method and returns the result.
    ! Optionally, the culumative integral is returned in the cum argument.
    !-------------------------------------------------------------------------
    implicit none 
    integer :: N,i
    real(kind=8):: trapz1
    real(kind=8):: X(N),Y(N)
    real(kind=8),optional::cum(N)
    real(kind=8),allocatable::cumInt(:)
    !-------------------------------------------------------------------------
    allocate(cumInt(N))
    cumInt(:)=0.d0
    trapz1=0.d0
    if (N.le.1) RETURN
    do i=2,N
       cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
    end do
    trapz1 = cumInt(N)
    if(present(cum)) cum=cumInt
    deallocate(cumInt)
    return    
  END FUNCTION trapz1
!*************************************************************************

  subroutine read_met_bins(SSPdir)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    character(1000) :: metsbin_file
    integer(kind=4) :: imet

    write(metsbin_file,'(a,a)') trim(SSPdir),"/metallicity_bins.dat"
    call test_stop(metsbin_file)
    open(unit=12,file=metsbin_file,status='old',form='formatted')
    read(12,'(i8)') num_mets
    allocate(metallicities(num_mets))
    do imet = 1,num_mets
       read(12,'(e14.6)') metallicities(imet)
    end do
    close(12)

    ! take log to make interpolations in log
    metallicities  = log10(metallicities)

    return

  end subroutine read_met_bins

!*******************************************************************************************************

  subroutine read_all_seds(SSPdir,lambdamin,lambdamax)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    real(kind=8),intent(in),optional :: lambdamin,lambdamax
    character(1000)           :: sed_file
    integer(kind=4)           :: dum,iage,imet
    integer(kind=4)           :: ilambdamin,ilambdamax
    real(kind=8),allocatable  :: sed(:)

    write(sed_file,'(a,a)') trim(SSPdir),"/all_seds.dat"
    call test_stop(sed_file)
    open(unit=12,file=sed_file,status='old',form='unformatted')
    read(12) nlambda,dum
    allocate(sed(nlambda))
    read(12) sed(:)

    if(present(lambdamin).and.present(lambdamax))then
       ! extract only desired wavelength range 
       call locatedb(sed,nlambda,lambdamin-10.,ilambdamin)  
       call locatedb(sed,nlambda,lambdamax+10.,ilambdamax)  
       nlambda = ilambdamax - ilambdamin + 1
       allocate(lambda(nlambda))
       lambda = sed(ilambdamin:ilambdamax)
    else
       allocate(lambda(nlambda))
       lambda(:) = sed(:)
       ilambdamin = 1
       ilambdamax = nlambda
    endif
    allocate(sed_list(nlambda,nagebins,num_mets))
    
    do imet = 1,num_mets
       do iage = 1,nagebins
          read(12) sed
          sed_list(:,iage,imet) = sed(ilambdamin:ilambdamax)
       end do
    end do
    close(12)

    deallocate(sed)
    
    return

  end subroutine read_all_seds

!*******************************************************************************************************

  subroutine read_age_bins(SSPdir)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    integer(kind=4) :: iage
    character(1000) :: filename
    
    write(filename,'(a,a)') trim(SSPdir),"/age_bins.dat"
    call test_stop(filename)
    open(unit=12,file=filename,status='old',form='formatted')
    read(12,'(i8)') nagebins
    allocate(agebins(nagebins))
    do iage = 1,nagebins
       read(12,'(e14.6)') agebins(iage)
       agebins(iage) = agebins(iage) * 1.e-9  ! convert from yr to Gyr
    end do
    close(12)
    
    ! define age boundaries where to use each SSP-SED
    allocate(mid_agebins(nagebins+1)) 
    mid_agebins(1)    = agebins(1) ! == 0 yr
    do iage = 2,nagebins
       mid_agebins(iage) = 0.5*(agebins(iage)+agebins(iage-1))
    end do
    mid_agebins(nagebins+1) = agebins(nagebins) ! == 20 Gyr

    return

  end subroutine read_age_bins

!*******************************************************************************************************

  subroutine test_stop(file)

    implicit none

    character(*) :: file
    logical      :: exist

    inquire(file=file,exist=exist)
    if (.not.exist) then
       write(*,*) '> Error: ',trim(file),' does not exist.'  
       stop
    endif

    return

  end subroutine test_stop

!*******************************************************************************************************

end module module_ssp_lib
