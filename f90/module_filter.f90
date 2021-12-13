module module_filter

  use module_constants, only:clight 

  private

  ! parameters in the [filter] section :
  character(1000) :: FilterFile   ! where filter response can be found.
  real(kind=8)    :: redshift=0d0 ! the redshift of the observation we wish to make

  logical         :: use_filter=.false.  ! set to true if used 

  integer(kind=4) :: filter_nlamb
  real(kind=8)    :: filter_lstart, filter_lend  ! start/end wavelengths in Angstom
  real(kind=8)    :: filter_dl
  real(kind=8),allocatable :: filter_transmission(:)

  public:: filter_response,load_filter,read_filter_params
  public:: use_filter

contains

  function filter_response(nu)

    implicit none 

    real(kind=8),intent(in) :: nu  ! input freq. in Hz
    real(kind=8)    :: lambda      ! input wavelength in Angstom
    real(kind=8)    :: filter_response
    integer(kind=4) :: ilam

    lambda = clight / nu * 1e8  ! Angstrom
    ilam = int((lambda-filter_lstart)/filter_dl) + 1
    if (ilam < 1 .or. ilam > filter_nlamb) then
       filter_response = 0.0d0
    else
       filter_response = filter_transmission(ilam)
    end if
    
    return
    
  end function filter_response
  
    
  subroutine load_filter()

    ! Here we load a filter's response curve which has been prepared for use.
    ! -> it is binned regularly. 
    
    implicit none
    integer(kind=4) :: i 

    open(unit=10,file=trim(FilterFile),status='old',form='formatted')
    read(10,*) ! skip header
    read(10,*) filter_lstart ! observer-frame
    read(10,*) filter_lend 
    ! move to rest-frame
    filter_lstart = filter_lstart / (1.0d0 + redshift)
    filter_lend   = filter_lend / (1.0d0 + redshift)
    read(10,*) filter_nlamb
    allocate(filter_transmission(filter_nlamb))
    filter_dl = (filter_lend - filter_lstart)/(filter_nlamb-1)
    do i = 1,filter_nlamb
       read(10,*) filter_transmission(i)
    end do
    close(10)

    use_filter = .true.
    
    return

  end subroutine load_filter


  subroutine read_filter_params(pfile)
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------
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
       if (line(1:8) == '[filter]') then
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
          case ('FilterFile')
             write(FilterFile,'(a)') trim(value)
          case ('redshift')
             read(value,*) redshift
          end select
       end do
    end if
    close(10)

    call load_filter()
    
    return

  end subroutine read_filter_params

end module module_filter
  
