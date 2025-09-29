program GenerateAMRmodel

  use module_domain
  use module_mesh
  use module_gas_composition
  use module_idealised_model

  implicit none

  type(domain)                             :: domaine_de_calcul
  type(domain),dimension(:),allocatable    :: domain_list
  type(mesh)                               :: domain_mesh

  character(2000)                          :: toto,meshroot,parameter_file,fichier, fichier2
  character(2000),dimension(:),allocatable :: domain_file_list, mesh_file_list
  real(kind=8)                             :: computdom_max_x, computdom_max_y, computdom_max_z
  real(kind=8)                             :: decompdom_max_x, decompdom_max_y, decompdom_max_z
  real(kind=8)                             :: computdom_min_x, computdom_min_y, computdom_min_z
  real(kind=8)                             :: decompdom_min_x, decompdom_min_y, decompdom_min_z
  real(kind=8)                             :: start, finish, rate
  integer(kind=4)                          :: i, narg, noctmax
  integer(kind=8)                          :: c1 , c2, cr
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [GenerateAMRmodel] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: DomDumpDir = 'test/'       ! directory to which outputs will be written

  ! --- computational domain  
  character(10)             :: comput_dom_type      = 'sphere'         ! shape type of domain  // default is a shpere.
  real(kind=8),dimension(3) :: comput_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: comput_dom_rsp       = 0.3              ! radius of sphere [code units]
  real(kind=8)              :: comput_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: comput_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: comput_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: comput_dom_thickness = 0.1              ! thickness of slab [code units]
  ! --- domain decomposition // Only works with a number of domains, sharing the same type.
  character(10)             :: decomp_dom_type = 'sphere'              ! shape type of domain  // default is one sphere.
  integer(kind=4)           :: decomp_dom_ndomain = 1                  ! nb of domains in decomposition
  real(kind=8),allocatable  :: decomp_dom_xc(:)                        ! x center of domain(s) [code units], default value 0.5
  real(kind=8),allocatable  :: decomp_dom_yc(:)                        ! y center of domain(s) [code units], default value 0.5
  real(kind=8),allocatable  :: decomp_dom_zc(:)                        ! z center of domain(s) [code units], default value 0.5
  real(kind=8),allocatable  :: decomp_dom_rsp(:)                       ! radius of sphere [code units]
  real(kind=8),allocatable  :: decomp_dom_rin(:)                       ! inner radius of shell [code units]
  real(kind=8),allocatable  :: decomp_dom_rout(:)                      ! outer radius of shell [code units]
  real(kind=8),allocatable  :: decomp_dom_size(:)                      ! size of cube [code units]
  real(kind=8),allocatable  :: decomp_dom_thickness(:)                 ! thickness of slab [code units]
  ! --- miscelaneous
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------
  
  call cpu_time(start)
  call system_clock(count_rate=cr)
  rate = float(cr)
  call system_clock(c1)
  
  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: GenerateAMRmodel params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_GenerateAMRmodel_params(parameter_file)
  if (verbose) call print_GenerateAMRmodel_params
  ! ------------------------------------------------------------

  ! idealised model do not work with multiple ions/elements
  if(element_number > 1)then
     print*,'element_number =',element_number
     print*,'idealised model works for at most one element/ion, check the config file...'
     stop
  endif

  
  ! Define the computational domain. This domain describes the volume in which photons fly.
  select case(comput_dom_type)
  case('sphere')
     call domain_constructor_from_scratch(domaine_de_calcul,comput_dom_type, &
          xc=comput_dom_pos(1),yc=comput_dom_pos(2),zc=comput_dom_pos(3),r=comput_dom_rsp)
     computdom_max_x = comput_dom_pos(1)+comput_dom_rsp
     computdom_max_y = comput_dom_pos(2)+comput_dom_rsp
     computdom_max_z = comput_dom_pos(3)+comput_dom_rsp
     computdom_min_x = comput_dom_pos(1)-comput_dom_rsp
     computdom_min_y = comput_dom_pos(2)-comput_dom_rsp
     computdom_min_z = comput_dom_pos(3)-comput_dom_rsp
  case('shell')
     call domain_constructor_from_scratch(domaine_de_calcul,comput_dom_type, &
          xc=comput_dom_pos(1),yc=comput_dom_pos(2),zc=comput_dom_pos(3),r_inbound=comput_dom_rin,r_outbound=comput_dom_rout)
     computdom_max_x = comput_dom_pos(1)+comput_dom_rout
     computdom_max_y = comput_dom_pos(2)+comput_dom_rout
     computdom_max_z = comput_dom_pos(3)+comput_dom_rout
     computdom_min_x = comput_dom_pos(1)-comput_dom_rout
     computdom_min_y = comput_dom_pos(2)-comput_dom_rout
     computdom_min_z = comput_dom_pos(3)-comput_dom_rout
  case('cube')
     call domain_constructor_from_scratch(domaine_de_calcul,comput_dom_type, & 
          xc=comput_dom_pos(1),yc=comput_dom_pos(2),zc=comput_dom_pos(3),size=comput_dom_size)
     computdom_max_x = comput_dom_pos(1)+comput_dom_size/2.0d0
     computdom_max_y = comput_dom_pos(2)+comput_dom_size/2.0d0
     computdom_max_z = comput_dom_pos(3)+comput_dom_size/2.0d0
     computdom_min_x = comput_dom_pos(1)-comput_dom_size/2.0d0
     computdom_min_y = comput_dom_pos(2)-comput_dom_size/2.0d0
     computdom_min_z = comput_dom_pos(3)-comput_dom_size/2.0d0
  case('slab')
     call domain_constructor_from_scratch(domaine_de_calcul,comput_dom_type, &
          xc=comput_dom_pos(1),yc=comput_dom_pos(2),zc=comput_dom_pos(3),thickness=comput_dom_thickness)
     computdom_max_x = 1.0d0 ; computdom_max_y = 1.0d0
     computdom_max_z = comput_dom_pos(3)+comput_dom_thickness/2.0d0
     computdom_min_x = 0.0d0 ; computdom_min_y = 0.0d0
     computdom_min_z = comput_dom_pos(3)-comput_dom_thickness/2.0d0
  end select
  
  
  
  ! domain decomposition 
  if (verbose) print *,'Building domains...'
  meshroot = 'domain_'
  allocate(domain_list(decomp_dom_ndomain))
  allocate(domain_file_list(decomp_dom_ndomain),mesh_file_list(decomp_dom_ndomain))
  decompdom_max_x = 0.0d0 ; decompdom_max_y = 0.0d0 ; decompdom_max_z = 0.0d0
  decompdom_min_x = 1.0d0 ; decompdom_min_y = 1.0d0 ; decompdom_min_z = 1.0d0
  do i = 1, decomp_dom_ndomain
     select case(decomp_dom_type)
     case('sphere')
        call domain_constructor_from_scratch(domain_list(i),decomp_dom_type, &
             xc=decomp_dom_xc(i),yc=decomp_dom_yc(i),zc=decomp_dom_zc(i),r=decomp_dom_rsp(i))
        decompdom_max_x = max(decompdom_max_x,decomp_dom_xc(i)+decomp_dom_rsp(i))
        decompdom_max_y = max(decompdom_max_y,decomp_dom_yc(i)+decomp_dom_rsp(i))
        decompdom_max_z = max(decompdom_max_z,decomp_dom_zc(i)+decomp_dom_rsp(i))
        decompdom_min_x = min(decompdom_min_x,decomp_dom_xc(i)-decomp_dom_rsp(i))
        decompdom_min_y = min(decompdom_min_y,decomp_dom_yc(i)-decomp_dom_rsp(i))
        decompdom_min_z = min(decompdom_min_z,decomp_dom_zc(i)-decomp_dom_rsp(i))
     case('shell')
        call domain_constructor_from_scratch(domain_list(i),decomp_dom_type, &
             xc=decomp_dom_xc(i),yc=decomp_dom_yc(i),zc=decomp_dom_zc(i),r_inbound=decomp_dom_rin(i),r_outbound=decomp_dom_rout(i))
        decompdom_max_x = max(decompdom_max_x,decomp_dom_xc(i)+decomp_dom_rout(i))
        decompdom_max_y = max(decompdom_max_y,decomp_dom_yc(i)+decomp_dom_rout(i))
        decompdom_max_z = max(decompdom_max_z,decomp_dom_zc(i)+decomp_dom_rout(i))
        decompdom_min_x = min(decompdom_min_x,decomp_dom_xc(i)-decomp_dom_rout(i))
        decompdom_min_y = min(decompdom_min_y,decomp_dom_yc(i)-decomp_dom_rout(i))
        decompdom_min_z = min(decompdom_min_z,decomp_dom_zc(i)-decomp_dom_rout(i))
     case('cube')
        call domain_constructor_from_scratch(domain_list(i),decomp_dom_type, & 
             xc=decomp_dom_xc(i),yc=decomp_dom_yc(i),zc=decomp_dom_zc(i),size=decomp_dom_size(i))
        decompdom_max_x = max(decompdom_max_x,decomp_dom_xc(i)+decomp_dom_size(i)/2.0d0)
        decompdom_max_y = max(decompdom_max_y,decomp_dom_yc(i)+decomp_dom_size(i)/2.0d0)
        decompdom_max_z = max(decompdom_max_z,decomp_dom_zc(i)+decomp_dom_size(i)/2.0d0)
        decompdom_min_x = min(decompdom_min_x,decomp_dom_xc(i)-decomp_dom_size(i)/2.0d0)
        decompdom_min_y = min(decompdom_min_y,decomp_dom_yc(i)-decomp_dom_size(i)/2.0d0)
        decompdom_min_z = min(decompdom_min_z,decomp_dom_zc(i)-decomp_dom_size(i)/2.0d0)
     case('slab')
        call domain_constructor_from_scratch(domain_list(i),decomp_dom_type, &
             xc=decomp_dom_xc(i),yc=decomp_dom_yc(i),zc=decomp_dom_zc(i),thickness=decomp_dom_thickness(i))
        decompdom_max_x = 1.0d0 ; decompdom_max_y = 1.0d0
        decompdom_min_x = 0.0d0 ; decompdom_min_y = 0.0d0
        decompdom_max_z = max(decompdom_max_z,decomp_dom_zc(i)+decomp_dom_thickness(i)/2.0d0)
        decompdom_min_z = min(decompdom_min_z,decomp_dom_zc(i)-decomp_dom_thickness(i)/2.0d0)
     end select
     write(toto,'(i8)') i
     write(toto,'(a)') adjustl(toto)  ! remove leading spaces
     write(domain_file_list(i),'(a,a,a)') trim(meshroot),trim(toto),'.dom'
     write(mesh_file_list(i),'(a,a,a)') trim(meshroot),trim(toto),'.mesh'
     if (verbose) write(*,'(a,i3,a,a,a,a)') '     \_',i,'  ',trim(domain_file_list(i)),' ',trim(mesh_file_list(i))
  end do

  ! The computational domain should be fully enclosed in the domain mesh.
  if ((computdom_max_x > decompdom_max_x).or.(computdom_max_y > decompdom_max_y).or.(computdom_max_z > decompdom_max_z).or.&
       (computdom_min_x < decompdom_min_x).or.(computdom_min_y < decompdom_min_y).or.(computdom_min_z < decompdom_min_z))then
     print*,'ERROR: computational domain should be fully enclosed in the data domains.'
     print*,computdom_max_x, decompdom_max_x
     print*,computdom_max_y, decompdom_max_y
     print*,computdom_max_z, decompdom_max_z
     print*,computdom_min_x, decompdom_min_x
     print*,computdom_min_y, decompdom_min_y
     print*,computdom_min_z, decompdom_min_z
     stop
  endif

  ! write master info
  fichier = "compute_domain.dom"
  call domain_write_file(trim(DomDumpDir)//trim(fichier),domaine_de_calcul)
  fichier2 = "domain_decomposition_params.dat"
  open(unit=10, file=trim(DomDumpDir)//trim(fichier2))
  write(10,*) 'computational_domain_file = ',trim(fichier)
  write(10,*) 'Ndomain = ',decomp_dom_ndomain
  do i=1,decomp_dom_ndomain
     write(10,*) 'domain_file = ',trim(domain_file_list(i))
     write(10,*) 'mesh_file   = ',trim(mesh_file_list(i))
     fichier = trim(DomDumpDir)//trim(domain_file_list(i))
     call domain_write_file(fichier,domain_list(i))
  end do
  close(10)
  
  
  
  ! building of the meshes
  do i = 1,decomp_dom_ndomain
     if(verbose)print*,' '

     noctmax = 40000000
     if (verbose) write(*,*)'Building the mesh from the model...'     
     call mesh_from_model(domain_list(i),domain_mesh,noctmax)

     fichier = trim(DomDumpDir)//trim(mesh_file_list(i))
     call dump_mesh(domain_mesh, fichier)
     call mesh_destructor(domain_mesh)

  enddo

  call cpu_time(finish)
  print '(" --> Done with GenerateAMRmodel. Total elapsed time = ",f12.3," seconds.")',finish-start
  call system_clock(c2)
  print '("                               system_clock time = ",f12.3," seconds.")',(c2-c1)/rate
  print*,' '
  
contains

  
  subroutine read_GenerateAMRmodel_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(10000)        :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present
    logical                 :: ndomain_present 
    
    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:18) == '[GenerateAMRmodel]') then
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
          case('comput_dom_type')
             write(comput_dom_type,'(a)') trim(value)
          case ('comput_dom_pos')
             read(value,*) comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
          case ('comput_dom_rsp')
             read(value,*) comput_dom_rsp
          case ('comput_dom_rin')
             read(value,*) comput_dom_rin
          case ('comput_dom_rout')
             read(value,*) comput_dom_rout
          case ('comput_dom_size')
             read(value,*) comput_dom_size
          case ('comput_dom_thickness')
             read(value,*) comput_dom_thickness
          case ('verbose')
             read(value,*) verbose
          case ('DomDumpDir')
             write(DomDumpDir,'(a)') trim(value)
          case('decomp_dom_type')
             write(decomp_dom_type,'(a)') trim(value)
          case ('decomp_dom_ndomain')
             ndomain_present = .true.
             read(value,*) decomp_dom_ndomain
             allocate(decomp_dom_xc(decomp_dom_ndomain),decomp_dom_yc(decomp_dom_ndomain),decomp_dom_zc(decomp_dom_ndomain))
             allocate(decomp_dom_rin(decomp_dom_ndomain),decomp_dom_rout(decomp_dom_ndomain))
             allocate(decomp_dom_rsp(decomp_dom_ndomain),decomp_dom_size(decomp_dom_ndomain))
             allocate(decomp_dom_thickness(decomp_dom_ndomain))
          case ('decomp_dom_xc')
             read(value,*) decomp_dom_xc(:)
          case ('decomp_dom_yc')
             read(value,*) decomp_dom_yc(:)
          case ('decomp_dom_zc')
             read(value,*) decomp_dom_zc(:)
          case ('decomp_dom_rsp')
             read(value,*) decomp_dom_rsp(:)
          case ('decomp_dom_rin')
             read(value,*) decomp_dom_rin(:)
          case ('decomp_dom_rout')
             read(value,*) decomp_dom_rout(:)
          case ('decomp_dom_size')
             read(value,*) decomp_dom_size(:)
          case ('decomp_dom_thickness')
             read(value,*) decomp_dom_thickness(:)
          end select
       end do
    end if
    close(10)

    if (.not. ndomain_present) then ! assign default values
       allocate(decomp_dom_rsp(1))
       allocate(decomp_dom_xc(1),decomp_dom_yc(1),decomp_dom_zc(1))
       decomp_dom_rsp(1) = min(comput_dom_rsp + 0.05d0, 0.5d0)
       decomp_dom_xc(1) = 0.5d0
       decomp_dom_yc(1) = 0.5d0
       decomp_dom_zc(1) = 0.5d0
    end if

    ! make sure DomDumpDir ends with a /
    DomDumpDir = trim(DomDumpDir)//"/"
    
    call read_mesh_params(pfile)
    ! JB- 
    ! force data domain to cover full box for idealised models. 
    decomp_dom_type    = 'cube'
    decomp_dom_ndomain = 1
    allocate(decomp_dom_size(1))
    decomp_dom_size(1) = 1.0d0
    if (.not. allocated(decomp_dom_xc)) allocate(decomp_dom_xc(1),decomp_dom_yc(1),decomp_dom_zc(1))
    decomp_dom_xc(1) = 0.5d0     
    decomp_dom_yc(1) = 0.5d0
    decomp_dom_zc(1) = 0.5d0
    ! -JB
    call idealised_model_read_params(pfile)

    return

  end subroutine read_GenerateAMRmodel_params

  
  subroutine print_GenerateAMRmodel_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[GenerateAMRmodel]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  DomDumpDir      = ',trim(DomDumpDir)
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  comput_dom_type      = ',trim(comput_dom_type)
       write(unit,'(a,3(ES10.3,1x))')'  comput_dom_pos       = ',comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
       select case(comput_dom_type)
       case ('sphere')
          write(unit,'(a,ES10.3)')      '  comput_dom_rsp       = ',comput_dom_rsp
       case ('shell')
          write(unit,'(a,ES10.3)')      '  comput_dom_rin       = ',comput_dom_rin
          write(unit,'(a,ES10.3)')      '  comput_dom_rout      = ',comput_dom_rout
       case ('cube')
          write(unit,'(a,ES10.3)')      '  comput_dom_size      = ',comput_dom_size
       case ('slab')
          write(unit,'(a,ES10.3)')      '  comput_dom_thickness = ',comput_dom_thickness
       end select
       write(unit,'(a)')             '# domain decomposition parameters'
       write(unit,'(a,a)')           '  decomp_dom_type      = ',trim(decomp_dom_type)
       write(unit,'(a,i5)')          '  decomp_dom_ndomain   = ',decomp_dom_ndomain
       write(fmt,'(a,i5,a)') '(a,',decomp_dom_ndomain,'(ES10.3,1x))'
       write(unit,fmt)               '  decomp_dom_xc        = ',decomp_dom_xc(:)
       write(unit,fmt)               '  decomp_dom_yc        = ',decomp_dom_yc(:)
       write(unit,fmt)               '  decomp_dom_zc        = ',decomp_dom_zc(:)
       select case(decomp_dom_type)
       case ('sphere')
          write(unit,fmt)               '  decomp_dom_rsp       = ',decomp_dom_rsp(:)
       case ('shell')
          write(unit,fmt)               '  decomp_dom_rin       = ',decomp_dom_rin(:)
          write(unit,fmt)               '  decomp_dom_rout      = ',decomp_dom_rout(:)
       case ('cube')
          write(unit,fmt)               '  decomp_dom_size      = ',decomp_dom_size(:)
       case ('slab')
          write(unit,fmt)               '  decomp_dom_thickness = ',decomp_dom_thickness(:)
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
       call idealised_model_print_params(unit)
       else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')         '[GenerateAMRmodel]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  DomDumpDir = ',trim(DomDumpDir)
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  comput_dom_type      = ',trim(comput_dom_type)
       write(*,'(a,3(ES10.3,1x))')'  comput_dom_pos       = ',comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
       select case(comput_dom_type)
       case ('sphere')
          write(*,'(a,ES10.3)')      '  comput_dom_rsp       = ',comput_dom_rsp
       case ('shell')
          write(*,'(a,ES10.3)')      '  comput_dom_rin       = ',comput_dom_rin
          write(*,'(a,ES10.3)')      '  comput_dom_rout      = ',comput_dom_rout
       case ('cube')
          write(*,'(a,ES10.3)')      '  comput_dom_size      = ',comput_dom_size
       case ('slab')
          write(*,'(a,ES10.3)')      '  comput_dom_thickness = ',comput_dom_thickness
       end select
       write(*,'(a)')             '# domain decomposition parameters'
       write(*,'(a,a)')           '  decomp_dom_type      = ',trim(decomp_dom_type)
       write(*,'(a,i5)')          '  decomp_dom_ndomain   = ',decomp_dom_ndomain
       write(fmt,'(a,i5,a)') '(a,',decomp_dom_ndomain,'(ES10.3,1x))'
       write(*,fmt)               '  decomp_dom_xc        = ',decomp_dom_xc(:)
       write(*,fmt)               '  decomp_dom_yc        = ',decomp_dom_yc(:)
       write(*,fmt)               '  decomp_dom_zc        = ',decomp_dom_zc(:)
       select case(decomp_dom_type)
       case ('sphere')
          write(*,fmt)               '  decomp_dom_rsp       = ',decomp_dom_rsp(:)
       case ('shell')
          write(*,fmt)               '  decomp_dom_rin       = ',decomp_dom_rin(:)
          write(*,fmt)               '  decomp_dom_rout      = ',decomp_dom_rout(:)
       case ('cube')
          write(*,fmt)               '  decomp_dom_size      = ',decomp_dom_size(:)
       case ('slab')
          write(*,fmt)               '  decomp_dom_thickness = ',decomp_dom_thickness(:)
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_mesh_params
       write(*,'(a)')             ' '
       call idealised_model_print_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_GenerateAMRmodel_params

  
end program GenerateAMRmodel
