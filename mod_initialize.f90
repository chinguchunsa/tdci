module initialize
  

  use global_variables
  implicit none


  !: contains subroutine read_input
  !: contains subroutine set_default
  !: contains subroutine read_tdcihead
  !: contains subroutine set_variables
  !: contains subroutine allocate_main( option )
  !: contains subroutine deallocate_main( option )
  !: contains subroutine read_hamdata
  !: contains subroutine write_input( option )
  !: contains subroutine writeme_modreadin( myroutine, option )
  !: contains subroutine erorrs_modreadin( myroutine, option )
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_INPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_input  
    

    !: reads in file 'input'    
    !: sets default parameters if file 'input' is not found
    implicit none

    
    integer(8), parameter :: myfile=100
    
    logical :: iamhere
    integer(8) :: i, mystat
    

    !: /DYNAMICS/         init_coeffs(0:10), init_states(0:10).  to change, mod_variables_global.f90
    !: /FIELD_strengths/  read_emax(20).  to change, mod_variables_global.f90
    !: /FIELD_directions/ read_*(100).  to change, mod_variables_global.f90
    namelist /DYNAMICS/         init_coeffs, init_states, restart      
    namelist /FIELD/            dirform, ellipt, envelope, ncyc, omega, phase, euler
    namelist /FIELD_units/      omega_units !:, angle_units
    namelist /FIELD_strengths/  nemax, read_emax                                                 
    namelist /FIELD_directions/ ndir, read_theta, read_phi, read_x, read_y, read_z 
    namelist /SYSTEM/           dt, eigmax, heuristic, ionization, jobtype, nstep, outstep
    namelist /SYSTEM_units/     dt_units, eigmax_units
    namelist /InOutFILES/       tdcidatfile, outputfile, restartbinfile
    namelist /DAVIDSON/         flag_davidson
    

    !: set default optional parameters
    call set_default

    
    !: in input file does not exist, set default parameters
    inquire( file=trim(inputfile), exist=iamhere )        
    if ( .not.iamhere ) call write_input(myfile)
    

    !: read input.  see namelists.  not all WARNINGs are fatal
    open( unit=10, file=trim(inputfile) )

    read( 10, nml=FIELD, iostat=mystat,err=40 ) 
40  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD, iostat= ',i0)") mystat 
    rewind(10)

    read( 10, nml=FIELD_units, iostat=mystat,err=41 )
41  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_units, iostat= ',i0)") mystat 
    rewind(10)

    read( 10, nml=FIELD_strengths, iostat=mystat, err=47 )
47  if ( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_strengths, iostat= ',i0)") mystat
    rewind(10)

    read( 10, nml=FIELD_directions, iostat=mystat, err=48 )
48  if ( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=FIELD_directions, iostat= ',i0)") mystat
    rewind(10)    
    
    read( 10, nml=SYSTEM, iostat=mystat, err=42 ) 
42  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=SYSTEM, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10, nml=SYSTEM_units, iostat=mystat, err=43 )
43  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=SYSTEM_units, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10,nml=DYNAMICS,iostat=mystat,err=44 )
44  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=DYNAMICS, iostat= ',i0)") mystat 
    rewind(10)
    
    read( 10, nml=InOutFILES, iostat=mystat, err=45 )
45  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=InOutFILES, iostat= ',i0)") mystat
    rewind(10)

    read( 10, nml=DAVIDSON, iostat=mystat, err=46 )
46  if( mystat.ne.0 ) write(100,"(' WARNING WARNING WARNING WARNING reading nml=DAVIDSON, iostat= ',i0)") mystat 
    continue
    
    close(10)


    !: officially opening outputfile.  
    open( iout,file=trim(outputfile) )


    call writeme_modreadin( 'read_input', 'greeting' ) 
    call writeme_modreadin( 'read_input', 'date' )
    
    if ( .not.iamhere ) write(iout,'(A)') " WARNING: Could not find file 'input'. Using default input settings"
    write(iout,'(A)') " finished reading Field and System variables from file 'input' "
    call write_input(iout)
    
    
    call write_header( 'read_input', 'initialize', 'leave' )
    

  end subroutine read_input
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SET_DEAFULTINPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine set_default
    

    !: sets default values for variables read in from file 'input'
    implicit none
    

    !: default for /DYNAMICS/
    !: default init staate = ground state
    !: # of superposition states = init_stats(0) = 1
    init_coeffs    = dcmplx( 0.d0,0.d0 )
    init_coeffs(1) = dcmplx( 1.d0,0.d0 )
    init_states(0) = 1  
    init_states(1) = 1
    restart = .False.

    !: default for /FIELD/
    dirform  = 'polar'
    ellipt   = 1.d0
    envelope = 'trap'
    ncyc     = 7
    omega    = 0.057d0
    phase    = 90.d0
    euler    = 0.d0
    
    !: default for /FIELD_units/
    omega_units = 'au'
    !angle_units = 'deg'
    
    !: default for /FIELD_strengths/
    nemax  = 3 
    read_emax = -1000.d0  !: will assign later with MO energies if not read-in
    
    !: default for /FIELD_directions/   
    ndir     = 62
    read_theta = -1000.d0
    read_phi   = -1000.d0
    
    !: default for /InOutFiles/
    tdcidatfile     = 'TDCI.dat'
    outputfile      = 'OUTPUT'
    restartbinfile  = 'none'

    !: default for /SYSTEM/
    dt         = 0.05d0
    eigmax     = 10.d0
    ionization = -1.0 
    jobtype    = flag_cis
    nstep      = 16000
    outstep    = 50

    !: default for /SYSTEM_units/
    dt_units     = 'au'
    eigmax_units = 'au'    
    
    !: default for davidson diagonalization
    flag_davidson = .False.

    
  end subroutine set_default
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_TDCIHEAD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_tdcihead
    

    !: read in 'TDCI.dat'.  Exits if file does not exist
    implicit none
    
    logical        :: iamhere
    character(100) :: cskip
  
  
    call write_header( 'read_tdcihead','initialize','enter' )
    
    
    !: stop if TDCI.dat file does not exist
    inquire( file=trim(tdcidatfile), exist=iamhere )
    if ( .not.iamhere ) call errors_modreadin( 'read_tdci1', 'no_tdcidat' )    
    
    !: read TDCI.dat
    open( unit=10, file=trim(tdcidatfile) )    
    read(10,*)     nbasis, nrorb, noa, nva, nob, nvb, unrstedflag, vabsflag
    read(10,'(A)') job_title
    read(10,*)     charge, mult, natoms    
    
    !: allocate coordinate arrays
    allocate( xcoord(natoms), &
         ycoord(natoms), &
         zcoord(natoms), &
         myatom(natoms) )    
    read(10,*) ( (myatom(iatom), &
         xcoord(iatom), &
         ycoord(iatom), &
         zcoord(iatom)) , iatom=1, natoms )
    read(10,*) cskip, cskip, dipx00, dipy00, dipz00
    
    if( restart ) then
       close(10) 
       write(iout,'(A)') " finished reading heading in '"//trim(tdcidatfile)//"'"
    else
       write(iout,'(A)') " finished reading heading in '"//trim(tdcidatfile)//"'"
       write(iout,'(A)') ' '//trim(tdcidatfile)//" still open for reading "
    end if


    call write_header( 'read_tdcihead','initialize','leave' )
    
    
  end subroutine read_tdcihead
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE SET_VARIABLES
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine set_variables


    !: set global variables.  convert units to atomic units
    implicit none
    
    integer(8) :: itheta, iphi, idir, nphi, ntheta
    real(8) :: norm, dumtheta, dumphi, dtheta, dphi
    logical :: iamhere, gen_direction
    
    
    call write_header( 'set_variables','initialize','enter' )
    
    
    !: check units, convert to AU
    if ( trim(omega_units).ne.'au' ) then
       if ( trim(omega_units).ne.'nm' ) call errors_modreadin( 'read_input', 'omega_error' )
       omega = convert_freq(omega)
    end if
    
    if ( trim(dt_units).ne.'au' ) then
       dt = convert_time( dt, trim(dt_units), 'au' )
       if ( dt.eq.-1.d20 ) call errors_modreadin( 'read_input', 'dt_error' )
    end if
    
    !if ( trim(angle_units).ne.'rad' ) then
    !   if ( trim(angle_units).ne.'deg' ) call errors_modreadin( 'read_input', 'euler_error' )
    !   euler = convert_angles( euler )
    !   phase = convert_angles( phase )
    !end if
    
    if ( trim(eigmax_units).ne.'au' ) then
       eigmax = convert_energy( eigmax, trim(eigmax_units), 'au' )
       if ( eigmax.eq.0.d0 ) call errors_modreadin( 'read_input', 'eigmax_error' )
    end if
    
    
    !: set up flag for linear or circular pulse
    linear = .True.
    if ( envelope.eq.'cirl' .or. envelope.eq.'cirr' ) linear = .False.
    
    !: if unreasonable envelope, set to default
    select case ( trim(envelope) ) 
    case( 'none' ) ; go to 150
    case( 'cos2' ) ; go to 150
    case( 'trap' ) ; go to 150
    case( 'gaus' ) ; go to 150
    case( 'stat' ) ; go to 150
    case( 'band' ) ; go to 150
    case( 'sin2' ) ; go to 150
    case( 'cirl' ) ; go to 150
    case( 'cirr' ) ; go to 150
    case default 
       envelope = 'trap'
       linear   = .true.
       call errors_modreadin( 'set_variables', 'no_envelope' )
    end select
    
150 continue

    
    !: if unreasonable eigmax, set eigmax to default
    if(eigmax.eq.0.d0) then
       eigmax = 10.d0
       call errors_modreadin( 'read_input', 'eigmax_error' )
    end if
    

    !: if unreasonable nstep, set nstep to 1.5 * number of cycles + 2
    period = int( 2.d0 * pi / omega / dt )    !: total number of timestep per period
    if ( nstep.le.1 ) then       
       nstep  = int( 1.5d0 * dble(ncyc) * dble(period) + 2.d0 ) !: total number of propagation steps
       call errors_modreadin( 'set_variables', 'nstep_error' )
    end if
    nstep = nstep + 1 
    field_duration = dble(ncyc) * period + 1
    

    if ( ionization.lt.0.d0 ) write(iout,"(' ionization < 0, setting energy cutoff automatically')")
    if ( ionization.gt.0.d0 ) write(iout,"(' ionization > 0, energy cutoff set to ', f10.4, 'au')") ionization
    
    write(iout,'(A)') ''
    if ( trim(dirform).eq.'polar' ) write(iout,"(' field/light polarization/prop directions in polar coordinates')")
    if ( trim(dirform).eq.'cart'  ) write(iout,"(' field/light polarization/prop directions in Cartesian coordinates')") 

    
    !: set default directions if not specified     
    gen_direction = .False.
    do idir=1, ndir 
       if ( read_theta(idir).lt.-100.0 .or. read_phi(idir).lt.-100.0 ) gen_direction = .True.
    end do

    if ( gen_direction ) then
       dirform = 'polar'
       read_theta = 0.d0
       read_phi   = 0.d0
       dtheta = 30.d0     ;    ntheta = int(180.d0/dtheta)
       dphi   = 30.d0     ;    nphi   = int(360.d0/dphi) 
       write(iout,*) ntheta, dphi
       write(iout,'(A)') ' WARNING:  setting ndir = 62 for automatic field generation' ; write(iout,'(A)') ''
       write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") 1, read_theta(1), 1, read_phi(1)
       idir = 1
       do itheta=2, ntheta
          dumtheta = dble(itheta-1)*dtheta
          do iphi=1, nphi
             idir = idir + 1
             read_theta(idir) = dumtheta
             read_phi(idir)   = dble(iphi-1)*dphi
             write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") idir, read_theta(idir), idir, read_phi(idir)
          end do
       end do
       idir = idir + 1 
       read_theta(idir) = 180.d0
       read_phi(idir)   = 0.d0
       write(iout,"(' setting read_theta(',i2,')= ',f7.3,5x,'read_phi(',i2,')= ',f7.3 )") idir, read_theta(idir), idir, read_phi(idir)
    end if
    
    
    if( read_emax(1).lt.-100.0 ) nemax = 3 

    !: outsteps
    ndata = int(nstep/outstep) 
    write(iout,'(A)') '' ; write(iout,"(' write norm and wavefunction out every ',i0,' steps, total of ' i0,' steps')") outstep, ndata
    
    
    !: DATA FROM TDCI.DAT    

    if( trim(jobtype) .eq. flag_ip ) Qalloc_indices = .True.
    if( vabsflag .eq. 0 ) Qalloc_vabs = .False.
    if( trim(jobtype) .eq. flag_soc ) then
       Qalloc_Zcomplex = .True.
       Qalloc_socmo    = .True.
       Qread_socx_ao   = .True.
       Qread_socy_ao   = .True.
       Qread_socz_ao   = .True.
    end if
    

    !: restart, set job route
    if( restart ) then
       if ( trim(restartbinfile) .eq. 'none' ) call errors_modreadin( 'set_variables','norestart' )
       Qallocate_main = .True.
       Qread_hamdata  = .False.
       Qget_field     = .True.
       Qget_ham0      = .False.
       Qget_nstuse    = .False.
       Qalloc_vabsmo  = .False.
       Qalloc_dipmo   = .False.       
       Qget_1eham     = .False.
       Qget_expVabs   = .False.
       Qread_binaries = .True.
       Qsave          = .False.
       Qdealloc       = .False.
       Qpropagate     = .True.
    end if

    !: davidson flag
    if ( flag_davidson ) then
       Qallocate_main = .False.
       Qread_hamdata  = .True.
       Qget_field     = .False.
       Qwrite_specifics1 = .False.
       Qwrite_ham0    = .False.
       Qwrite_specifics2 = .False.
       Qget_ham0      = .False.
       Qget_nstuse    = .False.
       Qget_1eham     = .False.
       Qget_expVabs   = .False.
       Qread_binaries = .False.
       Qsave          = .False.
       Qdealloc       = .False.
       Qpropagate     = .False.
       !: alloc variables
       Qalloc_vabsmo  = .False.
       Qalloc_dipmo   = .False.
       Qalloc_socmo   = .False.
       Qread_vabs_ao  = .False.
       Qread_dipx_ao  = .True.
       Qread_dipy_ao  = .True.
       Qread_dipz_ao  = .True.
       Qread_socx_ao  = .False.
       Qread_socy_ao  = .False.
       Qread_socz_ao  = .False.
       Qread_orben    = .True.
       Qread_cmo      = .True.
    end if
    call writeme_modreadin( 'set_variables','routine' )

    
    !: parameters for unrestricted or restricted
    select case( trim(jobtype) ) 
    case( flag_cis )
       noanva = noa*nva
       nobnvb = nob*nvb
       if(unrstedflag.ne.0) then
          unrestricted  = .True.
          norb          = noa + nva + nob + nvb
          nstates       = noanva + nobnvb + 1
       else
          unrestricted = .False.
          nstates      = noanva + 1
          norb         = nrorb
       end if
    case( flag_soc )
       unrestricted = .True.
       noanva = noa*nva
       nobnvb = nob*nvb
       noanvb = noa*nvb
       nobnva = nob*nva
       nstates = noanva + nobnvb + noanvb + nobnva + 1 
       norb = noa + nva + nob + nvb
    case( flag_tda )
       noanva = noa*nva
       nobnvb = nob*nvb
       if(unrstedflag.ne.0) then
          unrestricted  = .True.
          norb          = noa + nva + nob + nvb
          nstates       = noanva + nobnvb + 1
       else
          unrestricted = .False.
          nstates      = noanva + 1
          norb         = nrorb
       end if
    case ( flag_ip ) 
       unrestricted = .True.
       if ( .not. IP_alpha_beta ) then 
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = nob*(nob-1)/2*nvb + nob*noanva + nob
       else
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = nob*(nob-1)/2*nvb + nob*noanva + nob +  &
               noa*(noa-1)/2*nva + noa*nobnvb + noa           
       end if
    case( flag_cisd ) 
       if(unrstedflag.ne.0) then
          unrestricted = .True.
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates = noanva*nobnvb + (noa)*(noa-1)*(nva)*(nva-1)/4 + (nob)*(nob-1)*(nvb)*(nvb-1)/4 + &
               noanva + nobnvb + 1 !: for singles and ground state
       else
          unrestricted = .False.
          noanva       = noa*nva
          nobnvb       = nob*nvb
          norb         = noa + nva + nob + nvb
          nstates      = (noa)*(noa-1)*(nva)*(nva-1)/4 + noanva + 1
       end if
    end select
    
    
    !: initialize to CIS parameters
    nstuse  = nstates


    !: keep track of mem
    tot_use_mem = 0

    
    call write_header( 'set_variables','initialize','leave' )
    
    
  end subroutine set_variables
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ALLOCATE_MAIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine allocate_main( option )

    
    !: allocate main arrays
    implicit none
    
    character(*), intent(in) :: option
    integer(8), parameter :: npart = 2 ! nparticle
    

100 format(' allocated rank 1 array',a15,' of length = ',i11)


    call write_header( 'allocate_main','initialize','enter' )
    
    
    select case( trim(option) )
    case( 'main' )
       
       !: arrays for field 
       if( Qalloc_field ) then
          allocate( fvect1(nstep) )        ; fvect1 = 0.d0
          allocate( env(nstep) )           ; env = 0.d0
          write(iout,100) 'fvect1', nstep  ; call track_mem( nstep )
          write(iout,100) 'env',    nstep  ; call track_mem( nstep )
          !: circular stuff
          if ( .not.linear ) then
             allocate( fvect2(nstep) )        ; fvect2 = 0.d0
             write(iout,100) 'fvect2', nstep  ; call track_mem(nstep)
          end if
          allocate( tdciresults(nemax*ndir) ) 
          write(iout,"(' allocated tdciresults ')")
       end if

       !: for eigen-vectors and eigen-values
       select case( trim(jobtype) )
       case( flag_cis ) 
          allocate( cis_vec(nstates*nstates) )  ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )          ; cis_eig = 0.d0
       case( flag_soc ) 
          allocate( Zcis_vec(nstates*nstates) ) ; Zcis_vec = dcmplx(0.d0,0.d0)
          allocate( cis_eig(nstates) )          ; cis_eig = 0.d0 
       case( flag_ip ) 
          allocate( cis_vec(nstates*nstates) )  ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )          ; cis_eig = 0.d0
       case ( flag_cisd )
          allocate( cis_vec(nstates*nstates) )  ; cis_vec = 0.d0
          allocate( cis_eig(nstates) )          ; cis_eig = 0.d0
       end select

       if ( Qalloc_Zcomplex ) then
          write(iout,100) 'cis_vec', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates)
          write(iout,100) 'cis_eig', 2*nstates          ;  call track_mem( 2*nstates )
       else
          write(iout,100) 'cis_vec', nstates*nstates    ;  call track_mem( nstates*nstates)
          write(iout,100) 'cis_eig', nstates            ;  call track_mem( nstates )
       end if

       !: allocate arrays to keep track of indices .  optional
       if ( Qalloc_indices ) then
          select case( trim(jobtype) ) 
          case( flag_cis ) 
             allocate( hole_index(nstates,1) )      ; hole_index = 0
             allocate( part_index(nstates,1) )      ; part_index = 0
             write(iout,100) 'hole_index', nstates  ;  call track_mem( nstates )
             write(iout,100) 'part_index', nstates  ;  call track_mem( nstates )
          case( flag_soc ) 
             allocate( hole_index(nstates,1) )      ; hole_index = 0
             allocate( part_index(nstates,1) )      ; part_index = 0
             write(iout,100) 'hole_index', nstates  ;  call track_mem( nstates )
             write(iout,100) 'part_index', nstates  ;  call track_mem( nstates )
          case( flag_tda ) 
             allocate( hole_index(nstates,1) )      ; hole_index = 0
             allocate( part_index(nstates,1) )      ; part_index = 0
             write(iout,100) 'hole_index', nstates  ;  call track_mem( nstates )
             write(iout,100) 'part_index', nstates  ;  call track_mem( nstates )             
          case( flag_ip )
             allocate( hole_index(nstates,2) )        ; hole_index = 0
             allocate( part_index(nstates,1) )        ; part_index = 0
             write(iout,100) 'hole_index', 2*nstates  ;  call track_mem( 2*nstates )
             write(iout,100) 'part_index', nstates    ;  call track_mem( nstates )
          end select
       end if
              
       !: for Vabs elements 
       if( Qalloc_vabs ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Zabp(nstates*nstates) )         ; Zabp = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Zabp', 2*nstates*nstates ; call track_mem( 2*nstates*nstates )
          else 
             allocate( abp(nstates*nstates) )         ;  abp = 0.d0
             write(iout,100) 'abp',  nstates*nstates  ;  call track_mem( nstates*nstates )
          end if
       end if
       
       !: for [R][e^(-Vabs)][R]^T elements
       if( Qalloc_exp_vabs ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Zexp_abp(nstates*nstates) )          ; Zexp_abp = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Zexp_abp', 2*nstates*nstates  ; call track_mem( 2*nstates*nstates )
          else
             allocate( exp_abp(nstates*nstates) )        ; exp_abp = 0.d0
             write(iout,100) 'exp_abp', nstates*nstates  ; call track_mem( nstates*nstates )
          end if
       end if
       
       !: for x, y, z dipoles in CIS basis
       if( Qalloc_tdxyz ) then
          if ( Qalloc_Zcomplex ) then
             allocate( Ztdx(nstates*nstates) )  ;  Ztdx = dcmplx( 0.d0,0.d0 )
             allocate( Ztdy(nstates*nstates) )  ;  Ztdy = dcmplx( 0.d0,0.d0 )
             allocate( Ztdz(nstates*nstates) )  ;  Ztdz = dcmplx( 0.d0,0.d0 )
             write(iout,100) 'Ztdx', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
             write(iout,100) 'Ztdy', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
             write(iout,100) 'Ztdz', 2*nstates*nstates  ;  call track_mem( 2*nstates*nstates )
          else             
             allocate( tdx(nstates*nstates) )        ;  tdx = 0.d0
             allocate( tdy(nstates*nstates) )        ;  tdy = 0.d0
             allocate( tdz(nstates*nstates) )        ;  tdz = 0.d0
             write(iout,100) 'tdx', nstates*nstates  ;  call track_mem( nstates*nstates )
             write(iout,100) 'tdy', nstates*nstates  ;  call track_mem( nstates*nstates )
             write(iout,100) 'tdz', nstates*nstates  ;  call track_mem( nstates*nstates )
          end if
       end if
       
       
    case( 'mo' )
       
       !: Vabs MO elements
       if( Qalloc_vabsmo ) then
          allocate( vabsmoa(nrorb,nrorb) )            ;  vabsmoa = 0.d0
          write(iout,100) 'vabs_moa', nrorb*nrorb     ;  call track_mem( nrorb*nrorb )
          if( unrestricted ) then
             allocate( vabsmob(nrorb,nrorb) )         ;  vabsmob = 0.d0
             write(iout,100) 'vabs_mob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          end if
       end if

       !: dipx, dipy, dipz MO elements
       if( Qalloc_dipmo ) then
          allocate( dipxmoa(nrorb,nrorb) )        ;  dipxmoa = 0.d0
          allocate( dipymoa(nrorb,nrorb) )        ;  dipymoa = 0.d0
          allocate( dipzmoa(nrorb,nrorb) )        ;  dipzmoa = 0.d0
          write(iout,100) 'dipxmoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          write(iout,100) 'dipymoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          write(iout,100) 'dipzmoa', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          if( unrestricted ) then
             allocate( dipxmob(nrorb,nrorb) )        ;  dipxmob = 0.d0
             allocate( dipymob(nrorb,nrorb) )        ;  dipymob = 0.d0
             allocate( dipzmob(nrorb,nrorb) )        ;  dipzmob = 0.d0
             write(iout,100) 'dipxmob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
             write(iout,100) 'dipymob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
             write(iout,100) 'dipzmob', nrorb*nrorb  ;  call track_mem( nrorb*nrorb )
          end if
       end if

       if ( Qalloc_socmo ) then
          allocate( socmoAA(nrorb,nrorb) )  ;  socmoAA = 0.d0
          allocate( socmoBB(nrorb,nrorb) )  ;  socmoBB = 0.d0
          allocate( socmoAB(nrorb,nrorb) )  ;  socmoAB = 0.d0
          write(iout,100) 'socmoAA', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
          write(iout,100) 'socmoBB', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
          write(iout,100) 'socmoAB', nrorb*nrorb  ; call track_mem( nrorb*nrorb ) 
       end if

    end select
    
    
    call write_header( 'allocate_main','initialize','leave' )


  end subroutine allocate_main
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DEALLOCATE_MAIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine deallocate_main( option )
    
    
    use read_integrals
    
    implicit none

    character(*), intent(in) :: option
    

    call write_header( 'deallocate_main','initialize','enter' )

    
    select case( trim(option) )
    case( '2e_int' )

       if( allocated(dabcdBB) ) then
          deallocate( dabcdBB )
          write(iout,"( deallocated dabcdBB )" )  ;  call track_mem( -nvb3*nvb3 )
       end if
       if( allocated(dabcdAB) ) then
          deallocate( dabcdAB )
          write(iout,"(' deallocated dabcdAB ')")  ;  call track_mem( -nva2*nvb2 )
       end if
       if( allocated(dabcdAA) ) then
          deallocate( dabcdAA )  
          write(iout,"(' deallocated dabcdAA ')")  ;  call track_mem( -nva3*nva3 )
       end if
       if( allocated(diabcBB) ) then
          deallocate( diabcBB )
          write(iout,"(' deallocated diabcBB ')")  ;  call track_mem( -nobnvb*nvb3 )
       end if
       if( allocated(diabcBA) ) then
          deallocate( diabcBA )
          write(iout,"(' deallocated diabcBA ')")  ;  call track_mem( -nobnvb*nva2 )
       end if
       if( allocated(diabcAB) ) then
          deallocate( diabcAB )
          write(iout,"(' deallocated diabcAB ')")  ;  call track_mem( -noanva*nvb2 )
       end if
       if( allocated(diabcAA) ) then
          deallocate( diabcAA )
          write(iout,"(' deallocated diabcAA ')")  ;  call track_mem( -noanva*nva3 )
       end if
       if( allocated(dijkaBB) ) then
          deallocate( dijkaBB )
          write(iout,"(' deallocated dijkaBB ')")  ;  call track_mem( -nob3*nobnvb )
       end if
       if( allocated(dijkaBA) ) then
          deallocate( dijkaBA )
          write(iout,"(' deallocated dijkaBA ')")  ;  call track_mem( -nob2*noanva )
       end if
       if( allocated(dijkaAB) ) then
          deallocate( dijkaAB )
          write(iout,"(' deallocated dijkaAB ')")  ;  call track_mem( -noa2*nobnvb )
       end if
       if( allocated(dijkaAA) ) then
          deallocate( dijkaAA )
          write(iout,"(' deallocated dijkaAA ')")  ;  call track_mem( -noa3*noanva )
       end if
       if( allocated(diajbBB) ) then
          deallocate( diajbBB ) 
          write(iout,"(' deallocated diajbBB ')")  ;  call track_mem( -nobnvb*nobnvb )
       end if
       if ( allocated(dijklBB) ) then
          deallocate( dijklBB )
          write(iout,'(" deallocated dijklBB ")')  ;  call track_mem( -nob3*nob3 )
       end if
       if( allocated(diajbBA) ) then
          deallocate( diajbBA )
          write(iout,'(" deallocated diajbBA ")')  ;  call track_mem( -nva2*nob2 )
       end if
       if( allocated(diajbAB) ) then
          deallocate( diajbAB )
          write(iout,'(" deallocated diajbAB ")')  ;  call track_mem( -nvb*noa2 )
       end if
       if( allocated(dijklAB) ) then
          deallocate( dijklAB )
          write(iout,'(" deallocated dijklAB ")')  ;  call track_mem( -nob2*noa2 )
       end if
       if( allocated(diajbAA) ) then
          deallocate( diajbAA )
          write(iout,'(" deallocated diajbAA ")')  ;  call track_mem( -noanva*noanva )
       end if
       if( allocated(dijklAA) ) then
          deallocate( dijklAA )
          write(iout,'(" deallocated dijklAA ")')  ;  call track_mem( -noa3*noa3 )
       end if
       if( allocated(dijabBB) ) then
          deallocate( dijabBB )
          write(iout,'(" deallocated dijabBB ")')  ;  call track_mem( -nvb3*nob3 )
       end if
       if( allocated(dijabAB) ) then
          deallocate( dijabAB ) 
          write(iout,'(" deallocated dijabAB ")')  ;  call track_mem( -nobnvb*noanva )
       end if
       if( allocated(dijabAA) ) then
          deallocate( dijabAA )
          write(iout,'(" deallocated dijabAB ")')  ;  call track_mem( -nva3*noa3 )
       end if
       
    case( '1e_int' )
    
       ntt = nbasis*(nbasis+1)/2

       if( allocated(soczao)  ) then
          deallocate(socxao, socyao, soczao)
          write(iout,'(A)') ' deallocated soczao, socyao, socxao '  
          call track_mem( -3*ntt )
       end if
       if( allocated(cmo_b) ) then
          deallocate(cmo_b)
          write(iout,'(A)') ' deallocated cmo_b '  
          call track_mem( -nrorb*nbasis )
       end if
       if( allocated(cmo_a) ) then
          deallocate(cmo_a) 
          write(iout,'(A)') ' deallocated cmo_a '
          call track_mem( -nrorb*nbasis )
       end if
       if( allocated(orben) ) then
          deallocate(orben)
          write(iout,'(A)') ' deallocated orben '
          call track_mem( -norb )
       end if
       if( allocated(dipzao)) then
          deallocate(dipzao, dipyao, dipxao)
          write(iout,'(A)') ' deallocated dipzao, dipyao, dipxao '
          call track_mem( -3*ntt )
       end if
       if( allocated(vabsao)) then
          deallocate(vabsao)   
          write(iout,'(A)') ' deallocated vabsao '
          call track_mem( -ntt )
       end if
       if( allocated(dipzmob) ) then
          deallocate(dipzmob, dipymob, dipxmob) 
          write(iout,'(A)') ' deallocated dipzmob, dipymob, dipxmob '
          call track_mem( -3*nrorb*nrorb )
       end if
       if( allocated(dipzmoa) ) then
          deallocate(dipzmoa,dipymoa,dipxmoa) 
          write(iout,'(A)') ' deallocated dipzmoa, dipymoa, dipxmoa '
          call track_mem( -3*nrorb*nrorb )
       end if
       if( allocated(vabsmob) ) then
          deallocate(vabsmob)
          write(iout,'(A)') ' deallocated vabsmob '
          call track_mem( -nrorb*nrorb )
       end if
       if( allocated(vabsmoa) ) then
          deallocate(vabsmoa)
          write(iout,'(A)') ' deallocated vabsmoa '
          call track_mem( -nrorb*nrorb )
       end if
       if( allocated(socmoAB) ) then
          deallocate(socmoAB,socmoBB,socmoAA) 
          write(iout,'(A)') ' deallocated socmoAB, socmoBB, socmoAA'
          call track_mem( -3*nrorb*nrorb )
       end if
       if( allocated(env) ) then  
          deallocate(env)
          write(iout,'(A)') ' deallocated env'
       end if

    end select

    call write_header( 'deallocate_main','initialize','leave' )


  end subroutine deallocate_main
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE READ_HAMDATA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_hamdata
    

    !: read in rest of tdci.dat to build the Hamiltonian


    use read_integrals
    implicit none
    
    
    call write_header( 'read_hamdata','initialize','enter' )
    

    call cpu_time(start)    
        

    call read_vabs_ao( Qread_vabs_ao ) !: read CAP elements in AO basis
    call read_dipx_ao( Qread_dipx_ao ) !: read Dipole-X in AO basis        
    call read_dipy_ao( Qread_dipy_ao ) !: read Dipole-Y in AO basis        
    call read_dipz_ao( Qread_dipz_ao ) !: read Dipole-Z in AO basis     
    call read_socx_ao( Qread_socx_ao ) !: read del x terms in AO basis    
    call read_socy_ao( Qread_socy_ao ) !: read del y terms in AO basis    
    call read_socz_ao( Qread_socz_ao ) !: read del z terms in AO basis    
    call read_orben( Qread_orben )     !: read in MO orbital energies     
    call read_cmo( Qread_cmo )         !: read in mo-lcao coefficients 


    Qread_buckets(2) = .True.
    Qread_buckets(5) = .True.
    if( unrestricted ) Qread_buckets(10) = .True. 
    if( trim(jobtype).eq.flag_soc )  Qread_buckets(1:8)  = .True.
    if( trim(jobtype).eq.flag_ip )   Qread_buckets(1:14) = .True.
    if( trim(jobtype).eq.flag_cisd ) Qread_buckets(1:21) = .True.


    if( trim(jobtype).eq.flag_tda ) then
       Qread_buckets = .False.
       call read_custom_ham
    end if
    
    
    !: read in integrals
    if( Qread_buckets(2) )  call read_bucket_2  !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB
    if( Qread_buckets(5) )  call read_bucket_5  !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB       
    if( Qread_buckets(10) ) call read_bucket_10 !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB
    if( Qread_buckets(1) )  call read_bucket_1  !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB
    if( Qread_buckets(3) )  call read_bucket_3  !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB          
    if( Qread_buckets(7) )  call read_bucket_7  !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB       
    if( Qread_buckets(8) )  call read_bucket_8  !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB       
    if( Qread_buckets(4) )  call read_bucket_4  !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L,IJ.LE.KL  IJKL       
    if( Qread_buckets(6) )  call read_bucket_6  !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL       
    if( Qread_buckets(9) )  call read_bucket_9  !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL       
    if( Qread_buckets(11) ) call read_bucket_11 !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA       
    if( Qread_buckets(12) ) call read_bucket_12 !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA       
    if( Qread_buckets(13) ) call read_bucket_13 !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA       
    if( Qread_buckets(14) ) call read_bucket_14 !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA       
    if( Qread_buckets(15) ) call read_bucket_15 !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC       
    if( Qread_buckets(16) ) call read_bucket_16 !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC
    if( Qread_buckets(17) ) call read_bucket_17 !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC       
    if( Qread_buckets(18) ) call read_bucket_18 !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC       
    if( Qread_buckets(19) ) call read_bucket_19 !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD       
    if( Qread_buckets(20) ) call read_bucket_20 !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD       
    if( Qread_buckets(21) ) call read_bucket_21 !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD
    
    
    call cpu_time(finish)    
    
    
    close(10)
    write(iout,"(' TDCI.dat read time:',f12.4,' seconds')") finish-start
    write(iout,'(A)') ' closed file '//trim(tdcidatfile)

    
    call write_header( 'read_hamdata','initialize','leave' )
    
    
  end subroutine read_hamdata
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITE_INPUT
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine write_input(myoption)


    implicit none
    integer(8), intent(in) :: myoption

    integer(8) :: i
    

    if( myoption.ne.iout ) open(unit=myoption, file='input')

    write(myoption,'(A)') ' &DYNAMICS'
    write(myoption,"(' init_states(0) =',i2)") init_states(0)
    do i=1, init_states(0) 
       write(myoption,"(' init_states(', i1 , ') =',i2)") i, init_states(i)
    end do
    do i=1, init_states(0)
       write(myoption,"( ' init_coeffs(', i1, ') = ','(',f7.5,',',f7.5,')' )") i, real(init_coeffs(i)), aimag(init_coeffs(i)) 
    end do
    if ( restart )      write(myoption,'(A)') ' restart        = .True.  '
    if ( .not.restart ) write(myoption,'(A)') ' restart        = .False. '    
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') ' &FIELD_units'
    write(myoption,'(A)') " omega_units = "//"'"//trim(omega_units)//"'"
    !write(myoption,'(A)') " angle_units = "//"'"//trim(angle_units)//"'"
    write(myoption,'(A)') ' /'
    
    write(myoption,'(A)') ' &FIELD'
    write(myoption,"(' dirform  =', a8  )") adjustr("'"//trim(dirform)//"'")
    write(myoption,"(' ellipt   =', f7.3)") ellipt
    write(myoption,"(' envelope =', a7  )") adjustr("'"//trim(envelope)//"'")
    write(myoption,"(' ncyc     =', i7  )") ncyc
    write(myoption,"(' omega    =', f7.3)") omega
    write(myoption,"(' phase    =', f7.3)") phase
    write(myoption,"(' euler    =', f7.3)") euler
    write(myoption,'(A)') ' /'
    
    write(myoption,'(A)') '&FIELD_strengths '
    write(myoption,"(' nemax = ', i7)") nemax
    do i=1, nemax
       write(myoption,"(' read_emax(',i1,') = ',f10.4 )") i, read_emax(i)
    end do
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') '&FIELD_directions '
    write(myoption,"(' ndir = ', i7 )") ndir
    if( trim(dirform).eq.'cart' ) then
       do i=1, ndir
          write(myoption,"(' read_x(',i2,',') = ',f8.4 )") i, read_x(i)
          write(myoption,"(' read_y(',i2,',') = ',f8.4 )") i, read_y(i)
          write(myoption,"(' read_z(',i2,',') = ',f8.4 )") i, read_z(i)
       end do
    else 
       do i=1, ndir
          write(myoption,"(' read_theta(',i2,') = ',f10.4)") i, read_theta(i)
       end do
       do i=1, ndir
          write(myoption,"(' read_phi(',i2,')   = ',f10.4)") i, read_phi(i)
       end do
    end if
    
    write(myoption,'(A)') ' &SYSTEM_units'
    write(myoption,'(A)') " dt_units     = 'au' "
    write(myoption,'(A)') " eigmax_units = 'au' "
    write(myoption,'(A)') ' /'

    write(myoption,'(A)') ' &SYSTEM'
    write(myoption,"(' dt          =', f7.3)") dt
    write(myoption,"(' eigmax      =', f7.3)") eigmax
    write(myoption,"(' ionization  =', f7.3)") ionization
    write(myoption,"(' jobtype     =', a7 )") adjustr("'"//trim(jobtype)//"'")
    write(myoption,"(' nstep       =', i7  )") nstep
    write(myoption,"(' outstep     =', i7  )") outstep
    write(myoption,'(A)') ' /'        

    write(myoption,'(A)') ' &InOutFILES'
    write(myoption,'(A)') ' tdcidatfile     = '//"'"//trim(tdcidatfile)//"'"
    write(myoption,'(A)') ' outputfile      = '//"'"//trim(outputfile)//"'"
    write(myoption,'(A)') ' restartbinfile  = '//"'"//trim(restartbinfile)//"'"
    write(myoption,'(A)') ' /'    

    write(myoption,'(A)') '&Davidson'
    if( flag_davidson ) write(myoption,'(A)') " flag_davidson = .True. "
    if ( .not. flag_davidson ) write(myoption,'(A)') " flag_davidson = .False. "    
    write(myoption,'(A)') ' /'

    if( myoption.ne.iout ) close(myoption)
    
    
  end subroutine write_input
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_MODREADIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_modreadin(myroutine,option)
    
    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option

    
    select case(trim(myroutine))
    !: subroutine read_input
    case('read_input')
       select case( trim(option) ) 
       case( 'greeting' ) 
          write(iout,'(A)') ' THIS IS A WONDERFUL DAY, I HAVE NEVER SEEN THIS ONE BEFORE'
          write(iout,'(A)') ' - Maya Angelous'
          write(iout,'(A)') divide
       case( 'date' )          
          write(iout,'(A)') ' I WAS COMPILED ON Wed Nov 27 13:06:46 EST 2019 '
          write(iout,'(A)') ' I AM A REVISED CODE FOR CW PULSE GENERATION '
          write(iout,'(A)') ' RAMPING PARAMETER SET TO RAMP_STEP=16000, NOT NSTEP'
          call dnt(iout)          
          call write_header( 'read_input','initialize','enter' )
       end select

       !: subroutine set_variables
    case('set_variables')
       select case( trim(option) )
       case('routine') 
          write(iout,'(A)') ' job route'
          write(iout,"(' read_input            ',l1)") Qread_input
          write(iout,"(' read_tdcihead         ',l1)") Qread_tdcihead
          write(iout,"(' set_variables         ',l1)") Qset_variables
          write(iout,"(' allocate_main         ',l1)") Qallocate_main
          write(iout,"(' read_hamdata          ',l1)") Qread_hamdata
          write(iout,"(' get_field             ',l1)") Qget_field
          write(iout,"(' get_ham0              ',l1)") Qget_ham0
          write(iout,"(' get_nstuse            ',l1)") Qget_nstuse
          write(iout,"(' get_1eham             ',l1)") Qget_1eham
          write(iout,"(' get_expVabs           ',l1)") Qget_expVabs
          write(iout,"(' read_restart_binaries ',l1)") Qread_binaries
          write(iout,"(' save_ham              ',l1)") Qsave
          write(iout,"(' deallocate_main       ',l1)") Qdealloc
          write(iout,"(' propagate             ',l1)") Qpropagate
       end select
       
    end select
    
    flush(iout)
    
100 format(' allocated rank 1 array',a15,' of length = ',i11)

    
  end subroutine writeme_modreadin
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ERRORS_MODREADIN
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_modreadin(myroutine,option)

    implicit none
    
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option


    select case( trim(myroutine) )

    !: subroutine read_tdcihead
    case( 'read_tdci1' )
       select case( trim(option) )
       case('no_tdcidat')
          write(iout,'(A)') " ERROR: Cannot find setup file 'TDCI.dat' "
          go to 100
       end select

    !: subroutine set_variables
    case( 'read_input' )
       select case( trim(option) ) 
       case('omega_error') 
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         omega_units must be in 'au' or 'nm' "
          go to 100 
       case('dt_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         dt_units must be in 'au' or 'ps' or 'fs' or 'as' "
          go to 100
       case('euler_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         euler_units must be in 'rad' or 'deg' "
          go to 100
       case('phase_error')
          write(iout,'(A)') " ERROR:  wrong unit specification in file 'input' "
          write(iout,'(A)') "         euler_units must be in 'rad' or 'deg' "
          go to 100          
       case('eigmax_error')
          write(iout,'(A)') " WARNING:  OVERWRITING eigmax to 10.0 au "
          eigmax  = 10.d0
          go to 200
       end select

    !: subroutine set_variables
    case( 'set_variables' )
       select case( trim(option) )
       case('no_envelope')
          write(iout,'(A)') ' WARNING: envelope type not found'
          write(iout,'(A)') ' WARNING: overriding user-defined pulse shape to static pulse'
          go to 200 
       case('nstep_error') 
          write(iout,'(A)') ' WARNING: nstep is less than 1'
          write(iout,"(' WARNING: overriding user-defined total number of timesteps to: ',i0)") nstep
          go to 200
       case( 'norestart' )
          write(iout,'(A)') ' ERROR:  could not find restart binary file'
          go to 100
       end select
       
    end select
    
100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " DON'T BE DISMAYED BY GOODBYES - Richard Bach "
    stop
    
200 continue
    
  end subroutine errors_modreadin
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module initialize
