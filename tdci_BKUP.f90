!--------------------------------------------------------------------------------!
! FORTRAN code to propagate a molecule in an electric field                      !
! using Absorbing Boundary Potential to calculate ionization rates               !
!                                                                                !
! Initial version written by J.Sonk, Sept 2012                                   !
! Extensions/modifications by P.Krause, H.B.Schlegel                             !
! Restricted and unrestricted CIS, parallel implementation by H.B.Schlegel       !
!--------------------------------------------------------------------------------!
program main

  
  use global_variables !: inherits units, setup_variables, route_control
  use initialize
  use write_info
  use get_field
  use get_ham          !: inherits get_ham0
  use propagate
  use mod_davidson_ip  !: only to check Hamiltonian

  
  implicit none
  integer(8) :: i,j,a,b,ia, istate, iout2 !: temporary cisd variables
  real(8)  :: masterstart, masterfinish
  

  !: output opened in read_input with iout=42   
  call cpu_time(masterstart)
  
  
  !: read and set input parameters (module initialize) 
  if( Qread_input    ) call read_input  
  if( Qread_tdcihead ) call read_tdcihead
  if( Qset_variables ) call set_variables
  

  !: allocate main arrays before allocating temporary arrays (module initialize)
  if ( Qallocate_main ) then
     call allocate_main( 'main' )
     call allocate_main( 'mo' )
  end if
  
  
  !: read 1-electron and 2-electron matrix elements (module initialize)
  if( Qread_hamdata ) then
     call read_hamdata
     if ( Qwrite_mo_energies ) call write_mo_energies
  end if

  
  !: generate field 
  !: read polarization directions, dirform<0 in polar ; dirform>0, in Cartesian
  !: for circularly polarized pulses, two vectors perpendicular to propdirection   
  if( Qget_field ) then
     call check_emax
     call get_lindirection
     if (.not.linear) call get_circdirection
     call shape_field
     if ( Qwrite_fshape ) call write_field_shape
  end if


  !: write out useful information I
  if( Qwrite_specifics1 ) call write_specifics1
  
  
  stop


  !: form CIS Hamiltonian and diagonalize 
  !: matrix elements stored in cis_vec since cis_vec will be fed into dysev diagonalization
  !: dysev will spit out eigenvectors into cis_vec.
  if ( Qget_ham0 ) then
     select case ( trim(jobtype) ) 
     case( flag_cis ) 
        call get_cis
        cis_vec(1)=-500.d0
        call diagonalize
        cis_eig(1)=0.d0     
        call deallocate_main( '2e_int' )
     case( flag_tda ) 
        call get_cis_index
        cis_vec(1)=-500.d0
        call diagonalize
        cis_eig(1)=0.d0   
     case( flag_ip  ) 
        call get_ip_cisd
        call diagonalize
        call deallocate_main( '2e_int' )
     case( flag_cisd ) 
        call get_cisd_index
        
        !: tmp
        iout2=780
        open( unit=iout2, file='cisd_hamiltonian.out' )
        

        istate = 1         
        !: with ground
        call get_cisd_0( cis_eig )
        i=9
        !write( i, "('# -----',3(i0,1x))" ) istate, 0, 0
        !call write_cisd_hamiltonian( cis_eig, i )
        call write_cisd_hamiltonian( cis_eig, iout2 )
        !: singles alpha
        do i=1, noa
           do a=1, nva
              istate = istate + 1 
              !write( istate*10, "('# -----',3(i0,1x))" ) istate, -i, -a
              call get_cisd_ia_AA( i, a, cis_eig )
              !call write_cisd_hamiltonian( cis_eig, istate*10 )
              call write_cisd_hamiltonian( cis_eig, iout2 )
           end do
        end do
        !: singles beta
        do i=1, nob
           do a=1, nvb
              istate = istate + 1
              !write( istate*10, "('# -----',3(i0,1x))" ) istate, i, a
              call get_cisd_ia_BB( i, a, cis_eig )
              !call write_cisd_hamiltonian( cis_eig, istate*10 )
              call write_cisd_hamiltonian( cis_eig, iout2 )
           end do
        end do
        !: doubles alpha beta
        do i=1, noa
           do j=1, nob
              do a=1, nva
                 do b=1, nvb
                    istate = istate + 1 !cisd_indices(-a,b,-i,j)
                    !write( istate*10, "('# -----',5(i0,1x))" ) istate, -i, j, -a, b
                    call get_cisd_ijab_ABAB( i,j,a,b,cis_eig )
                    !call write_cisd_hamiltonian( cis_eig, istate*10 )
                    call write_cisd_hamiltonian( cis_eig, iout2 )
                 end do
              end do
           end do
        end do
        !: doubles alpha alpha
        do i=2, noa
           do j=1, (i-1)
              do a=2, nva
                 do b=1, (a-1)
                    istate = istate + 1 !cisd_indices(-a,-b,-i,-j)
                    !write( istate*10, "('# -----',5(i0,1x))" ) istate, -i, -j, -a, -b
                    call get_cisd_ijab_AAAA( i,j,a,b,cis_eig )
                    !call write_cisd_hamiltonian( cis_eig, istate*10 )
                    call write_cisd_hamiltonian( cis_eig, iout2 )
                 end do
              end do
           end do
        end do
        !: doubles beta beta
        do i=2, nob
           do j=1, (i-1)
              do a=2, nvb
                 do b=1, (a-1)
                    istate = istate + 1 !cisd_indices(a,b,i,j)
                    !write( istate*10, "('# -----',5(i0,1x))" ) istate, i, j, a, b
                    call get_cisd_ijab_BBBB( i,j,a,b,cis_eig )
                    !call write_cisd_hamiltonian( cis_eig, istate*10 )
                    call write_cisd_hamiltonian( cis_eig, iout2 )
                 end do
              end do
           end do
        end do
   
        close(iout2)

        stop
        call diagonalize         
     end select
     if( Qwrite_ham0 ) call write_ham0     
  end if

  
  !: davidson
  if ( flag_davidson ) then
     select case( trim(jobtype) ) 
     case( flag_ip )
        allocate( hole_index(nstates,2) ) ; hole_index = 0
        allocate( part_index(nstates,1) ) ; part_index = 0
        call get_ip_index
        call davidson_ip
     case( flag_cisd )
        write(iout,'(A)') " Getting there "
        stop
        !call get_cisd_index !: assigns state_index
        !call davidson_cisd
     end select
  end if
        

  !: set NstUse = number of states to use <= nstates
  if( Qget_nstuse ) call get_nstuse  
  

  !: Transform ONE electron integrals from AO to MO basis (VabsAO,Dip*AO-->VabsMO,Dip*MO)
  !: Transform from MO basis to CIS states (abp,TDX,TDY,TDZ)    
  !: Or read in the binaries if restart = .True.
  if( Qget_1eham ) then
     call get_1eham
     if( Qwrite_mo_elements ) call write_mo_elements
  end if
  if( Qread_binaries ) call read_restart_bin


  !: write out useful information II
  if( Qwrite_specifics2 ) call write_specifics2

  
  !: compute exp(-Vabs *dt/2) if restart = .False.
  if( Qget_expVabs ) call get_expVabs
  

  !: save for later restart.  Need to change Qsave in mod_variables_control to False if not.
  if( Qsave ) call save_restart_bin


  !: deallocate un-used arrays
  if( Qdealloc ) call deallocate_main( '1e_int' )

  
  !: PROPAGATE 
  if ( Qpropagate ) then
     if( linear ) then
        call trotter_linear
     else
        call trotter_circular
     end if
     call write_summary
  end if
  

  call cpu_time(masterfinish)  

  
  write(iout,'(A)') divide
  write(iout,'(A)') ' '
  write(iout,"(' ---> total elapsed time:',f12.4,' s')") masterfinish-masterstart
  write(iout,'(A)') ' '
  call dnt(iout)
  write(iout,'(A)') ' '
  write(iout,'(A)') " MY MISSION IN LIFE IS NOT MERELY TO SURVIVE, BUT TO THRIVE;"
  write(iout,'(A)') " AND TO DO SO WITH SOME PASSION,"
  write(iout,'(A)') " SOME COMPASSION, SOME HUMOR, AND SOME STYLE"
  write(iout,'(A)') ' -Maya Angelou'
  flush(iout)

  close(iout)
  

end program main
