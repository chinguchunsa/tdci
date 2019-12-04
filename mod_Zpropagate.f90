module Zpropagate
  
  use global_variables
  use analyze_psi
  use util
  
  implicit none


contains
  !==================================================================!
  !==================================================================!
  subroutine Ztrotter_linear
    
    use omp_lib    
    implicit none
    
    !: shared variables
    integer(8) :: nstuse2 
    complex(8) :: exphel(nstuse), psi0(nstuse)
    real(8)    :: tdvals1(nstuse)  !: eigenvectors of hp1
    complex(8) :: hp1(nstuse*nstuse)
    complex(8), allocatable :: hp2(:), tdvals2(:)  
    
    !: temporary arrays to be deallocated
    real(8) :: norm0
    real(8),    allocatable :: pop0(:)
    complex(8), allocatable :: psi_det0(:)
    

    !: private variables
    integer(8) :: i, j, ii, jj, k, kk
    integer(8) :: itime, idir, iemax
    integer(8) :: ithread, idata
    
    !: field info
    real(8) :: dirx1, diry1, dirz1, emax1, efield1
    real(8) :: dirx2, diry2, dirz2, emax2, efield2
    real(8) :: temp, temp1, temp2

    !: reults and psi stuff
    real(8)    :: norm, rate, mux, muy, muz
    complex(8) :: psi_j, psi_k
    complex(8) :: psi(nstuse), psi1(nstuse) 
    complex(8), allocatable :: psi2(:)

    !: file stuff
    integer(8)   :: funit1, funit2
    character(4) :: dirstr, emaxstr
    character(100) :: cifile, datafile
    !: lapack stuff and misc
    integer(8) :: info1, info2, lscratch
    real(8)    :: start1, start2, finish1, finish2
    real(8)  :: rwork(3*nstuse-2)
    complex(8) :: scratch(nstuse*nstuse) 
    
    
    call write_header( 'Ztrotter_linear','propagate','enter' )    
    call writeme_propagate( 'trot_lin', 'equation' ) 
    

    !: nstuse2
    nstuse2 = nstuse * nstuse
    
    !: initialize psi0
    psi0 = dcmplx(0.d0,0.d0)
    do i=1, init_states(0)
       psi0( init_states(i) ) = init_coeffs(i)
    end do
    

    !: normalize
    call get_norm( norm0, nstuse, psi0 )
    psi0 = psi0 / norm0 


    !: write psi0
    call writeme_propagate( 'trot_lin', 'psi0' )    
    
    
    !: get initial population
    allocate( pop0(norb), psi_det0(nstates) )    
    
    norm0 = 1.d0
    

    call get_Zpsid( nstuse, nstates, Zcis_vec, norm0, psi0, psi_det0 )
    

    select case ( trim(jobtype) ) 
    case( flag_cis ) 
       if ( unrestricted ) then
          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
       else
          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
       end if
    case( flag_tda ) 
       if ( unrestricted ) then
          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
       else
          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
       end if
    case( flag_soc ) 
       call pop_soc( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb)
    case( flag_ip) 
       call pop_ip( hole_index, noa,nob,norb,nstates,nva,nvb, part_index(:,1), pop0, psi_det0 )
    end select
    call writeme_propagate( 'trot_lin','pop0', pop0, norb )
    
    deallocate( pop0, psi_det0 )

    
    !: exphel = exp(-iH*dt/2)
    do i=1, nstuse
       temp = -0.5d0 * cis_eig(i) * dt
       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
    end do
    

    !: Start loop over directions.  Counters need to be passed in as non-derived datatype

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
    !$OMP cifile, datafile, dirstr, emaxstr, finish1, finish2, funit1, funit2, info1, info2, ithread, lscratch, start1, start2,  &
    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, psi_j, psi_k, temp, temp1, temp2,            &
    !$OMP norm, mux, muy, muz, rate,                                &
    !$OMP hp1, hp2, psi, psi1, scratch, rwork, tdvals1, tdvals2 ),  &
    !$OMP SHARED( au2fs, dt, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
    !$OMP Zabp, Zcis_vec, Zexp_abp, exphel, fvect1, fvect2, psi0, tdciresults, Ztdx, Ztdy, Ztdz)
    
    !$OMP DO  

    dir_loop : do idir=1, ndir
       
       ithread = omp_get_thread_num()
       call cpu_time( start1 )
       

       !: get directions stored in TDCItdciresults
       dirx1 = tdciresults(idir)%x0  
       diry1 = tdciresults(idir)%y0  
       dirz1 = tdciresults(idir)%z0  
       

       !: get mu dot e-vector matrix elements in CIS basis.  diagonalize
       hp1(:) = dirx1*Ztdx(1:nstuse2) + diry1*Ztdy(1:nstuse2) + dirz1*Ztdz(1:nstuse2)
       
       !: diagonalize mu dot E
       info1 = 10  ;  lscratch = nstuse2
       !:   zheev(   'v'   ,  'u'  , N    , A ,  LDA ,   W   ,WORK   ,LWORK,   RWORK,INFO)
       call zheev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,rwork,info1)
       
       !: hp = W * exp(-Vabs dt/2)
       call Zmatmult( hp1, Zexp_abp(1:nstuse2), nstuse )
       
       call cpu_time(finish1)
       
       !: loop over intensities
       emax_loop : do iemax=1, nemax
          
          !: for writing out files later
          write( emaxstr, '(i0)' ) iemax
          write( dirstr, '(i0)' )  idir

          !: get emax
          emax1 = tdciresults(1+(iemax-1)*ndir)%fstrength0
          

          !: cifile binary
          funit1  = iemax*100 + idir
          cifile ='CI-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
          
          open( unit=funit1, file=trim(cifile), form='unformatted' )
          write(funit1) ndata, nstuse, nstates
          
          !:            time,efld1,efld2,E1x(t),E1y(t),E1z(t),E2x(t),E2y(t),E2z(t), norm2
          write(funit1) 0.d0, 0.d0, 0.d0, dirx1, diry1, dirz1, 0.d0, 0.d0, 0.d0, 1.d0
          write(funit1) real(psi0)
          write(funit1) aimag(psi0) 
          flush(funit1)          
          
          !: datafile
          funit2 = 1000*iemax + idir
          datafile = 'RESULTS-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.out'
          open( unit=funit2,file=trim(datafile) )
          write( funit2, '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'
          write( funit2, "( '#',  20(f16.10,1x) )" ) emax1, 0.d0, &
               tdciresults(idir+(iemax-1)*ndir)%theta0, tdciresults(idir+(iemax-1)*ndir)%phi0,  0.d0,0.d0,  0.d0,0.d0, &
               dirx1, diry1, dirz1,  0.d0,0.d0,0.d0,  0.d0,0.d0,0.d0
          write( funit2,"(a5,3(a10,1x),3(a10,1x),2(1x,a15),3(1x,a15) )" ) '#','time(fs)','field1','field2', 'x', 'y', 'z', &
               'norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
          
          !$OMP CRITICAL
          !: all thread execute the code, but only one at a time
          write( iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0)" ) idir, iemax, ithread
          flush( iout )
          !$OMP END CRITICAL
                    
          !: initialize psi
          psi = psi0
          

          !: begin Looping over time
          call cpu_time( start2 )
          timestep_loop : do itime=1, nstep-1
             
             !: modified midpoint
             efield1 = 0.5d0 * emax1 * ( fvect1(itime) + fvect1(itime+1) )
             
             
             !: exp(-iHel dt/2 ) * psi
             do i=1, nstuse
                psi(i) = exphel(i) * psi(i)
             end do

             
             !: W * exp(-Vabs dt/2) * psi
             psi1 = dcmplx(0.d0,0.d0)
             do j = 1, nstuse
                jj = ( j-1) * nstuse  
                psi_j = psi(j)
                do k=1, nstuse
                   psi1(k) = psi1(k) + hp1( jj+k ) * psi_j
                end do
             end do

             
             !: exp(-E(t+dt/2)*mu*dt) * psi
             do j = 1, nstuse
                temp = dt * efield1 * tdvals1(j)
                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
             end do
             
             
             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W * psi
             do j = 1, nstuse
                jj = nstuse * ( j-1 )
                psi_j = dcmplx( 0.d0, 0.d0 )
                do k=1, nstuse
                   psi_j  = psi_j + dconjg( hp1(jj+k) ) * psi1(k)
                end do
                psi(j) = exphel(j) * psi_j
             end do
             
             write(*,*) itime, psi(1)

             analysis : if ( mod(itime,outstep).eq.0 ) then                
                
                idata = int( itime/outstep)  
                
                call get_norm( norm,nstuse, psi )                
                call get_Zexpectation( nstuse, norm, psi, Zabp, rate) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdx, mux ) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdy, muy ) !: rate expectation value
                call get_Zexpectation( nstuse, norm, psi, Ztdz, muz ) !: rate expectation value
                rate = -2.d0 * rate * norm**2


                write( funit2,"( i5,f10.4,1x,2(f10.7,1x),3(f10.7,1x),2(1x,f15.10),3(1x,f15.10) )") &
                     idata, dble(itime)*dt*au2fs, efield1, 0.d0, &
                     efield1*dirx1, efield1*diry1, efield1*dirz1, &
                     norm**2, rate/au2fs, mux, muy, muz
                flush(funit2)
                

                write(funit1) dble(itime)*dt*au2fs, efield1, 0.d0, dirx1, diry1, dirz1, 0.d0, 0.d0, 0.d0, norm**2

                write(funit1) real( psi )
                write(funit1) aimag( psi )
                
                flush(funit1)
                
             end if analysis
             
          end do timestep_loop
          call cpu_time(finish2)

          close(funit1)
          close(funit2)
          
          !$OMP CRITICAL                    
          !: record data at last timestep
          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz
          
          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread, idir, iemax
          write(iout,"(12x,'dir = (',f8.5,',',f8.5,',',f8.5,')    emax = ',f8.5,' au')")           dirx1, diry1, dirz1, emax1
          write(iout,"(12x,'(iemax=1) TD diag and TDvec*exp_abp time: ',f12.4,' s')")              finish1 - start1
          write(iout,"(12x,'(iemax=1) LAPACK dysev TD diagonalization INFO=',i0)")                 info1         
          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2 - start2  
          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
          flush(iout)

          !$OMP END CRITICAL
          
       end do emax_loop
    end do dir_loop

    !$OMP END DO
    !$OMP END PARALLEL
    
    call write_header( 'Ztrotter_linear','propagate','leave' )

    
  end subroutine Ztrotter_linear


!  !: NEED TO FIX FOR *Z*
!
!  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!  ! SUBROUTINE TROTTER_CIRCULAR
!  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!  
!  subroutine trotter_circular
!
!
!    ! <C> propagation using circularly polarized lights.  Takes in fvect1 and fvect2
!    ! <C> C(t+dt) = exp(-iHel dt/2)exp(-Vabs dt/2) * W1exp(-iE1(t) dt/2)W1* W2*exp(-iE2(t+dt/2)mu dt/2)*W2
!    ! <C>           *W1exp(-iE1(t) dt/2)W1 * exp(-Vabs dt/2)exp(-iHel dt/2)*C(t)
!
!
!    use omp_lib
!    implicit none
!
!    !: shared variables
!    integer(8) :: nstuse2 
!    complex(8) :: exphel(nstuse), psi0(nstuse)
!    real(8) :: hp1(nstuse*nstuse), tdvals1(nstuse)
!    real(8) :: hp2(nstuse*nstuse), tdvals2(nstuse)
!    
!    !: temporary arrays to be deallocated
!    real(8)    :: norm0
!    real(8),    allocatable :: pop0(:)
!    complex(8), allocatable :: psi_det0(:)    
!    
!
!    !: private variables
!    integer(8) :: i, j, ii, jj, k, kk 
!    integer(8) :: itime, idir, iemax
!    integer(8) :: ithread, idata
!
!    !: field info
!    real(8) :: dirx1, diry1, dirz1, emax1, efield1
!    real(8) :: dirx2, diry2, dirz2, emax2, efield2
!    real(8) :: temp, temp1, temp2
!
!    !: reults and psi stuff
!    real(8)    :: norm, rate, mux, muy, muz
!    complex(8) :: psi_j, psi_k
!    complex(8) :: psi(nstuse), psi1(nstuse) 
!    complex(8) :: psi2(:)    
!
!    !: file stuff
!    integer(8)   :: funit1, funit2
!    character(4) :: dirstr, emaxstr
!    character(100) :: cifile, datafile
!    !: lapack stuff and misc
!    integer(8) :: info1, info2, lscratch
!    real(8)    :: start1, start2, finish1, finish2
!    real(8) :: scratch(nstuse*nstuse) 
!
!
!
!    call write_header( 'trotter_circular','propagate','enter' )    
!    call writeme_propagate( 'trot_cir', 'equation' ) 
!
!
!    !: nstuse
!    nstuse2 = nstuse*nstuse
!
!
!    !: initialize psi0
!    psi0 = dcmplx(0.d0,0.d0)
!    do i=1, init_states(0)
!       psi0( init_states(i) ) = init_coeffs(i)
!    end do
!    
!
!    !: normalize
!    call get_norm( norm0, nstuse, psi0 )
!    psi0 = psi0 / norm0 
!    
!    !: write psi0
!    call writeme_propagate( 'trot_lin', 'psi0' )    
!    
!    
!    !: get initial population
!    allocate( pop0(norb), psi_det0(nstates) )    
!    
!    norm0 = 1.d0
!    
!    call get_psid( nstuse, nstates, cis_vec, norm0, psi0, psi_det0 )
!    select case ( trim(jobtype) ) 
!    case( flag_cis ) 
!       if ( unrestricted ) then
!          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
!       else
!          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
!       end if
!    case( flag_tda ) 
!       if ( unrestricted ) then
!          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0, nob,nvb )
!       else
!          call pop_cis( hole_index(:,1), noa,norb,nstates,nva, part_index(:,1), pop0, psi_det0 )
!       end if
!    case( flag_ip) 
!       call pop_ip( hole_index, noa,nob,norb,nstates,nva,nvb, part_index(:,1), pop0, psi_det0 )
!    end select
!    call writeme_propagate( 'trot_lin','pop0', pop0, norb )
!    
!    deallocate( pop0, psi_det0 )
!
!
!    !: exphel = exp(-iH*dt/2)
!    do i=1, nstuse
!       temp = -0.5d0 * cis_eig(i) * dt
!       exphel(i) = dcmplx( dcos(temp) , dsin(temp) )
!    end do
!    
!
!    !: Start loop over directions
!
!    !$OMP PARALLEL DEFAULT(NONE),&
!    !$OMP PRIVATE( i, idata, idir, iemax, ii, itime, j, jj, k, kk, &
!    !$OMP cifile, datafile, dirstr, emaxstr, finish1, finish2, funit1, funit2, info1, info2, ithread, lscratch, start1, start2,  &
!    !$OMP dirx1, dirx2, diry1, diry2, dirz1, dirz2, efield1, efield2, emax1, emax2, psi_j, psi_k, temp, temp1, temp2,            &
!    !$OMP norm, mux, muy, muz, rate,                               &
!    !$OMP hp1, hp2, psi, psi1, scratch, tdvals1, tdvals2 ),        &
!    !$OMP SHARED( dt, au2fs, iout, ndata, ndir, nemax, nstates, nstep, nstuse, nstuse2, outstep, &
!    !$OMP abp, cis_vec, exp_abp, exphel, fvect1, fvect2, psi0, tdciresults, tdx, tdy, tdz)
!    
!    !$OMP DO
!    dir_loop : do idir=1, ndir 
!
!       ithread = omp_get_thread_num()
!       call cpu_time( start1 ) 
!
!
!       !: get directions stored in TDCItdciresults
!       dirx1 = tdciresults(idir)%x1 ; dirx2 = tdciresults(idir)%x2
!       diry1 = tdciresults(idir)%y1 ; diry2 = tdciresults(idir)%y2
!       dirz1 = tdciresults(idir)%z1 ; dirz2 = tdciresults(idir)%z2
!
!
!       !: Form the transition dipole in (dirx1,diry1,dirz1) and (dirx2,diry2,dirz2) directions,
!       hp1(:) = dirx1*tdx(1:nstuse2) + diry1*tdy(1:nstuse2) + dirz1*tdz(1:nstuse2)
!       hp2(:) = dirx2*tdx(1:nstuse2) + diry2*tdy(1:nstuse2) + dirz2*tdz(1:nstuse2)
!
!
!       !: diagonalize hp1 and hp2
!       info1=10  ;  info2=10  ;  lscratch = nstuse2
!       call dsyev('vectors','upper',nstuse,hp1,nstuse,tdvals1,scratch,lscratch,info1)
!       call dsyev('vectors','upper',nstuse,hp2,nstuse,tdvals2,scratch,lscratch,info2)
!
!       
!       !: hp2 = hp2 * hp1 ; hp2 = W2 * W1
!       call matmult( hp2,hp1,nstuse )
!       !: hp1 = W1 * exp(-Vabs dt/2)
!       call matmult( hp1,exp_abp(1:nstuse*nstuse),nstuse )
!
!       
!       call cpu_time( finish1 )
!
!
!       !: loop over intensities
!       emax_loop : do iemax=1, nemax
!
!
!          !: for writing out files later
!          write( emaxstr, '(i0)' ) iemax
!          write( dirstr, '(i0)' ) idir
!
!          !: get emax ; emax1==emax2 in this version
!          emax1 = tdciresults((iemax-1)*ndir+1)%fstrength1
!          emax2 = tdciresults((iemax-1)*ndir+1)%fstrength2
!          
!
!          !: cifile binary
!          funit1  = iemax*100 + idir
!          cifile ='CI-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.bin'
!          
!          open( unit=funit1, file=trim(cifile), form='unformatted' )
!          write(funit1) ndata, nstuse, nstates
!          
!          write(funit1) 0.d0, 0.d0, dirx1, diry1, dirz1, dirx2, diry2, dirz2, 1.d0
!          write(funit1) real(psi0)
!          write(funit1) aimag(psi0) 
!          flush(funit1)          
!          
!          !: datafile
!          funit2 = 1000*iemax + idir
!          datafile = 'RESULTS-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.out'
!          open( unit=funit2,file=trim(datafile) )
!          write( funit2, '(A)' ) '# emax1 emax2 theta0 phi0 theta1 phi1 theta2 phi2 dirx0 diry0 dirz0 dirx1 diry1 dirz1 dirx2 diry2 dirz2'
!          write( funit2, "( '#',  20(f16.10,1x) )" ) emax1, emax2, &
!               tdciresults(1+(iemax-1)*ndir)%theta0, tdciresults(1+(iemax-1)*ndir)%phi0, &
!               tdciresults(1+(iemax-1)*ndir)%theta1, tdciresults(1+(iemax-1)*ndir)%phi1, &
!               tdciresults(1+(iemax-1)*ndir)%theta2, tdciresults(1+(iemax-1)*ndir)%phi2, &
!               tdciresults(1+(iemax-1)*ndir)%x0, tdciresults(1+(iemax-1)*ndir)%y0, tdciresults(1+(iemax-1)*ndir)%z0, &               
!               dirx1, diry1, dirz1,  dirx2, diry2, dirz2
!          write( funit2,"(a5,3(a10,1x),2(1x,a15),3(1x,a15) )" ) '#','time(fs)','field1','field2','norm2','rate(fs-1)', 'mu_x(au)','mu_y(au)','mu_z(au)'
!
!
!          !$OMP CRITICAL
!          write(iout,"(' start propagation for direction',i4,' intensity',i4,' thread # ',i0 )") idir, iemax, ithread
!          flush(iout)
!          !$OMP END CRITICAL
!          
!
!          !: initialize psi
!          psi  = psi0
!
!          !: begin Looping over time
!          call cpu_time( start2 )
!          timestep_loop : do itime=1, nstep-1
!             
!
!             efield1 = 0.5d0 * emax1*(fvect1(itime)+fvect1(itime+1)) !: modified midpoint
!             efield2 = 0.5d0 * emax2*(fvect2(itime)+fvect2(itime+1)) !: modified midpoint
!             temp1 = dt * efield1 * 0.5d0 
!             temp2 = dt * efield2
!             
!             
!             !: exp(-iHel dt/2 ) * psi
!             do j = 1, nstuse
!                psi(j) = exphel(j) * psi(j)
!             end do
!             
!
!             !: W1 * exp(-Vabs dt/2) * psi
!             psi1 = dcmplx( 0.d0, 0.d0 )
!             do j = 1, nstuse
!                jj = (j-1) * nstuse
!                psi_j = psi(j)
!                do k = 1 , nstuse
!                   psi1(k) = psi1(k) + hp1(jj+k) * psi_j
!                end do
!             end do
!
!             
!             !: exp( iE1(t+dt/2)*mu * dt/2) * psi
!             do j = 1, nstuse
!                temp = temp1 * tdvals1(j)
!                psi1(j) = dcmplx( dcos(temp),dsin(temp) ) * psi1(j)
!             end do
!             
!             
!             !: W2*W1 * psi 
!             psi = dcmplx( 0.d0, 0.d0 )
!             do j = 1, nstuse
!                jj = (j-1)*nstuse
!                psi_j  = psi1(j)
!                do k=1, nstuse
!                   psi(k) = psi(k) + hp2(jj+k) * psi_j
!                end do
!             end do
!             
!             
!             !: exp( iE2(t+dt/2)*mu* dt ) * psi
!             do j = 1, nstuse
!                temp = temp2 * tdvals2(j)
!                psi(j) = dcmplx(dcos(temp),dsin(temp)) * psi(j)
!             end do
!
!
!             !: W2*W1 * psi
!             do j = 1, nstuse
!                jj = nstuse*(j-1)
!                psi_j = dcmplx( 0.d0, 0.d0 )
!                do k=1, nstuse 
!                   psi_j = psi_j + hp2(jj+k) * psi(k)
!                end do
!                psi1(j) = psi_j
!             end do
!             
!
!             !: exp( iE1(t+dt/2)*mu*dt/2) * psi
!             do j = 1, nstuse
!                temp = temp1 * tdvals1(j)
!                psi1(j) = dcmplx(dcos(temp),dsin(temp)) * psi1(j)
!             end do
!             
!
!             !: exp(-iHel dt/2) * exp(-Vabs dt/2) * W1 * psi
!             do j = 1, nstuse
!                jj = nstuse*(j-1)
!                psi_j = dcmplx( 0.d0, 0.d0 )
!                do k=1, nstuse
!                   psi_j = psi_j + hp1(jj+k) * psi1(k)
!                end do
!                psi(j) = exphel(j) * psi_j
!             end do
!             
!             
!
!             analysis : if ( mod(itime,outstep).eq.0 ) then
!
!                idata = int(itime/outstep)
!                call get_norm( norm,nstuse, psi )                
!                call get_expectation( nstuse, norm, psi, abp, rate) !: rate expectation value
!                call get_expectation( nstuse, norm, psi, tdx, mux ) !: rate expectation value
!                call get_expectation( nstuse, norm, psi, tdy, muy ) !: rate expectation value
!                call get_expectation( nstuse, norm, psi, tdz, muz ) !: rate expectation value
!                rate = -2.d0 * rate * norm**2
!
!                
!                write( funit2,"( i5,f10.4,1x,2(f10.7,1x),2(1x,f15.10),3(1x,f15.10) )") &
!                     idata, dble(itime)*dt*au2fs, efield1, efield2, norm**2, rate/au2fs, mux, muy, muz
!                flush(funit2)
!                                
!                write(funit1) dble(itime)*dt*au2fs, efield1, efield2, dirx1, diry1, dirz1, dirx2, diry2, dirz2, norm**2
!
!                write(funit1) real( psi )
!                write(funit1) aimag( psi )
!                
!                flush(funit1)
!                
!             end if analysis
!
!          end do timestep_loop
!          call cpu_time(finish2)
!
!          close(funit1)
!          close(funit2)
!
!
!          !$OMP CRITICAL
!          !: record data at last timestep
!          tdciresults( idir + (iemax-1)*ndir)%norm0 = norm**2
!          tdciresults( idir + (iemax-1)*ndir)%dipx  = mux
!          tdciresults( idir + (iemax-1)*ndir)%dipy  = muy
!          tdciresults( idir + (iemax-1)*ndir)%dipz  = muz
!
!          write(iout,"(' thread # ',i0,' propagation done for direction',i4,' and intensity',i4)") ithread,idir,iemax
!          write(iout,"(12x,'dir1 = (',f8.5,',',f8.5,',',f8.5,')  emax1 = ',f8.5,' au')") dirx1,diry1,dirz1,emax1
!          write(iout,"(12x,'dir2 = (',f8.5,',',f8.5,',',f8.5,')  emax2 = ',f8.5,' au')") dirx2,diry2,dirz2,emax2
!          write(iout,"(12x,'(iemax=1) TD diag and TDvec*exp_abp time: ',f12.4,' s')") finish1 - start1
!          write(iout,"(12x,'(iemax=1) LAPACK dysev TD diagonalization INFO=',i0)") info1
!          write(iout,"(12x,'propagation time:',f12.4,' s')") finish2-start2
!          write(iout,"(12x,'final norm = ',f10.5)")          norm**2
!          flush(iout)
!          !$OMP END CRITICAL
!
!
!       end do emax_loop
!    end do dir_loop
!
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    call write_header( 'trotter_circular','propagate','leave' )
!
!
!  end subroutine trotter_circular
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_PROPAGATE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_propagate(myroutine,option,pop,norb)

    implicit none

    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option

    integer(8), optional, intent(in) :: norb
    real(8),    optional, intent(in) :: pop(norb)
    
    integer(8) :: i, iorb
   
    
    select case( trim(myroutine) )
    case( 'trot_lin' ) 
       trotter_linear : select case( trim(option) ) 
       case('equation') 
          write(iout,'(A)') " C(t+dt) = exp[-iHel dt/2] * exp[-Vabs dt/2] * "
          write(iout,'(A)') "           W^T * exp[iE(t+dt/2)*mu dt/2] * W * "
          write(iout,'(A)') "           exp[-Vabs dt/2] * exp[-iHel dt/2] * C(t)"                 
          write(iout,"(' start propagation for ',i10,' timesteps')")        nstep
          write(iout,"('                   for ',i10,' field directions')") ndir
          write(iout,"('                   for ',i10,' field strengths')")  nemax    
       case('psi0')
          write(iout,'(A)') ' '
          write(iout,"(A)") " INITIALIZED STATE: "
          write(iout,"(A)",advance="no") '     |psi(0)> = '
          write(iout,100) ( '(',init_coeffs(i),')', init_states(i), i=1, init_states(0) )
       case('pop0')
          write(iout,'(A)') ' '
          write(iout,"(A)") " INITIALIZED MO POPULATION"
          if( unrestricted ) then
             write(iout,101) ( pop(iorb), iorb=1, noa )
             write(iout,300) ( pop(iorb), iorb=(noa+1), nrorb )
             write(iout,200) ( pop(iorb), iorb=(nrorb+1),(nrorb+nob) )
             write(iout,400) ( pop(iorb), iorb=(nrorb+nob+1),norb )
             write(iout,'(A)') ' '
          else
             write(iout,500) ( pop(iorb), iorb=1, noa )
             write(iout,600) ( pop(iorb), iorb=(noa+1),nrorb )
             write(iout,'(A)') ' '
          end if
       end select trotter_linear
    case( 'trot_circ' )
       trotter_circ : select case( option ) 
       case('equation')
          write(iout,'(A)') " C(t+dt) = exp[-iHel dt/2] * exp[-Vabs dt/2] * W1^T * exp[iE1(t+dt/2)*mu dt/2] * W1 *"
          write(iout,'(A)') "           W2^T * exp[iE2(t+dt/2)*mu dt] * W2 * "
          write(iout,'(A)') "           W1^T * exp[iE1(t+dt/2)*mu dt/2]  * W1 * exp[-Vabs dt/2] * exp[-iHel dt/2] * C(t)"
          
          write(iout,"(' start propagation for ',i10,' timesteps')")  nstep
          write(iout,"('                   for ',i10,' directions')") ndir
          write(iout,"('                   for ',i10,' strengths')")  nemax
       end select trotter_circ
    end select
    

    flush(iout)

100 format( 10(a1,f7.5,','f7.5,a1,'|',i0,'>  ') )
101 format( '  occ_a:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
300 format( '  vir_a:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
200 format( '  occ_b:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
400 format( '  vir_b:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
500 format( '    occ:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )
600 format( '    vir:', 14(1x,f7.5) / 100( 15(1x,f7.5) /) )


  end subroutine writeme_propagate
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module Zpropagate

