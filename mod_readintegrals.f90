module read_integrals


  use global_variables
  implicit none

  
  integer(8) :: noa2, nva2, nob2, nvb2, noa3, nob3, nva3, nvb3, ntt
  real(8)    :: rskip
  character(1000) :: cskip
  character(1000), parameter :: form1e = "(1x,'allocated rank 1 array ',a8,' of size ',i11)"
  character(1000), parameter :: form2e = "(1x,'allocated rank 2 array ',a8,' of size ',i11,i11)"


contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_size_variables
    

    implicit none

    
    noa2 = noa*(noa+1)/2
    nva2 = nva*(nva+1)/2
    nob2 = nob*(nob+1)/2
    nvb2 = nvb*(nvb+1)/2
    noa3 = noa*(noa-1)/2
    nob3 = nob*(nob-1)/2
    nva3 = nva*(nva-1)/2
    nvb3 = nvb*(nvb-1)/2
    ntt  = nbasis*(nbasis+1)/2    


  end subroutine get_size_variables
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_vabs_ao( Qstore )
    

    implicit none    
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables    
    

    if ( Qstore ) then

       allocate( vabsao(ntt) )
       write(iout,form1e) 'vabsao', ntt  ;   call track_mem( ntt )
       
       read(10, '(A)') cskip
       read(10,*) ( vabsao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    

  end subroutine read_vabs_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipx_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables

    
    if ( Qstore ) then

       allocate( dipxao(ntt) )
       write(iout,form1e) 'dipxao', ntt  ;  call track_mem( ntt )       
       
       read(10,'(A)') cskip
       read(10,*) ( dipxao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    

  end subroutine read_dipx_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipy_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables


    if ( Qstore ) then

       allocate( dipyao(ntt) )
       write(iout,form1e) 'dipyao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       read(10,*) ( dipyao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else

       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if

    
  end subroutine read_dipy_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_dipz_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables

    
    if ( Qstore ) then

       allocate( dipzao(ntt) )
       write(iout,form1e) 'dipzao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       read(10,*) ( dipzao(i) , i=1, ntt )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
    else

       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    

  end subroutine read_dipz_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socx_ao( Qstore )
    

    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables


    if ( Qstore ) then

       !allocate( socxao(ntt), r_array(ntt) )
       allocate( socxao(nbasis,nbasis) ) 
       write(iout,form1e) 'socxao', ntt  ;  call track_mem( ntt )       

       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             socxao(i,j) = dcmplx( 0.d0, rdum )
             socxao(j,i) = dconjg( socxao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !socxao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
       !deallocate( r_array )

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    
    
  end subroutine read_socx_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socy_ao( Qstore ) 


    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables
    

    if ( Qstore ) then

       !allocate( socyao(ntt), r_array(ntt) )
       allocate( socyao(nbasis,nbasis) )
       write(iout,form1e) 'socyao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             socyao(i,j) = dcmplx( 0.d0, rdum )
             socyao(j,i) = dconjg( socyao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !socyao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
       
       !deallocate( r_array )

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    

  end subroutine read_socy_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_socz_ao( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    real(8), allocatable :: r_array(:)
    real(8) :: rdum
    integer(8) :: i, j


    call get_size_variables


    if ( Qstore ) then

       !allocate( soczao(ntt), r_array(ntt) )
       allocate( soczao(nbasis,nbasis) )
       write(iout,form1e) 'soczao', ntt  ;  call track_mem( ntt )
       
       read(10,'(A)') cskip
       do i=1, nbasis
          do j=1, i
             read(10,*) rdum
             soczao(i,j) = dcmplx( 0.d0, rdum )
             soczao(j,i) = dconjg( soczao(i,j) )
          end do
       end do
       !read(10,*) ( r_array(i), i=1, ntt )
       !soczao(:) = dcmplx( 0.d0, r_array(:) )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

       !deallocate( r_array ) 

    else
       
       read(10, '(A)') cskip
       read(10,*) ( rskip , i=1, ntt )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if

    
  end subroutine read_socz_ao
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_orben( Qstore )
    
    
    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i

    
    call get_size_variables

    
    if ( Qstore ) then

       allocate( orben(norb) )
       write(iout,form1e) 'orben', norb  ;  call track_mem( norb )
       
       read(10,'(A)') cskip
       read(10,*) ( orben(i) , i=1, norb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else
       
       read(10,'(A)') cskip
       read(10,*) ( rskip , i=1, norb )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))
       
    end if
    
    
  end subroutine read_orben
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_cmo( Qstore )


    implicit none
    logical, intent(in) :: Qstore
    integer(8) :: i


    call get_size_variables


    if ( Qstore ) then

       allocate( cmo_a(nrorb*nbasis) )
       write(iout,form1e) 'cmo_a', nrorb*nbasis  ;  call track_mem( nrorb*nbasis )
       
       if ( unrestricted ) then
          allocate( cmo_b(nrorb*nbasis) ) 
          write(iout,form1e) 'cmo_b', nrorb*nbasis  ;  call track_mem( nrorb*nbasis )
       end if
       
       read(10,'(A)') cskip
       read(10,*) ( cmo_a(i) , i=1, nbasis*nrorb )       
       if (unrestricted) read(10,*) ( cmo_b(i) , i=1, nbasis*nrorb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    else
       
       read(10,'(A)') cskip
       read(10,*) ( rskip , i=1, nbasis*nrorb )       
       if (unrestricted) read(10,*) ( rskip, i=1, nbasis*nrorb )
       write(iout,'(A)') ' finished skipping '//trim(adjustl(cskip))

    end if
    
  end subroutine read_cmo
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_custom_ham

    implicit none
    integer(8) :: i, j


    !: write out warning
    write(iout,'(A)') ' directly reading in the Hamiltonian from TDCI.dat '
    write(iout,'(A)') ' WARNING:  change subroutine read_custom_ham if necessary'
    write(iout,'(A)') ' '
    flush(iout)

    
    !: change if not TDA elements read from Gaussian output
    cis_vec(1) = 0.d0 


    !: read from TDCI.dat
    read(10,'(A)') cskip
    do i=2, nstates
       read(10,*) ( cis_vec( (i-1)*nstates+j ), j=2,nstates )
    end do


    !: make sure Ham is Hermetian
    do i=1, nstates
       do j=(i+1), nstates
          cis_vec( (i-1)*nstates+j ) = 0.5d0 * ( cis_vec((i-1)*nstates+j) + cis_vec((j-1)*nstates+i) )
          cis_vec( (j-1)*nstates+i ) = cis_vec( (i-1)*nstates+j )
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_custom_ham
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_1
    !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( dijabAA(nva3,noa3) )            ; dijabAA = 0.d0 ! Bucket 2    
    write(iout,form2e) 'dijabAA', nva3, noa3  ; call track_mem( nva3*noa3 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( dijabAA(j,i), j=1, nva3 ), i=1, noa3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_1
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_2
    !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( dijabAB(nobnvb,noanva) )            ;  dijabAB = 0.d0 ! Bucket 2    
    write(iout,form2e) 'dijabAB', nobnvb, noanva  ;  call track_mem( nobnvb*noanva )
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( dijabAB(j,i), j=1, nobnvb ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_3
    !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB
    
    implicit none
    integer(8) :: i, j 
    

    call get_size_variables

    
    allocate( dijabBB(nvb3,nob3) )            ;  dijabBB = 0.d0 ! Bucket 2    
    write(iout,form2e) 'dijabBB', nvb3, nob3  ;  call track_mem( nvb3*nob3 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( dijabBB(j,i), j=1, nvb3 ), i=1, nob3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_3
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_4
    !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L,IJ.LE.KL  IJKL
    

    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( dijklAA(noa3,noa3) )            ;  dijklAA = 0.d0 
    write(iout,form2e) 'dijklAA', noa3, noa3  ;  call track_mem( noa3*noa3 )
    

    read(10,'(A)') cskip
    do i=1, noa3
       do j=i, noa3
          read(10,*) dijklAA(j,i)
          dijklAA(i,j) = dijklAA(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_4
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_5
    !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( diajbAA(noanva,noanva) )            ;  diajbAA = 0.d0 ! Bucket 5 
    write(iout,form2e) 'diajbAA', noanva, noanva  ;  call track_mem( noanva*noanva )
    

    read(10,'(A)') cskip
    read(10,*) ( ( diajbAA(j,i), j=1, noanva ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_5
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_6
    !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( dijklAB(nob2,noa2) )            ;  dijklAB = 0.d0 
    write(iout,form2e) 'dijklAB', nob2, noa2  ;  call track_mem( nob2*noa2 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( dijklAB(j,i), j=1, nob2 ), i=1, noa2 )
    write(iout,'(A)') ' finishd reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_6
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_7
    !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( diajbAB(nvb2,noa2) )            ;  diajbAB = 0.d0 
    write(iout,form2e) 'diajbAB', nvb2, noa2  ;  call track_mem( nvb*noa2 )


    read(10,'(A)') cskip
    read(10,*) ( ( diajbAB(j,i), j=1, nvb2 ), i=1, noa2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_7
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_8
    !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( diajbBA(nva2,nob2) )            ;  diajbBA = 0.d0 
    write(iout,form2e) 'diajbBA', nva2, nob2  ;  call track_mem( nva2*nob2 )

    
    read(10,'(A)') cskip
    read(10,*) ( ( diajbBA(j,i), j=1, nva2 ), i=1, nob2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_8
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_9
    !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( dijklBB(nob3,nob3) )            ;  dijklBB = 0.d0 
    write(iout,form2e) 'dijklBB', nob3, nob3  ;  call track_mem( nob3*nob3 )


    read(10,'(A)') cskip
    do i=1, nob3
       do j=i, nob3
          read(10,*) dijklBB(j,i)
          dijklBB(i,j) = dijklBB(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_9
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!  
  subroutine read_bucket_10
    !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( diajbBB(nobnvb,nobnvb) )            ;  diajbBB = 0.d0 ! Bucket 10
    write(iout,form2e) 'diajbBB', nobnvb, nobnvb  ;  call track_mem( nobnvb*nobnvb ) 

    
    if(unrestricted) then
       read(10,'(A)') cskip
       read(10,*) ( ( diajbBB(j,i), j=1, nobnvb ), i=1, nobnvb )
       write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    end if

    
  end subroutine read_bucket_10
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_11
    !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( dijkaAA(noanva,noa3))             ;  dijkaAA = 0.d0 
    write(iout,form2e) 'dijkaAA', noa3, noanva  ;  call track_mem( noa3*noanva )


    read(10,'(A)') cskip
    read(10,*) ( ( dijkaAA(j,i), j=1, noanva ), i=1, noa3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_11
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_12
    !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( dijkaAB(nobnvb,noa2))             ;  dijkaAB = 0.d0 
    write(iout,form2e) 'dijkaAB', nobnvb, noa2  ;  call track_mem( noa2*nobnvb )
    

    read(10,'(A)') cskip
    read(10,*) ( ( dijkaAB(j,i), j=1, nobnvb ), i=1, noa2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_12
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_13
    !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( dijkaBA(noanva,nob2))             ;  dijkaBA = 0.d0 
    write(iout,form2e) 'dijkaBA', noanva, nob2  ;  call track_mem( nob2*noanva )
    

    read(10,'(A)') cskip
    read(10,*) ( ( dijkaBA(j,i), j=1, noanva ), i=1, nob2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_13
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_14
    !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( dijkaBB(nobnvb,nob3))             ;  dijkaBB = 0.d0 
    write(iout,form2e) 'dijkaBB', nobnvb, nob3  ;  call track_mem( nob3*nobnvb )


    read(10,'(A)') cskip
    read(10,*) ( ( dijkaBB(j,i), j=1, nobnvb ), i=1, nob3 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_14
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_15
    !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables

    
    allocate( diabcAA(nva3,noanva))             ;  diabcAA = 0.d0 
    write(iout,form2e) 'diabcAA', noanva, nva3  ;  call track_mem( noanva*nva3 )
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( diabcAA(j,i), j=1, nva3 ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_15
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_16
    !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC


    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( diabcAB(nvb2,noanva))             ;  diabcAB = 0.d0
    write(iout,form2e) 'diabcAB', nvb2, noanva  ;  call track_mem( noanva*nvb2 ) 
    
    
    read(10,'(A)') cskip
    read(10,*) ( ( diabcAB(j,i), j=1, nvb2 ), i=1, noanva )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    
    
  end subroutine read_bucket_16
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_17
    !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables    

    
    allocate( diabcBA(nva2,nobnvb))             ;  diabcBA = 0.d0 
    write(iout,form2e) 'diabcBA', nva2, nobnvb  ;  call track_mem( nobnvb*nva2 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( diabcBA(j,i), j=1, nva2 ), i=1, nobnvb )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_17
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_18
    !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC

    
    implicit none
    integer(8) :: i, j 

    
    call get_size_variables
    

    allocate( diabcBB(nvb3,nobnvb))              ;  diabcBB = 0.d0 
    write(iout,form2e)  'diabcBB', nvb3, nobnvb  ;  call track_mem( nobnvb*nvb3 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( diabcBB(j,i), j=1, nvb3 ), i=1, nobnvb )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_18
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_19
    !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD


    implicit none
    integer(8) :: i, j 

    
    call get_size_variables


    allocate( dabcdAA(nva3,nva3) )            ;  dabcdAA = 0.d0 
    write(iout,form2e) 'dabcdAA', nva3, nva3  ;  call track_mem( nva3*nva3 )
    

    read(10,'(A)') cskip
    do i=1, nva3
       do j=i, nva3
          read(10,*) dabcdAA(j,i)
          dabcdAA(i,j) = dabcdAA(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_19
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_20
    !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD

    implicit none
    integer(8) :: i, j 


    call get_size_variables

    
    allocate( dabcdAB(nvb2,nva2) )            ;  dabcdAB = 0.d0 
    write(iout,form2e) 'dabcdAB', nvb2, nva2  ;  call track_mem( nva2*nvb2 )
    

    read(10,'(A)') cskip
    read(10,*) ( ( dabcdAB(j,i), j=1, nvb2 ), i=1, nva2 )
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))
    

  end subroutine read_bucket_20
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine read_bucket_21
    !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD


    implicit none
    integer(8) :: i, j 


    call get_size_variables


    allocate( dabcdBB(nvb3,nvb3) )            ;  dabcdBB = 0.d0
    write(iout,form2e) 'dabcdBB', nvb3, nvb3  ;  call track_mem( nvb3*nvb3 )
    
    
    read(10,'(A)') cskip
    do i=1, nvb3
       do j=i, nvb3
          read(10,*) dabcdBB(j,i)
          dabcdBB(i,j) = dabcdBB(j,i)
       end do
    end do
    write(iout,'(A)') ' finished reading '//trim(adjustl(cskip))

    
  end subroutine read_bucket_21
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  !: GET INTEGRALS FOR EACH BUCKET.  LOOK get last function INTEGER(8) GET_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabAA( i, j, a, b )
    !: Bucket  1        (IJ//AB)       AAAA    I.lt.J,A.lt.B     IJAB

    integer(8), intent(in) :: i, j, a, b
    integer(8) :: ij, ab

    !: <ij||ab> DijabAA( ab, ij )

    ij = get_index( i,j, noa, 'lt' )
    ab = get_index( a,b, nva, 'lt' )

    get_dijabAA = dijabAA(ab,ij)

  end function get_dijabAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabBB( i, j, a, b )
    !: Bucket  3        (IJ//AB)       BBBB    I.lt.J,A.lt.B     IJAB

    integer(8), intent(in) :: i, j, a, b
    integer(8) :: ij, ab

    !: <IJ||AB> DijabBB( AB, IJ )

    IJ = get_index( I, J, nob, 'lt' )
    AB = get_index( A, B, nvb, 'lt' )
    
    get_dijabBB = dijabBB(AB,IJ)
    
  end function get_dijabBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijabAB( i, j, a, b )
    !: Bucket  2        (IJ//AB)       ABAB        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <iJ||aB> DijabAB( JB, ia )

    ia = get_index( i,a, nva, 'all' )
    jb = get_index( j,b, nvb, 'all' )

    get_dijabAB = dijabAB(jb,ia)

  end function get_dijabAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklAA( i, j, k, l )
    !: Bucket  4        (IJ//KL)       AAAA    I<J,K<L  IJKL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: kl, ij
    real(8) :: sign

    !: <ij||kl> DijklAA( k<l, i<j )

    sign = 1.d0
    if ( l.lt.k ) sign = - 1.d0 * sign
    if ( j.lt.i ) sign = - 1.d0 * sign

    kl = get_index( k,l, noa, 'lt' )
    ij = get_index( i,j, noa, 'lt' )

    get_dijklAA = sign * dijklAA(kl,ij)

  end function get_dijklAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbAA( i, a, j, b )
    !: Bucket  5        (IA//JB)       AAAA        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <ia||jb> DiajbAA( jb, ia )

    ia = get_index( i,a, nva, 'all' )
    jb = get_index( j,b, nva, 'all' )

    get_diajbAA = diajbAA(jb,ia)

  end function get_diajbAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklAB( i, j, k, l )
    !: Bucket  6        (IJ//KL)       ABAB    I.le.K,J.le.L     IKJL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: jl, ik

    !: <iJ||kL> DijklAB( J<=L,i<=k )

    jl = get_index( j, l, nob, 'le' )
    ik = get_index( i, k, noa, 'le' )

    get_dijklAB = dijklAB(jl,ik)

  end function get_dijklAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbAB( i, a, j, b )
    !: Bucket  7        (IA//JB)       ABAB    I.le.J,A.le.B     IJAB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ab, ij

    !: <iA||jB> DiajbAB( A<=B, i<=j )

    ab = get_index( a,b, nvb, 'le' )
    ij = get_index( i,j, noa, 'le' )

    get_diajbAB = diajbAB(ab,ij)

  end function get_diajbAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbBA( i, a, j, b )
    !: Bucket  8        (IA//JB)       BABA    I.le.J,A.le.B     IJAB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ab, ij

    !: <Ia||Jb> DiajbBA( a<=b, I<=J )

    ab = get_index( a,b, nva, 'le' )
    ij = get_index( i,j, nob, 'le' )

    get_diajbBA = diajbBA(ab,ij)

  end function get_diajbBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijklBB( i, j, k, l )
    !: Bucket  9        (IJ//KL)       BBBB    I<J,K<L,IJ.LE.KL  IJKL

    integer(8), intent(in) :: i, j, k, l
    integer(8) :: kl, ij
    real(8) :: sign

    !: <IJ||KL> DijklBB( K<L,I<J )

    sign = 1.d0
    if ( l.lt.k ) sign = -1.d0 * sign
    if ( j.lt.i ) sign = -1.d0 * sign

    kl = get_index( k,l, nob, 'lt' )
    ij = get_index( i,j, nob, 'lt' )

    get_dijklBB = sign * dijklBB(kl,ij)

  end function get_dijklBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diajbBB( i, a, j, b )
    !: Bucket 10        (IA//JB)       BBBB        ALL           IAJB

    integer(8), intent(in) :: i, a, j, b
    integer(8) :: ia, jb

    !: <IA||JB> DiajbBB( JB, IA )

    ia = get_index( i,a, nvb, 'all' )
    jb = get_index( j,b, nvb, 'all' )

    get_diajbBB = diajbBB(jb,ia)

  end function get_diajbBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaAA( i, j, k, a )
    !: Bucket 11        (IJ//KA)       AAAA     I<J, ALL KA      IJKA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ka, ij
    real(8) :: sign

    !: <ij||ka> DijkaAA( ka, i<j )

    sign = 1.d0
    if ( j.lt.i ) sign = -1.d0 * sign

    ka = get_index( k, a, nva, 'all' )
    ij = get_index( i, j, noa, 'lt' )

    get_dijkaAA = sign * dijkaAA(ka,ij)

  end function get_dijkaAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaAB( i, j, k, a )
    !: Bucket 12        (IJ//KA)       ABAB    I.LE.K, ALL JA    IKJA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ja, ik

    !: <iJ||kA> DijkaAB( JA, i<=k )

    ja = get_index( j,a, nvb, 'all' )
    ik = get_index( i,k, noa, 'le' )

    get_dijkaAB = dijkaAB(ja,ik)

  end function get_dijkaAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaBA( i, j, k, a )
    !: Bucket 13        (IJ//KA)       BABA    I.LE.K, ALL JA    IKJA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ja, ik

    !: <Ij||Ka> DijkaBA( ja, I<=K )

    ja = get_index( j, a, nva, 'all' )
    ik = get_index( i, k, nob, 'le' )

    get_dijkaBA = dijkaBA(ja,ik)

  end function get_dijkaBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dijkaBB( i, j, k, a )
    !: Bucket 14        (IJ//KA)       BBBB     I<J, ALL KA      IJKA

    integer(8), intent(in) :: i, j, k, a
    integer(8) :: ka, ij
    real(8) :: sign

    !: <IJ||KA> DijkaBB( KA, I<J )

    sign = 1.d0
    if ( j.lt.i ) sign = -1.d0  * sign

    ka = get_index( k, a, nvb, 'all' )
    ij = get_index( i, j, nob, 'lt' )

    get_dijkaBB = sign * dijkaBB(ka,ij)

  end function get_dijkaBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcAA( i, a, b, c )
    !: Bucket 15        (IA//BC)       AAAA     ALL IA, B<C      IABC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: bc, ia
    real(8) :: sign

    !: <ia||bc> DiabcAA( b<c, ia )

    sign = 1.d0
    if( c.lt.b ) sign = -1.d0 * sign

    bc = get_index( b, c, nva, 'lt'  )
    ia = get_index( i, a, nva, 'all' )

    get_diabcAA = sign * diabcAA(bc,ia)

  end function get_diabcAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcAB( i, a, b, c )
    !: Bucket 16        (IA//BC)       ABAB    ALL IB, A.LE.C    IBAC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: ib, ac

    !: <iA||bC> DiabcAB( A<=C, ib)

    ac = get_index( a, c, nvb, 'le'  )
    ib = get_index( i, b, nva, 'all' )

    get_diabcAB = diabcAB(ac,ib)

  end function get_diabcAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcBA( i, a, b, c )
    !: Bucket 17        (IA//BC)       BABA    ALL IB, A.LE.C    IBAC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: ac, ib

    !: <Ia||Bc> DiabcBA( a<=c, IB )

    ac = get_index( a, c, nva, 'le' )
    ib = get_index( i, b, nvb, 'all' )

    get_diabcBA = diabcBA(ac,ib)

  end function get_diabcBA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_diabcBB( i, a, b, c )
    !: Bucket 18        (IA//BC)       BBBB     ALL IA, B<C      IABC

    integer(8), intent(in) :: i, a, b, c
    integer(8) :: bc, ia
    real(8) :: sign

    !: <IA||BC> DiabcBB( B<C, IA )

    sign = 1.d0
    if ( c.lt.b ) sign = -1.d0 * sign

    bc = get_index( b, c, nvb, 'lt'  )
    ia = get_index( i, a, nvb, 'all' )

    get_diabcBB = sign * diabcBB(bc,ia)

  end function get_diabcBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdAA( a, b, c, d )
    !: Bucket 19        (AB//CD)       AAAA    A<B,C<D,AB.LE.CD  ABCD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: ab, cd
    real(8) :: sign

    !: <ab||cd> DabcdAA( c<d, a<b )

    sign = 1.d0
    if ( d.lt.c ) sign = -1.d0 * sign
    if ( b.lt.a ) sign = -1.d0 * sign

    ab = get_index( a, b, nva, 'lt' )
    cd = get_index( c, d, nva, 'lt' )

    get_dabcdAA = sign * dabcdAA(cd,ab)

  end function get_dabcdAA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdAB( a, b, c, d )
    !: Bucket 20        (AB//CD)       ABAB    A.le.C,B.le.D     ACBD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: ac, bd

    !: <aB||cD> DabcdAB( B<=D, a<=c )

    bd = get_index( b, d, nvb, 'le' )
    ac = get_index( a, c, nva, 'le' )

    get_dabcdAB = dabcdAB(bd,ac)

  end function get_dabcdAB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  real(8) function get_dabcdBB( a, b, c, d )
    !: Bucket 21        (AB/CD)        BBBB    A<B,C<D,AB.LE.CD  ABCD

    integer(8), intent(in) :: a, b, c, d
    integer(8) :: cd, ab
    real(8) :: sign

    !: <AB||CD> DabcdBB( C<D, A<B )

    sign = 1.d0
    if ( d.lt.c ) sign = -1.d0 * sign
    if ( b.lt.a ) sign = -1.d0 * sign

    cd = get_index( c, d, nvb, 'lt' )
    ab = get_index( a, b, nvb, 'lt' )

    get_dabcdBB = sign * dabcdBB(cd,ab)

  end function get_dabcdBB
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  integer(8) function get_index( i, j, nn, opt )

    integer(8), intent(in) :: i,j, nn
    character(*), intent(in) :: opt

    integer(8) :: big, small

    if ( trim(opt).eq.'lt' ) then
       big   = max( i,j )
       small = min( i,j )
       get_index = (small-1) * nn - small*(small-1)/2 + (big-small)
       return
    else if ( trim(opt).eq.'le' ) then
       big   = max( i,j )
       small = min( i,j )
       get_index = (small-1) * nn - (small-1)*(small-2)/2 + (big-small) + 1
       return
    else if ( trim(opt).eq.'all' ) then
       get_index = (i-1) * nn + j
       return
    end if

  end function get_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module read_integrals
