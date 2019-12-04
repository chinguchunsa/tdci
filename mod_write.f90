module write_info 


  use global_variables
  implicit none


  !: contains subroutine write_specifics1
  !: contains subroutine write_specifics2
  !: contains subroutine write_summary
  !: contains subroutine write_mo_energies
  !: contains subroutine write_field_shape
  !: contains subroutine write_ham0
  !: contains subroutine write_mo_elements
  !: contains subroutine save_restart_bin

  
contains
  !: --------------------------- :!
  !: subroutine write_specifics1 :!
  !: --------------------------- :!
  subroutine write_specifics1
    

    implicit none

    real(8)    :: fstrnth
    integer(8) :: i, ifield, idir, ij
    

    call write_header( 'write_specifics1','write_info','enter' )
    

    !: job title taken from TDCI.dat
    write(iout,'(A)') ' JOB TITLE: '//trim(job_title)

    !: molecule info
    write(iout,"(5x,'charge = ',i0,'  multiplicity = ',i0,'  natoms = ',i0)") charge, mult, natoms
    write(iout,"(5x,'coordinates in Angstroms')")

    do iatom=1, natoms
       write(iout,"(5x,a4,2x,3(f17.11,2x))") &
            myatom(iatom),          &
            xcoord(iatom)*bohr2ang, &
            ycoord(iatom)*bohr2ang, &
            zcoord(iatom)*bohr2ang
    end do
    
    write(iout,'(A)') '     ground state dipole in au:'
    write(iout,"(5x,'xdip = ',f10.4,'  ydip = ',f10.4,'  zdip = ',f10.4)") dipx00, dipy00, dipz00
    write(iout,'(A)') '     ground state dipole in Debye:'
    write(iout,"(5x,'xdip = ',f10.4,'  ydip = ',f10.4,'  zdip = ',f10.4)") &
         dipx00*audip2debye, dipy00*audip2debye, dipz00*audip2debye

    !: jobtype info
    if( unrestricted )      write(iout,"(A)") ' unrestricted'
    if( .not.unrestricted ) write(iout,"(A)") ' restricted'

    select case( trim(jobtype) )
    case( flag_cis)   ; write(iout,'(A)') ' CIS'
    case( flag_tda)   ; write(iout,'(A)') ' TDA'
    case( flag_ip)    ; write(iout,'(A)') ' CISD-IP'
    case( flag_cisd ) ; write(iout,'(A)') ' CISD'
    end select
    
    !: system size info   
    write(iout,"(5x,'nbasis  = ',i0)") nbasis
    write(iout,"(5x,'noa     = ',i0)") noa
    write(iout,"(5x,'nva     = ',i0)") nva
    write(iout,"(5x,'nob     = ',i0)") nob
    write(iout,"(5x,'nvb     = ',i0)") nvb
    write(iout,"(5x,'norb    = ',i0)") norb
    write(iout,"(5x,'nrorb   = ',i0)") nrorb
    write(iout,"(5x,'noanva  = ',i0)") noanva
    write(iout,"(5x,'nobnvb  = ',i0)") nobnvb
    write(iout,"(5x,'nstates = ',i0)") nstates
    
    !: field info
    write( iout, '(A)' )  ' FIELD VARIABLES'
    do i=1, nemax
       fstrnth = tdciresults((i-1)*ndir+1)%fstrength0
       write(iout,100) i, fstrnth, fstrnth**2 * auflux2si
    end do
100 format( 5x,'max field strength',i3,' = ',f9.5,'  au',1x,'laser intensity ',es12.4,' W/cm2' )

    if ( envelope.eq.'cirl' .or. envelope.eq.'cirr' ) then
       write(iout,"(5x, 'nonlinear light')") 
       write(iout,"(5x, 'EULER angle euler    = ',f10.4)") euler
       write(iout,"(5x, 'ellipticity, E1/E2 = ',f10.4)")  ellipt
    end if
    
    if ( envelope.ne.'none' .and. envelope.ne.'stat' ) then
       write(iout,"(5x,'laser frequency       = ',f10.4,' au',es12.4,' Hz')") omega, omega*s2autime
       write(iout,"(5x,'laser wavelength      = ',f10.4,' nm')")              light *2.d0*pi / omega * autime2s * m2nm
       write(iout,"(5x,'pulse duration        = ',f10.4,' au',f12.4,' fs')")  field_duration * dble(dt), field_duration * dble(dt) *autime2s * s2fs

       select case ( trim(envelope) )
       case( 'none' ) ; write(iout,"(5x,'envelope none ',7x,' =   null pulse')")      
       case( 'cos2' ) ; write(iout,"(5x,'envelope cos2 ',7x,' =   cosine pulse')")    
       case( 'gaus' ) ; write(iout,"(5x,'envelope gaus ',7x,' =   gaussian pulse')")  
       case( 'stat' ) ; write(iout,"(5x,'envelope stat ',7x,' =   static pulse')")    
       case( 'band' ) ; write(iout,"(5x,'envelope band ',7x,' =   Bandrauk pulse')")  
       case( 'sin2' ) ; write(iout,"(5x,'envelope sin2 ',7x,' =   squared pulse')")   
       case( 'trap' ) ; write(iout,"(5x,'envelope trap ',7x,' =   trapezoidal pulse')")  
       case( 'cirl' ) ; write(iout,"(5x,'envelope cirl ',7x,' =   left circular sine squared pulse')") 
       case( 'cirr' ) ; write(iout,"(5x,'envelope cirr ',7x,' =   right circular sine squared pulse')")  
       end select

       write(iout,"(5x,'# of laser cycles     = ',i10,' cycles')")               ncyc
       write(iout,"(5x,'# timesteps/cycle     = ',i10)")                         int(period)
       write(iout,"(5x,'# timesteps for pulse = ',i10)")                         int(field_duration)
       write(iout,"(5x,'write results every     ',i10,  ' step   (total ',i0,' every ',f6.4,' fs)')") outstep, ndata, dble(outstep)*dt*autime2s* s2fs
       write(iout,"(5x,'timestep size         = ',f10.3,' au',f12.4,' fs')")     dt, dt * autime2s* s2fs
       write(iout,"(5x,'TOTAL SIMULATION TIME = ',f10.3,' au',f12.4,' fs (nstep=',i0,')')") &
            dt*dble(nstep), dt*dble(nstep)*autime2s*s2fs, nstep
    else
       write(iout,"(5x,'pulse duration        = ',f10.4,' au',f12.4,' fs')") field_duration * dble(dt), field_duration * dble(dt) *autime2s * s2fs
       write(iout,"(5x,'timestep size         = ',f10.3,' au',f12.4,' fs')")    dt, dt * autime2s* s2fs
       write(iout,"(5x,'TOTAL SIMULATION TIME = ',f10.3,' au',f12.4,' fs (nstep=',i0,')')") &
            dt*dble(nstep), dt*dble(nstep)*autime2s*s2fs, nstep
    end if

    !: field direction
    write( iout,'(A)' ) ' FIELD DIRECTIONS'
    if ( envelope.ne.'cirl' .and. envelope.ne.'cirr' ) then

       write(iout,50) '#' , 'E_0(au)' , 'x0' , 'y0' , 'z0' , 'theta0' , 'phi0'
       do i=1, nemax*ndir
          write(iout,60) i, tdciresults(i)%fstrength0, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%theta0, tdciresults(i)%phi0
       end do
       
    else
       write(iout,70) '#' , 'E_0(au)' , 'E_1(au)' , 'E_2(au)' ,           &
            'x0' , 'y0' , 'z0' , 'x1' , 'y1' , 'z1' , 'x2' , 'y2' , 'z2' , &
            'theta0' , 'phi0' , 'theta1' , 'phi1' , 'theta2' , 'phi2'
       do i=1, nemax*ndir
          write(iout,80) i, &
               tdciresults(i)%fstrength0, tdciresults(i)%fstrength1, tdciresults(i)%fstrength2, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%x1, tdciresults(i)%y1, tdciresults(i)%z1, &
               tdciresults(i)%x2, tdciresults(i)%y2, tdciresults(i)%z2, &
               tdciresults(i)%theta0, tdciresults(i)%phi0, &
               tdciresults(i)%theta1, tdciresults(i)%phi1, &
               tdciresults(i)%theta2, tdciresults(i)%phi2
       end do

    end if
       
50  format( 1x, a4, a8,   2x, ,a9  ,1x,a9,  1x,a9,  1x,2x,a9  ,1x,a9,  1x )
60  format( 1x, i4, f8.4, 2x, ,f9.4,1x,f9.4,1x,f9.4,1x,2x,f9.4,1x,f9.4,1x )
70  format( 1x,a4,a8,  2x,a8,  2x,a8,  2x,3(a9  ,1x,a9  ,1x,a9  ,1x,2x), 3(a9,  1x,a9,  1x,2x) )
80  format( 1x,i4,f8.4,2x,f8.4,2x,f8.4,2x,3(f9.4,1x,f9.4,1x,f9.4,1x,2x), 3(f9.4,1x,f9.4,1x,2x) )
       

    call write_header( 'write_specifics1','write_info','leave' )
    

  end subroutine write_specifics1
  !: ----------------------------- :!
  !: subroutine write_specific2    :!
  !: ----------------------------- :!
  subroutine write_specifics2


    implicit none
    integer(8) :: i, ii
    
    
    if ( trim(jobtype).eq.flag_soc )  go to 78
    

    call write_header( 'write_specifics2','write_info','enter' )


    write(iout,'(A)') ' NSTUSE' 
    write(iout,"(5x,'NstUse          = ',i12)") nstuse
    write(iout,"(5x,'Nstates         = ',i12)") nstates
    write(iout,"(5x,'cis_eig(nstuse) = ',f12.6,' au',f12.6,' eV')") cis_eig(nstuse), cis_eig(nstuse)*au2eV


    write(iout,'(A)') ' ENERGIES and EXPECTATION VALUES'
    write(iout,40) '#', '<i|H0|i> eV', '<i|mu_x|i> au', '<i|mu_y|i> au', '<i|mu_z|i> au', '<i|Vabs|i> au'
    do i=1, 20
       ii = (i-1)*nstuse + i
       write(iout,50) i, cis_eig(i)*au2eV, tdx(ii), tdy(ii), tdz(ii), abp(ii)
    end do

    write(iout,'(A)') ' ENERGIES AND TRANSITION MATRIX ELEMENTS'
    write(iout,40) '#', '<i|H0|i> eV', '<i|mu_x|0> au', '<i|mu_y|i> au', '<i|mu_z|i> au', '<i|Vabs|i> au'
    do i=2, 20
       write(iout,50) i, cis_eig(i)*au2eV, tdx(i), tdy(i), tdz(i), abp(i)
    end do

40  format( a5,5(1x,a15) )
50  format( i5,5(1x,f15.10) )

    call write_header( 'write_specifics2','write_info','leave' )
    return


78  continue
    call Zwrite_specifics2
    
    
  end subroutine write_specifics2
  !: ----------------------------- :!
  !: subroutine Zwrite_specific2    :!
  !: ----------------------------- :!
  subroutine Zwrite_specifics2


    implicit none
    integer(8) :: i, ii
    
    call write_header( 'Zwrite_specifics2','write_info','enter' )


    write(iout,'(A)') ' NSTUSE' 
    write(iout,"(5x,'NstUse          = ',i12)") nstuse
    write(iout,"(5x,'Nstates         = ',i12)") nstates
    write(iout,"(5x,'cis_eig(nstuse) = ',f12.6,' au',f12.6,' eV')") cis_eig(nstuse), cis_eig(nstuse)*au2eV
    

    write(iout,'(A)') ' ENERGIES and EXPECTATION VALUES'
    write(iout,40) '#', '<i|H0|i> eV', &
         'R<i|mu_x|i> au', 'I<i|mu_x|i> au', &
         'R<i|mu_y|i> au', 'I<i|mu_y|i> au', &
         'R<i|mu_z|i> au', 'I<i|mu_z|i> au', &
         'R<i|Vabs|i> au', 'I<i|mu_z|i> au'
    do i=1, 20
       ii = (i-1)*nstuse + i
       write(iout,50) i, cis_eig(i)*au2eV, Ztdx(ii), Ztdy(ii), Ztdz(ii), Zabp(ii)
    end do

    write(iout,'(A)') ' ENERGIES AND TRANSITION MATRIX ELEMENTS'
    write(iout,40) '#', '<i|H0|i> eV', &
         'R<i|mu_x|0> au', 'I<i|mu_x|0> au', &
         'R<i|mu_y|i> au', 'I<i|mu_y|i> au', &
         'R<i|mu_z|i> au', 'I<i|mu_z|i> au', &
         'R<i|Vabs|i> au', 'I<i|Vabs|i> au'
    do i=2, 20
       write(iout,50) i, cis_eig(i)*au2eV, Ztdx(i), Ztdy(i), Ztdz(i), Zabp(i)
    end do

40  format( a5,9(1x,a15) )
50  format( i5,9(1x,f15.10) )

    call write_header( 'Zwrite_specifics2','write_info','leave' )
    
    
  end subroutine Zwrite_specifics2
  !: ----------------------------- :!
  !: SUBROUTINE WRITE_SUMMARY      :!
  !: ----------------------------- :!
  subroutine write_summary

    implicit none

    integer(8) :: ifield
    character(1000) :: cformat

    call write_header( 'write_summary', 'write_info', 'enter' )
    
    write(iout,'(A)') ' SUMMARY'
    cformat = '    #'//'   e_max'//'         x'//'         y'//'         z' //'      norm'//'      dipx'//'      dipy'//'      dipz'

    write(iout,'(A)') trim(cformat)
    do ifield=1, nemax*ndir
       write(iout,"(i5,1x,f7.5,100(f10.4))") ifield, &
            tdciresults(ifield)%fstrength0,       &
            tdciresults(ifield)%x0, &
            tdciresults(ifield)%y0, &
            tdciresults(ifield)%z0, &
            tdciresults(ifield)%norm0, &
            tdciresults(ifield)%dipx,  &
            tdciresults(ifield)%dipy,  &
            tdciresults(ifield)%dipz
    end do
    

    call write_header( 'write_summary', 'write_info', 'leave' )


  end subroutine write_summary
  !: ----------------------------- :!
  !: subroutine write_mo_energies  :! 
  !: ----------------------------- :!
  subroutine write_mo_energies


    implicit none


    character(100) :: myout 
    integer(8) :: i

    
    call write_header( 'write_mo_energies','write_info','enter' )

    
    myout = trim(outputfile)//'_MO_ENERGIES'
    open( unit=100,file=trim(myout) )


    !: restricted MOs
    if ( .not.unrestricted )  then
       write(100,"(' HOMO =', i0, 5x, ' LUMO =', i0)") noa, noa+1
       write(100,"(a10,a20,1x,a20)") ' # MO' , 'restricted (au)' , 'restricted (eV)'
       do i=1, nrorb
          write(100,"(i10,2(f20.10,1x))") i, orben(i), orben(i)*au2eV
       end do
       
    !: unrestricted MOs
    else if ( unrestricted ) then
       write(100,"(' alpha HOMO = ', i0, 5x, ' alpha LUMO = ', i0)") noa, noa+1
       write(100,"(' beta  HOMO = ', i0, 5x, ' beta  LUMO = ', i0)") nob, nob+1
       write(100,"(a10,4(a20,1x))") ' # MO','alpha (au)','alpha (eV)','beta (au)','beta (eV)'
       do i=1, nrorb
          write(100,"(i10,4(f20.10,1x))") i, orben(i), orben(i)*au2eV, orben(i+nrorb), orben(i+nrorb)*au2eV
       end do
    end if

    close(100)
    write(iout,'(A)') " MO orbital energies written out to file "//"'"//trim(myout)//"'"


    call write_header( 'write_mo_energies','write_info','leave' )

    
  end subroutine write_mo_energies
  !: ----------------------------- :!
  !: subroutine write_field_shape  :!
  !: ----------------------------- :!
  subroutine write_field_shape
  
    
    implicit none
    
    
    character(100) :: myout_shape 
    character(100) :: myout_dir   
    integer(8) :: i, istp


    call write_header( 'write_field_shape','write_info','enter' )


    myout_shape = trim(outputfile)//'_FIELD_SHAPE'
    myout_dir   = trim(outputfile)//'_FIELD_DIR'
    
    field_shape : if ( .not.linear ) then

       open( unit=100,file=trim(myout_shape) )

       write(100,"( a10,a10,a10,a20,a20,a20 )") 'step#','time(au)','time(fs)','envelope','E1(t)','E2(t)'
       do istp=1, nstep
          write(100,100) istp, dble(istp)*dt, dble(istp)*dt*autime2s*s2fs, env(istp), &
               1.d0/dsqrt(1.d0+ellipt**2)*fvect2(istp), ellipt/dsqrt(1.d0+ellipt**2)*fvect1(istp)
       end do

       close(100)
       
    else
       
       open( unit=100,file=trim(myout_shape) )

       write(100,"( a10,a10,a10,a20,a20 )") 'step#','time(au)','time(fs)','envelope','E1(t)'       
       do istp=1, nstep
          write(100,200) istp, dble(istp)*dt, dble(istp)*dt*autime2s*s2fs, env(istp), fvect1(istp)
       end do

       close(100)
              
    end if field_shape

100 format( i10,f10.2,f10.4,e20.4,e20.4,e20.4 )
200 format( i10,f10.2,f10.4,e20.4,e20.4 )
    
    
    write(iout, '(A)') " pulse shape written to file "//"'"//trim(myout_shape)//"'"
    

    !: field direction
    open( unit=101, file=trim(myout_dir) )

    field_dir : if ( envelope.ne.'cirl' .and. envelope.ne.'cirr' ) then

       write(101,50) '#' , 'E_0(au)' , 'x0' , 'y0' , 'z0' , 'theta0' , 'phi0'
       do i=1, nemax*ndir
          write(101,60) i, &
               tdciresults(i)%fstrength0, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%theta0, tdciresults(i)%phi0
       end do
       
    else

       write(101,70) '#' , 'E_0(au)' , 'E_1(au)' , 'E_2(au)' ,           &
            'x0' , 'y0' , 'z0' , 'x1' , 'y1' , 'z1' , 'x2' , 'y2' , 'z2' , &
            'theta0' , 'phi0' , 'theta1' , 'phi1' , 'theta2' , 'phi2'
       do i=1, nemax*ndir
          write(101,80) i, &
               tdciresults(i)%fstrength0, tdciresults(i)%fstrength1, tdciresults(i)%fstrength2, &
               tdciresults(i)%x0, tdciresults(i)%y0, tdciresults(i)%z0, &
               tdciresults(i)%x1, tdciresults(i)%y1, tdciresults(i)%z1, &
               tdciresults(i)%x2, tdciresults(i)%y2, tdciresults(i)%z2, &
               tdciresults(i)%theta0, tdciresults(i)%phi0, &
               tdciresults(i)%theta1, tdciresults(i)%phi1, &
               tdciresults(i)%theta2, tdciresults(i)%phi2
       end do
       
    end if field_dir

    
50  format( 1x, a4, a8,   2x, '(',a9  ,1x,a9,  1x,a9,  1x,')',2x,'(',a9  ,1x,a9,  1x,')' )
60  format( 1x, i4, f8.4, 2x, '(',f9.4,1x,f9.4,1x,f9.4,1x,')',2x,'(',f9.4,1x,f9.4,1x,')' )
70  format( 1x,a4,a8,  2x,a8,  2x,a8,  2x,3('(',a9  ,1x,a9  ,1x,a9  ,1x,')',2x), 3(' (',a9,  1x,a9,  1x,')',2x) )
80  format( 1x,i4,f8.4,2x,f8.4,2x,f8.4,2x,3('(',f9.4,1x,f9.4,1x,f9.4,1x,')',2x), 3(' (',f9.4,1x,f9.4,1x,')',2x) )

    
    write(iout, '(A)') " field direction written to file "//"'"//trim(myout_dir)//"'"
    

    call write_header( 'write_field_shape','write_info','leave' )


  end subroutine write_field_shape
  !: ----------------------- :!
  !: SUBROUTINE WRITE_HAM0   :! 
  !: ----------------------- :!
  subroutine write_ham0

    implicit none

    real(8), parameter :: thres  = 0.00d0
    real(8), parameter :: thres1 = thres * dsqrt(2.d0)
    
    integer(8)     :: ia, ii, aa, xx, jj, bb, istate, i, j, a, b
    real(8), allocatable    :: coeffs(:,:)
    complex(8), allocatable :: Zcoeffs(:,:)
    character(5)   :: cii, caa, cxx, aorb, xorb, cjj, cbb
    character(100) :: myout
    
    
    call write_header( 'write_ham0','write_info','enter' )

    myout = trim(outputfile)//'_HAM0'

    open( unit=100,file=trim(myout) )

    if ( trim(jobtype).eq.flag_soc ) then
       allocate( Zcoeffs(nstates,nstates) )
       Zcoeffs = reshape( Zcis_vec, (/ nstates, nstates/) )
    else
       allocate( coeffs(nstates,nstates) )
       coeffs = reshape( cis_vec,(/ nstates, nstates /) )
    end if

    select case ( trim(jobtype) )
    case( flag_cis , flag_tda )
       
       if ( .not.unrestricted ) then
          restricted : do i=1, nstates
             write(100,100) i, cis_eig(i)*au2eV
             do ia=1, nstates
                if( abs(coeffs(ia,i)).ge.thres1 ) &
                     write(100,200) abs(hole_index(ia,1)), abs(part_index(ia,1)), coeffs(ia,i)/sqrt(2.d0)
             end do
          end do restricted
       else
          unrestricted : do i=1, nstates
             write(100,100) i, cis_eig(i)*au2eV
             do ia=1, nstates
                if ( abs(coeffs(ia,i)).ge.thres ) &
                     write(100,200) hole_index(ia,1), part_index(ia,1), coeffs(ia,i)
             end do
          end do unrestricted
       end if
       
    case( flag_soc )
       do i=1, nstates
          write(100,100) i, cis_eig(i)*au2eV
          do ia=1, nstates
             if ( abs(Zcoeffs(ia,i)).ge.thres ) &
                  write(100,200) hole_index(ia,1), part_index(ia,1), Zcoeffs(ia,i)
          end do
       end do
       
    case( flag_ip )

       ip : do i=1, nstates
          write(100,100) i, cis_eig(i)*au2eV
          do ia=1, nstates
             if ( abs(coeffs(ia,i)).ge.thres ) then
                xx = hole_index(ia,1) 
                ii = hole_index(ia,2) 
                aa = part_index(ia,1) 
                write(100,310) xx, ii, aa, coeffs(ia,i)
             end if
          end do
       end do ip

    case( flag_cisd ) 
       
       do ia=1, nstates


          write(100,100) ia, cis_eig(ia)*au2eV

          if( abs(coeffs(1,ia).ge.thres) ) write(100,350) ii, aa, jj, bb, coeffs(1,ia)
          
          istate = 1 
          
          alpha_singles : do i=1, noa
             do a=1, nva
                istate = istate + 1 
                if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) -i, -a, 0, 0, coeffs(istate,ia)
             end do
          end do alpha_singles

          beta_singles : do i=1, nob
             do a=1, nvb
                istate = istate + 1 
                if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) i, a, 0, 0, coeffs(istate,ia)
             end do
          end do beta_singles

          alpha_beta_doubles : do i=1, noa
             do j=1, nob
                do a=1, nva
                   do b=1, nvb
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) j, b, -i, -a, coeffs(istate,ia)
                   end do
                end do
             end do
          end do alpha_beta_doubles

          alpha_alpha_doubles : do i=1, noa
             do j=(i+1), noa
                do a=1, nva
                   do b=(a+1), nva
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) -i, -a, -j, -b, coeffs(istate,ia)
                   end do
                end do
             end do
          end do alpha_alpha_doubles
          
          beta_beta_doubles : do i=1, nob
             do j=(i+1), nob
                do a=1, nvb
                   do b=(a+1), nvb
                      istate = istate + 1 
                      if ( abs(coeffs(istate,ia)).ge.thres ) write(100,350) i, a, -j, -b, coeffs(istate,ia)
                   end do
                end do
             end do
          end do beta_beta_doubles
          
       end do
       
    end select

    close(100)
    
100 format(' Excited state ',i10,' : ',f20.10,' eV ')
200 format( 7x,i3,' -> ',i3,2x,f15.10,2x,f15.10 )
300 format( 7x,i5,' -> ',i5,2x,f10.7 )
310 format( 7x,i5,' -> ','inf',2x,i5,' ->',i5,2x,f10.7 )
350 format( 7x,i5,' -> ',i5,2x,i5,' ->',i5,2x,f10.7 )
    
    write(iout,'(A)') ' finished writing eigen-info to file '//"'"//trim(myout)//"'"

    call write_header( 'write_ham0','write_info','leave' )
    

  end subroutine write_ham0
  !: ---------------------------- :!
  !: SUBROUTINE WRITE_MO_ELEMENTS :!
  !: ---------------------------- :!
  subroutine write_mo_elements
  
    implicit none

    integer(8)     :: i, j
    character(100) :: myout, myout2
    character(1)   :: ov1, ov2 !: occupied or virtual

    
    call write_header( 'write_mo_elements','write_info','enter' )


    if( .not. unrestricted ) then
       myout = trim(outputfile)//'_MO_ELEMENTS'
       open( unit=100,file=trim(myout) )
    else
       myout  = trim(outputfile)//'_MO_ELEMENTS_ALPHA'
       myout2 = trim(outputfile)//'_MO_ELEMENTS_BETA'
       open( unit=100,file=trim(myout) )
       open( unit=200,file=trim(myout2) )
    end if

    
    write(100,50) 'i','j', '<i|mu_x|j>(au)', '<i|mu_y|j>(au)', '<i|mu_z|j>(au)', '<i|Vabs|j>(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=i, (noa+nva)
          ov2 = 'v' ; if ( j.le.noa ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, dipxmoa(j,i), dipymoa(j,i), dipzmoa(j,i), vabsmoa(j,i) 
       end do
    end do
    close(100)
    
    write(iout,'(A)') " MO matrix elements written out to file "//"'"//trim(myout)//"'"

    if ( unrestricted ) then

       write(200,50) 'i','j', '<i|mu_x|j>(au)', '<i|mu_y|j>(au)', '<i|mu_z|j>(au)', '<i|Vabs|j>(au)'
       do i=1, (nob+nvb)
          ov1 = 'v' ; if ( i.le.nob ) ov1 = 'o'
          do j=i, (nob+nvb) 
             ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
             write(200,100) i, ov1, j, ov2, dipxmob(j,i), dipymob(j,i), dipzmob(j,i), vabsmob(j,i)
          end do
       end do
       close(200)
       write(iout,'(A)') " MO matrix elements written out to file "//"'"//trim(myout2)//"'"
       
    end if


    if( trim(jobtype).eq.flag_soc ) call write_soc_elements

             
50  format( 2(a7),4(a20) )
100 format( 2(i5,a2,1x),4(f20.10) )
    
    call write_header( 'write_mo_elements','write_info','leave' )
    

  end subroutine write_mo_elements
  !:-----------------------------:!
  !: SUBROUTINE WRITE_SOC_ELEMENTS
  !:-----------------------------:!
  subroutine write_soc_elements

    implicit none

    integer(8)     :: i, j
    character(100) :: myout
    character(1)   :: ov1, ov2 !: occupied or virtual

    
    call write_header( 'write_soc_elements','write_info','enter' )
    

    myout = trim(outputfile)//'_SOC_ELEMENTS_AA'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_alpha', 'j_alpha', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (noa+nva)
          ov2 = 'v' ; if ( j.le.noa ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoAA(j,i)
       end do
    end do
    close(100)

    write(iout,'(A)') " alpha alpha SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"


    myout = trim(outputfile)//'_SOC_ELEMENTS_BB'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_beta', 'j_beta', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (nob+nvb)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (nob+nvb)
          ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoBB(j,i)
       end do
    end do
    close(100)
    
    write(iout,'(A)') " beta  beta  SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"


    myout = trim(outputfile)//'_SOC_ELEMENTS_AB'
    open( unit=100, file=trim(myout) )
    write(100,50) 'i_alpha', 'j_beta', 'Re_socAA(au)', 'Im_socAA(au)'
    do i=1, (noa+nva)
       ov1 = 'v' ; if ( i.le.noa ) ov1 = 'o'
       do j=1, (nob+nvb)
          ov2 = 'v' ; if ( j.le.nob ) ov2 = 'o'
          write(100,100) i, ov1, j, ov2, socmoAB(j,i)
       end do
    end do
    close(100)

    write(iout,'(A)') " alpha beta  SOC MO matrix elements written out to file "//"'"//trim(myout)//"'"    


50  format( 2(a7,1x),2(a20) )
100 format( 2(i5,a2,1x),2(f20.10) )
    

  end subroutine write_soc_elements
  !: --------------------------- :!
  !: SUBROUTINE SAVE_RESTART_BIN
  !: --------------------------- :!
  subroutine save_restart_bin


    implicit none

    integer(8), parameter :: zero = 0

    character(100) :: myout
    integer(8) :: i, j, ij, iflag


    call write_header( 'save_restart_bin','write_info','enter' )


    myout = trim(outputfile)//'_RESTART.bin'
    open( unit=50,file=trim(myout),form='unformatted' )


    !: restart not available for soc_cis yet
    if ( trim(jobtype).eq.flag_soc ) then
       write(iout,'(A)') ' WARNING:  RESTART NOT AVAILABLE FOR SOC_CIS'
       call write_header( 'save_restart_bin','write_info','leave' )
       return
    end if
    
    
    write(50) nstates, nstuse
    write(50) cis_vec
    write(50) cis_eig
    write(50) tdx
    write(50) tdy
    write(50) tdz
    write(50) abp
    write(50) exp_abp
    
    close(50)
    write( iout,'(A)' ) ' restart binary file written out to file '//"'"//trim(myout)//"'"
    write( iout,'(A)' ) ' to read binary file '
    write( iout,"(5x,'nstates, nstuse')" ) 
    write( iout,"(5x,'cis_vec, cis_eig, tdx, tdy, tdz, abp, exp_abp')") 
    

    myout = trim(outputfile)//'_RESTART_MO.bin'
    open( unit=50, file=trim(myout), form='unformatted' )

    if( unrestricted ) write(50) noa, nva, nob, nvb
    if( .not.unrestricted ) write(50) noa, nva, zero, zero
    write(50) ( vabsmoa(:,i), i=1, nrorb )
    write(50) ( dipxmoa(:,i), i=1, nrorb )
    write(50) ( dipymoa(:,i), i=1, nrorb )
    write(50) ( dipzmoa(:,i), i=1, nrorb )
    if ( unrestricted ) then
       write(50) ( vabsmob(:,i), i=1, nrorb )
       write(50) ( dipxmob(:,i), i=1, nrorb )
       write(50) ( dipymob(:,i), i=1, nrorb )
       write(50) ( dipzmob(:,i), i=1, nrorb )
    end if
    
    close(50)
    write( iout,'(A)' ) ' restart binary file written out to file '//"'"//trim(myout)//"'"
    write( iout,'(A)' ) ' to read binary file '
    write( iout,"(5x,'noa, nva, nob, nvb (nob=nvb=0 for restricted)' )" ) 
    write( iout,"(5x,'vabsmoa, dipxmoa, dipymoa, dipzmoa')")
    write( iout,"(5x,'vabsmob, dipxmob, dipymob, dipzmob')")
    
    call write_header( 'save_restart_bin','write_info','leave' )


  end subroutine save_restart_bin
  !: ------------------------ :!
  !: END MODULE WRITE_INFO    :!
  !: ------------------------ :!
end module write_info
