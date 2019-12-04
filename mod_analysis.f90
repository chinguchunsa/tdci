module analyze_psi
  
  !use global_variables ! <C> use at your risk, beware of private & shared OMP variables
  implicit none

  
contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_NORM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_norm( norm, nstuse, psi )
    
    ! <C> Get norm of a complex vector
    
    real(8), intent(inout) :: norm    
    integer(8), intent(in) :: nstuse
    complex(8), intent(in) :: psi(nstuse)
    
    integer(8) :: i
    real(8)    :: rdum, cdum
    
    
    norm = dot_product( psi, psi )
    norm = dsqrt(norm)

    
  end subroutine get_norm
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_EXPECTATION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_expectation(nstuse,norm,psi,obsv,expect_value)
    
    ! <C> get expectation value of obsv

    implicit none

    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: norm
    real(8)   , intent(in) :: obsv(nstuse*nstuse)
    complex(8), intent(in) :: psi(nstuse)
    real(8), intent(inout) :: expect_value
    
    integer(8) :: i, ii, j, ij
    complex(8) :: psi_i

    
    expect_value = 0.d0
    
    
    ! <C> off-diagonal
    do i = 2, nstuse
       ii = (i-1)*nstuse
       psi_i = dconjg( psi(i) )
       do j = 1, (i-1)
          ij = ii + j
          expect_value = expect_value + obsv(ij) * real( psi_i*psi(j) )
       end do
    end do

    expect_value = 2.d0 * expect_value
    
    ! <C> diagonal
    do i = 1, nstuse
       ii = (i-1)*nstuse + i
       expect_value = expect_value + obsv(ii) * dconjg(psi(i))*psi(i) 
    end do
    
    
    expect_value = expect_value / norm**2
    
    
  end subroutine get_expectation
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_EXPECTATION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_Zexpectation(nstuse,norm,psi,Zobsv,expect_value)
    
    ! <C> get expectation value of obsv

    implicit none

    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: Zobsv(nstuse*nstuse)
    complex(8), intent(in) :: psi(nstuse)
    real(8), intent(inout) :: expect_value
    
    integer(8) :: i, ii, j, ij
    complex(8) :: psi_i, cdum

    
    expect_value = 0.d0
    cdum = dcmplx( 0.d0, 0.d0 )
    
    do i=1, nstuse
       ii = (i-1)*nstuse
       psi_i = dconjg( psi(i) )
       do j=1, nstuse
          ij = ii + j
          cdum = cdum + Zobsv(ij) * psi(j) * psi_i
       end do
    end do
    
    
    expect_value = dble(cdum) / norm**2

    
  end subroutine get_Zexpectation
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET PSI_DET
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_psid(nstuse,nstates,cis_vec,norm,psi,psi_det)

    implicit none

    integer(8), intent(in) :: nstuse, nstates
    real(8),    intent(in) :: cis_vec(nstates*nstates), norm
    complex(8), intent(in) :: psi(nstuse)
    complex(8), intent(inout) :: psi_det(nstates)

    integer(8) :: k, kk
    complex(8) :: psi_k


    psi_det = dcmplx(0.d0,0.d0)


    do k = 1, nstuse
       kk = nstates * (k-1)
       psi_k = psi(k) 
       psi_det(:) = psi_det(:) + cis_vec(kk+1:kk+nstates)*psi_k
    end do

    psi_det = psi_det / norm

    
  end subroutine get_psid
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET PSI_DET
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_Zpsid(nstuse,nstates,Zcis_vec,norm,psi,psi_det)

    implicit none

    integer(8), intent(in) :: nstuse, nstates
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: psi(nstuse), Zcis_vec(nstates*nstates)
    complex(8), intent(inout) :: psi_det(nstates)

    integer(8) :: k, kk
    complex(8) :: psi_k


    psi_det = dcmplx(0.d0,0.d0)


    do k = 1, nstuse
       kk = nstates * (k-1)
       psi_k = psi(k) 
       psi_det(:) = psi_det(:) + Zcis_vec(kk+1:kk+nstates)*psi_k
    end do

    psi_det = psi_det / norm

    
  end subroutine get_Zpsid
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_CIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_cis( hole,noa,norb,nstates,nva,part,pop,psi_det,nob,nvb)
    
    ! <C> Calculate MO occupancies for a CIS vector

    implicit none
    
    integer(8), intent(in) :: noa, norb, nstates, nva
    integer(8), intent(in) :: hole(nstates), part(nstates)
    complex(8), intent(in) :: psi_det(nstates)
    real(8), intent(inout) :: pop(norb)
    integer(8), optional, intent(in) :: nob, nvb

    integer(8) :: i, a, ia, ii, aa
    real(8)    :: psi2
    
    
    pop = 0.d0
    
    if ( present(nob) ) then
       pop(1:noa) = 1.d0
       pop(noa+nva+1:noa+nva+nob) = 1.d0
    else
       pop(1:noa) = 2.d0
    end if
    
    do ia = 2, nstates
       ii   = hole(ia) 
       aa   = part(ia)
       psi2 = dconjg(psi_det(ia))*psi_det(ia)
       if ( ii.lt.0 ) then
          i = -ii
          a = -aa + noa          
          pop(i) = pop(i) - psi2
          pop(a) = pop(a) + psi2
       else
          i = ii + noa + nva
          a = aa + noa + nva + nob
          pop(i) = pop(i) - psi2
          pop(a) = pop(a) + psi2
       end if
    end do

    
  end subroutine pop_cis
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_SOC
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_soc( hole,noa,norb,nstates,nva,part,pop,psi_det,nob,nvb)
    
    ! <C> Calculate MO occupancies for a CIS vector

    implicit none
    
    integer(8), intent(in) :: noa, norb, nstates, nva, nob, nvb
    integer(8), intent(in) :: hole(nstates), part(nstates)
    complex(8), intent(in) :: psi_det(nstates)
    real(8), intent(inout) :: pop(norb)    

    integer(8) :: i, a, ia, ii, aa
    real(8)    :: psi2
    
    
    pop = 0.d0
    pop(1:noa) = 1.d0
    pop(noa+nva+1:noa+nva+nob) = 1.d0
    

    do ia = 2, nstates
       ii   = hole(ia) 
       aa   = part(ia)
       psi2 = dconjg(psi_det(ia))*psi_det(ia)
       if ( ii.lt.0 ) then
          i = -ii
          pop(i) = pop(i) - psi2
       else
          i = ii + noa + nva
          pop(i) = pop(i) - psi2
       end if
       if ( aa.lt.0 ) then
          a = -aa + noa
          pop(a) = pop(a) + psi2
       else
          a = aa + noa + nva + nob
          pop(a) = pop(a) + psi2
       end if
    end do

    
  end subroutine pop_soc
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_IP
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ip( hole,noa,nob,norb,nstates,nva,nvb,part,pop,psi_det )
    
    implicit none
    
    integer(8), intent(in) :: noa, nob, norb, nstates, nva, nvb
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    real(8), intent(inout) :: pop(norb)
    complex(8), intent(in) :: psi_det(nstates)

    integer(8) :: ia, ii, aa, xx, i, a, x, j
    integer(8) :: nrorb
    real(8)    :: psi2
    

    pop   = 0.d0
    nrorb = noa + nva
    

    pop(1:noa) = 1.d0
    pop( (nrorb+1):(nrorb+nob) ) = 1.d0

    do ia=1, nstates
       
       xx   = hole(ia,1) 
       ii   = hole(ia,2)
       aa   = part(ia)

       x = xx + nrorb       
       i = ii + nrorb
       a = aa + nrorb + nob 
       if ( xx.lt.0 ) x = -xx
       if ( ii.lt.0 ) i = -ii
       if ( aa.lt.0 ) a = -aa + noa
       
       psi2 = dcmplx(psi_det(ia)) * psi_det(ia)
       
       S_or_D : if ( ii.eq.0 ) then
          
          pop(x) = pop(x) - psi2
          
       else

          pop(x) = pop(x) - psi2
          pop(i) = pop(i) - psi2
          pop(a) = pop(a) + psi2
          
       end if S_OR_D
    
    end do


  end subroutine pop_ip
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_UION
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ion( noa,nva,norb,nstates,psi_det,vabsmoa,pop,nob,nvb,vabsmob)
    
    implicit none
    
    integer(8), intent(in) :: noa, nva, norb, nstates
    real(8),    intent(in) :: vabsmoa(noa+nva,noa+nva)
    complex(8), intent(in) :: psi_det(nstates)
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: vabsmob(nob+nvb,nob+nvb)
    real(8),    intent(inout) :: pop(norb)
    
    integer(8) :: i, a, b, ia, ib
    real(8)    :: psi2, rdum
    complex(8) :: psi0
    
    
    pop = 0.d0
    
    ! <C> diagonal <ia|V|ia> alpha
    do i=1, noa
       rdum = 0.d0
       do a=1, nva
          ia = (i-1)*nva + a + 1
          psi0 = dconjg(psi_det(ia))
          rdum = rdum + dconjg(psi0)*psi0 * vabsmoa(noa+a,noa+a)
          do b=(a+1), nva
             ib = (i-1)*nva + b + 1 
             rdum = rdum + 2.d0*real(psi0*psi_det(ib)) * vabsmoa(noa+a,noa+b) 
          end do
       end do
       pop(i) = rdum
    end do
    

    ! <C> diagonal <ia|V|ia> beta
    if ( present(nob) ) then
       do i=1, nob
          rdum = 0.d0
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             psi0 = dconjg(psi_det(ia))
             rdum = rdum + dconjg(psi0)*psi0 * vabsmob(nob+a,nob+a)
             do b=(a+1), nvb
                ib = (i-1)*nvb + b + noa*nva + 1
                rdum = rdum + 2.d0*real(psi0*psi_det(ib)) * vabsmob(nob+a,nob+b)
             end do
          end do
          pop(noa+nva+i) = rdum
       end do
    end if
    
    
  end subroutine pop_ion
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE POP_ION2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine pop_ion2( noa,nva,norb,nstates,psi_det,vabsmoa,pop,nob,nvb,vabsmob )
    
    implicit none

    integer(8), intent(in) :: noa,nva,norb,nstates
    real(8)   , intent(in) :: vabsmoa(noa+nva,noa+nva)
    complex(8), intent(in) :: psi_det(nstates)
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: vabsmob(nob+nvb,nob+nvb)
    real(8), intent(inout) :: pop(norb)
    
    integer(8) :: i,j,k,a,ia,b,jb
    real(8) :: tmp(norb,norb), dens(norb,norb), const
    complex(8) :: psi2, psi0
    
    
    dens = 0.d0 ; pop = 0.d0
    tmp(1:noa+nva,1:noa+nva) = vabsmoa(:,:)

    
    if ( present(nob) ) then
       tmp( noa+nva+1:norb, noa+nva+1:norb ) = vabsmob(:,:)
       const = 1.d0 
       do i=1, noa
          dens(i,i) = 1.d0
       end do
       do i=1, nob
          dens(noa+nva+i,noa+nva+i) = 1.d0
       end do
    else
       const = dsqrt(2.d0)
       do i=1, noa
          dens(i,i) = 2.d0
       end do
    end if
    

    ! <C> get one-electron density matrix < 0 | ia >_(N-1) alpha
    psi0 = dconjg(psi_det(1))
    do i=1, noa
       do a=1, nva
          ia = (i-1)*nva + a + 1 
          dens(i,noa+a) = const * real( psi0 * psi_det(ia) )
          dens(noa+a,i) = const * real( psi0 * psi_det(ia) )
       end do
    end do
    
    ! <C> < ia | jb > alpha
    do ia=1, noa*nva
       i = (ia-1)/nva + 1 
       a = ia - (i-1)*nva + noa
       psi2 = dconjg(psi_det(ia+1))
       do jb=1, noa*nva
          j = (jb-1)/nva + 1 
          b = jb - (j-1)*nva + noa
          if( i.eq.j ) dens(b,a) = dens(b,a) + real(psi2*psi_det(jb+1))
          if( a.eq.b ) dens(j,i) = dens(j,i) - real(psi2*psi_det(jb+1))
       end do
    end do

    
    ! <C> < ia | jb > beta
    if ( present(nob) ) then
       ! <C> get one-electron density matrix < 0 | ia >_(N-1) beta
       psi0 = dconjg(psi_det(1))
       do i=1, nob
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             dens(noa+nva+i,noa+nva+nob+a) = real( psi0 * psi_det(ia) )
             dens(noa+nva+nob+a,noa+nva+i) = real( psi0 * psi_det(ia) )
          end do
       end do
       ! <C> < ia | jb > beta
       do ia=1, nob*nvb
          i = (ia-1)/nvb + 1 
          a = ia - (i-1)*nvb + noa + nva + nob
          i = i + noa + nva
          psi2 = dconjg(psi_det(ia+1+noa*nva))
          do jb=1, nob*nvb
             j = (jb-1)/nvb + 1 
             b = jb - (j-1)*nvb + noa + nva + nob
             j = j + noa + nva
             if( i.eq.j ) dens(b,a) = dens(b,a) + real(psi2*psi_det(jb+1+noa*nva))
             if( a.eq.b ) dens(j,i) = dens(j,i) - real(psi2*psi_det(jb+1+noa*nva))
          end do
       end do
    end if
    

    tmp = matmul( dens, tmp ) 
    
    do i=1, norb
       pop(i) = tmp(i,i) 
    end do
    

  end subroutine pop_ion2
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_DYSON
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_dyson( noa,nva,norb,torbitals,nstuse,nstates,nbasis,cis_vec,cmo_a,psi,psi0,dirstr,emaxstr,nob,nvb,cmo_b )
    
    ! <C> use at your risk.  need to be double-checked

    implicit none
    
    integer(8),  intent(in) :: noa, nva, norb, torbitals, nstates, nstuse, nbasis
    real(8),     intent(in) :: cis_vec(nstates*nstates), cmo_a(nbasis,torbitals)
    complex(8),  intent(in) :: psi(nstuse), psi0(nstuse)
    character(4),intent(in) :: dirstr, emaxstr
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: cmo_b(nbasis,torbitals)
    
    integer(8) :: i,j,k,a,ia,b,jb, ibasis, jbasis, iorb, jorb, i0, den0
    real(8)    :: const, c, norm, norm0, dens_mo(torbitals,torbitals), dens_ao(nbasis*(nbasis+1)/2)
    complex(8) :: psi_det(nstates), psi_det0(nstates), psi2
    character(100) :: outfile
    
    
    norm0 = 1.d0
    norm = dot_product( psi,psi ) ; norm = dsqrt(norm)
    call get_psid(nstuse,nstates,cis_vec,norm0,psi0,psi_det0)
    call get_psid(nstuse,nstates,cis_vec,norm,psi,psi_det)
    
   
    dens_mo = 0.d0 

    den0 = 2.d0 ; const = dsqrt(2.d0)
    if ( present(nob) ) then
       den0 = 1.d0
       const = 1.d0
    end if

    i0 = torbitals - norb/2

    ! <C> ALPHA DENSITY 
    do i=1, i0
       dens_mo(i,i) = den0
    end do
    do i=1, noa
       dens_mo(i+i0,i+i0) = den0
    end do
    
    ! <C>  < 0 | ia >_(N-1) alpha
    psi2 = dconjg(psi_det0(1))
    do i=1, noa
       do a=1, nva
          ia = (i-1)*nva + a + 1 
          dens_mo(i+i0,noa+a+i0) = const * real( psi2 * psi_det(ia) )
          dens_mo(noa+a+i0,i+i0) = const * real( psi2 * psi_det(ia) )
       end do
    end do
    
    ! <C> < ia | jb > alpha
    do ia=1, noa*nva
       i = (ia-1)/nva + 1 
       a = ia - (i-1)*nva + noa + i0
       i = i + i0
       psi2 = dconjg(psi_det0(ia+1))
       do jb=1, noa*nva
          j = (jb-1)/nva + 1 
          b = jb - (j-1)*nva + noa + i0
          j = j + i0
          if( i.eq.j ) dens_mo(b,a) = dens_mo(b,a) + real(psi2*psi_det(jb+1))
          if( a.eq.b ) dens_mo(j,i) = dens_mo(j,i) - real(psi2*psi_det(jb+1))
       end do
    end do
    
    dens_ao = 0.d0 ; k = 0
    ! <C> get rho_ij to rho_uv
    do ibasis=1, nbasis
       do jbasis=1, ibasis
          k = k + 1 ; c = 0.d0
          do iorb=1, torbitals
             do jorb=1, torbitals
                c = c + cmo_a(ibasis,iorb) * cmo_a(jbasis,jorb) * dens_mo(jorb,iorb)
             end do
          end do
          if ( abs(c).lt.1.d-10 ) c = 0.d0
          dens_ao(k) = c
       end do
    end do


    
    ! <C> BETA DENSITY
    beta_den : if( present(nob) ) then

       dens_mo = 0.d0
       do i=1, i0
          dens_mo(i,i) = den0
       end do
       do i=1, nob
          dens_mo(i+i0,i+i0) = den0
       end do
       
       ! <C> get one-electron dens_moity matrix < 0 | ia >_(N-1) beta
       psi2 = dconjg(psi_det0(1))
       do i=1, nob
          do a=1, nvb
             ia = (i-1)*nvb + a + noa*nva + 1 
             dens_mo(i+i0,nob+a+i0) = real( psi2 * psi_det(ia) )
             dens_mo(nob+a+i0,i+i0) = real( psi2 * psi_det(ia) )
          end do
       end do
       
       ! <C> < ia | jb > beta
       do ia=1, nob*nvb
          i = (ia-1)/nvb + 1 
          a = ia - (i-1)*nva + nob + i0
          i = i + i0
          psi2 = dconjg(psi_det(ia+1+noa*nva))
          do jb=1, nob*nvb
             j = (jb-1)/nvb + 1 
             b = jb - (j-1)*nvb + nob + i0
             j = j + i0
             if( i.eq.j ) dens_mo(b,a) = dens_mo(b,a) + real(psi2*psi_det(jb+1+noa*nva))
             if( a.eq.b ) dens_mo(j,i) = dens_mo(j,i) - real(psi2*psi_det(jb+1+noa*nva))
          end do
       end do
       
       k = 0
       ! <C> get rho_ij to rho_uv
       do ibasis=1, nbasis
          do jbasis=1, ibasis
             k = k + 1 ; c = 0.d0
             do iorb=1, torbitals
                do jorb=1, torbitals
                   c = c + cmo_b(ibasis,iorb) * cmo_b(jbasis,jorb) * dens_mo(jorb,iorb)
                end do
             end do
             if ( abs(c).lt.1.d-10 ) c = 0.d0
             dens_ao(k) = dens_ao(k) + c
          end do
       end do
       
    end if beta_den

    outfile='DYSON'//'-e'//trim(emaxstr)//'-d'//trim(dirstr)//'.out'
    open( unit=434,file=trim(outfile) )
    write(434,"(5es16.8)") ( dens_ao(k), k=1, nbasis*(nbasis+1)/2 )
    close(434)


    
  end subroutine get_dyson
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_MODENSITY
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_modensity( noa,nva,norb,nstates,hole,part,psi_det,norm,dens,nob,nvb )

    implicit none 
    integer(8), intent(in) :: noa, nva, norb, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: norm
    complex(8), intent(in) :: psi_det(nstates)

    real(8), intent(inout) :: dens(norb,norb)
    integer(8), optional, intent(in) :: nob, nvb
    
    integer(8) :: ia, jb, i, a, j, b, ii, aa, jj, bb
    real(8)    :: const
    complex(8) :: c0, cdum

    
    const = dsqrt(2.d0)
    if( present(nob) ) const = 1.d0

    dens = 0.d0
    do i=1, noa
       dens(i,i) = norm * norm
    end do

    if ( present(nob) ) then
       do i = (noa+nva+1) , (noa+nva+nob )
          dens(i,i) =  norm * norm
       end do
    end if

    
    ! <C> 0,ia
    c0 = dconjg( psi_det(1) ) 
    do ia = 2, nstates
       
       ii = hole(ia)
       aa = part(ia)
       
       i = -ii
       a = -aa + noa
       if ( ii.gt.0 ) i = ii + noa + nva
       if ( aa.gt.0 ) a = aa + noa + nva + nob
       
       cdum = c0 * psi_det(ia)
       dens(a,i) = abs( cdum ) **2
       dens(i,a) = abs( cdum ) **2
       
    end do
    
    
    do ia = 2, nstates

       ii = hole(ia)
       aa = part(ia)
       i = -ii
       a = -aa + noa
       if ( ii.gt.0 ) i = ii + noa + nva
       if ( aa.gt.0 ) a = aa + noa + nva + nob

       c0 = dconjg(psi_det(ia))
       do jb = 2, nstates
          
          jj = hole(jb)
          bb = part(jb)
          j = -jj
          b = -bb + noa 
          if ( jj.gt.0 ) j = jj + noa + nva
          if ( bb.gt.0 ) b = bb + noa + nva  + nob
          
          if ( i.eq.j ) dens(a,b) = dens(a,b) + abs(c0 * psi_det(jb))**2
          if ( a.eq.b ) dens(j,i) = dens(j,i) - abs(c0 * psi_det(jb))**2
          
       end do
    end do


  end subroutine get_modensity
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FOURIER TRANSFORM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine fourier_transform( ntime, dt, corrt, outfile )
  
    implicit none

    real(8), parameter :: pi = 2.d0*acos(0.d0)      ! pi
    real(8), parameter :: au2eV = 27.21138602d0     ! hartrees to eV
    real(8), parameter :: au2fs = 2.418884326509d-2 ! au time to fs 
    
    integer(8), intent(in) :: ntime
    real(8),    intent(in) :: dt
    complex(8), intent(in) :: corrt(0:ntime)
    character(100) , intent(in) :: outfile
    

    ! <C> blackman window
    real(8), parameter :: a0 = 7938.d0/18608.d0
    real(8), parameter :: a1 = 9240.d0/18608.d0
    real(8), parameter :: a2 = 1430.d0/18608.d0

    integer(8) :: it, iw, i
    real(8)    :: pic, window(0:ntime), dw
    complex(8) :: corrw(0:ntime), c, expk0, expk
    

    corrw = dcmplx(0.d0,0.d0)
    dw = 2.d0 * pi / dt / dble(ntime)
    pic   = 2.d0 * pi / dble(ntime)

    ! <C> assign window function
    !do it=0, ntime
    !   window(it) = a0 - a1*cos(pic*dble(it)) + a2*cos(2.d0*pic*dble(it))
    !end do

    
    do iw=0, int(ntime/2)
       c = corrt(0)
       !expk0 = exp( eye * dble(iw) * pic )
       expk0 = dcmplx( cos(-dble(iw)*pic), sin(-dble(iw)*pic) )
       expk  = dcmplx( 1.d0,0.d0 )
       do it=1, ntime
          expk = expk * expk0
          c = c + corrt(it) * expk
       end do
       corrw(iw) = c
    end do

    open( unit=100,file=trim(outfile) )
    write(100,"( 2a20,'|',2a20 )") 'omega(eV)', 'corr(w)', 'time(fs)' ,'corr(t)'
    do i=0, ntime
       write(100,"(2f20.10,'|',2f20.10)")  dw*dble(i)*au2eV, abs( corrw(i) ), dt*dble(i)*au2fs, abs( corrt(i) )
    end do
    close(100)
    
    
  end subroutine fourier_transform
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITE_DATA
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine write_data(nstep,outstep,norb,ndata,dt,noa,nva,nob,nvb,normdata,ratedata,&
       popdata,iondata,mudata,normfile,popfile,ionfile,dipolefile,tdx,tdy,tdz,        &
       iout,writenorm,writepop,writeion,writemu,unrestricted)
    

    implicit none

    real(8), parameter :: au2fs = 2.418884326509d-2

    integer(8), intent(in) :: nstep, outstep, norb, ndata, noa, nva, nob, nvb
    integer(8), intent(in) :: iout
    real(8),    intent(in) :: dt, tdx, tdy, tdz
    real(8),    intent(in) :: normdata(ndata), ratedata(ndata)
    real(8),    intent(in) :: popdata(norb,ndata), iondata(norb,ndata), mudata(3,ndata)
    character(60), intent(in) :: normfile, popfile, ionfile, dipolefile
    logical, intent(in) :: writenorm, writepop, writeion, writemu, unrestricted

    integer(8) :: idata, j
    
    
    if( writenorm ) then
       open( 18, file=trim(normfile) )
       write(18,"( 2(1x,a9),2(1x,a12) )") 'time(au)','time(fs)','norm', 'rate(fs)'
       write(18,"( 2(1x,f9.3),2(1x,f13.10) )") 0.d0, 0.d0, 1.d0, 0.d0
       do idata = 1, ndata
          write(18,"( 2(1x,f9.3),2(1x,f13.10) )")   &
               dble(outstep)*dble(idata)*dt,       &
               dble(outstep)*dble(idata)*dt*au2fs, &
               normdata(idata), ratedata(idata)/au2fs
       end do
       close(18)
       write(iout,100) 'norm written to                         ',trim(normfile)
    end if
    

    if ( writepop ) then
       open( 19, file=trim(popfile) )
       do idata = 1, ndata
          write(19,"(' TIME(fs)= ',f9.2)") (outstep)*dble(idata)*dt*au2fs
          if ( unrestricted ) then
             write(19,101) ( popdata(j,idata), j=1,noa )
             write(19,300) ( popdata(j,idata), j=noa+1,noa+nva )
             write(19,200) ( popdata(j,idata), j=noa+nva+1,noa+nva+nob )
             write(19,400) ( popdata(j,idata), j=noa+nva+nob+1,noa+nva+nob+nvb )
          else
             write(19,500) ( popdata(j,idata), j=1,noa )
             write(19,600) ( popdata(j,idata), j=noa+1,noa+nva )
          end if
       end do
       close(19)
       write(iout,100) 'MO populations written to               ',trim(popfile)
    end if

101 format( '   occ_a:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
200 format( '   occ_b:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
300 format( '   vir_a:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
400 format( '   vir_b:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
500 format( '     occ:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )
600 format( '     vir:', 14(1x,f8.5) / 100( 15(1x,f8.5) /) )

    if ( writeion ) then
       open( 20,file=trim(ionfile) )
       do idata = 1, ndata
          write(20,"(' TIME(fs)=',f9.2,5x,'NORM=',f13.10,5x,'RATE(1/fs)=',f13.10,5x,'CHECK(1/fs)=',f13.10)") &
               (outstep)*dble(idata)*dt*au2fs, normdata(idata), ratedata(idata)/au2fs, sum(iondata(:,idata))/au2fs
          if ( unrestricted ) then
             write(20,101) ( iondata(j,idata)/au2fs, j=1,noa )
             write(20,300) ( iondata(j,idata)/au2fs, j=noa+1,noa+nva )
             write(20,200) ( iondata(j,idata)/au2fs, j=noa+nva+1,noa+nva+nob )
             write(20,400) ( iondata(j,idata)/au2fs, j=noa+nva+nob+1,noa+nva+nob+nvb )
          else
             write(20,500) ( iondata(j,idata)/au2fs, j=1,noa )
             write(20,600) ( iondata(j,idata)/au2fs, j=noa+1,noa+nva )
          end if
       end do
       close(20)
       write(iout,100) 'Instantaneous rate & ion pop written to ',trim(ionfile)
       write(iout,'(A)') '                                integrated MO populations for the ion'
    end if


    if ( writemu ) then
       open( 28,file=trim(dipolefile) )
       write(28,"(a9,3(1x,a12))") 'time(fs)','mu_x (au)','mu_y (au)','mu_z (au)'
       write(28,"(f9.2,3(1x,f13.10))") 0.d0, tdx, tdy, tdz
       do idata = 1, ndata
          write(28,"(f9.2,3(1x,f13.10))")           &
               dble(outstep)*dble(idata)*dt*au2fs, &
               mudata(1,idata), mudata(2,idata), mudata(3,idata)
       end do
       close(28)
       write(iout,100) 'Dipole moments written to               ',trim(dipolefile)
    end if
    

100 format(12x,a39,a18)    


    flush(iout)


  end subroutine write_data
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end module analyze_psi
