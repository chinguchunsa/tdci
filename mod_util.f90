module util

  
  implicit none
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DET00
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine det00(noa,nva,mo_matAA,myval,nob,nvb,mo_matBB)

    ! <C> get < 0 | 1e_operator | 0 > matrix element

    implicit none    
    
    real(8), intent(inout) :: myval
    integer(8), intent(in) :: noa, nva
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva)
    ! <C> if unrestricted 
    integer(8), optional, intent(in) :: nob, nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)

    
    integer(8) :: i 
    
    
    myval = 0.d0
    do i=1, noa 
       myval = myval + mo_matAA(i,i)
    end do
    
    if ( present(nob) ) then
       do i=1, nob
          myval = myval + mo_matBB(i,i)
       end do
    end if
       
        
  end subroutine det00
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE AO2MO
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine ao2mo(nbasis,nrorb,ao_mat,mo_mat,rotmat)  !nrorbB,mo_matBB,rotmatBB)
    

    ! <C> Transform one electron matrix from AO (lower triangle) to MO (full) 
    ! <C> mo_mat = [rotmat] [ao_mat] [rotmat]^T

    implicit none
    
    integer(8), intent(in) :: nbasis, nrorb
    real(8),    intent(in) :: ao_mat(nbasis*nbasis)
    real(8),    intent(in) :: rotmat(nbasis,nrorb)
    real(8), intent(inout) :: mo_mat(nrorb,nrorb)
    !integer(8), optional, intent(in) :: nrorbB
    !real(8), optional, intent(in)    :: rotmatBB(nbasis,nrorbB)
    !real(8), optional, intent(inout) :: mo_matBB(nrorbB,nrorbB)
    
    
    integer(8) :: iao, jao, kmo, lmo, ij, i
    real(8)    :: aodum, rotmatdum, rdum
    real(8)    :: tmp(nrorb,nbasis)
    
    
    mo_mat = 0.d0  
    tmp    = 0.d0
    
    ! <C> ao is stored as a 1D array, consisting of the lower half of the square matrix
    do iao=1, nbasis    
       do jao=1, nbasis 
          
          ij  = iao*(iao-1)/2 + jao
          if (jao.gt.iao) ij = jao*(jao-1)/2 + iao

          aodum = ao_mat(ij)
          tmp(:,jao) = tmp(:,jao) + rotmat(iao,:) * aodum
          
       end do
    end do
    
    mo_mat = matmul( tmp,rotmat )
    
    
  end subroutine ao2mo
  !: ---------------------------!
  !: SUBROUTINE AO2MO_COMPLEX   !
  !: ---------------------------!
  subroutine ao2mo_complex(nbasis,nrorb,ao_mat,mo_mat,rotmatAA,rotmatBB)  
    

    !: used for SOC ao-->mo transformation
    
    implicit none
    
    integer(8), intent(in) :: nbasis, nrorb
    complex(8), intent(in) :: ao_mat(nbasis,nbasis)
    real(8),    intent(in) :: rotmatAA(nbasis,nrorb), rotmatBB(nbasis,nrorb)
    complex(8), intent(inout) :: mo_mat(nrorb,nrorb)    
    
    integer(8) :: iao, jao, kmo, lmo, ij, i
    complex(8) :: aodum, cdum
    complex(8) :: tmp(nrorb,nbasis)
    
    
    mo_mat = dcmplx( 0.d0, 0.d0 )
    tmp    = dcmplx( 0.d0, 0.d0 )
    
    ! <C> ao is stored as a 1D array, consisting of the lower half of the square matrix
    do iao=1, nbasis    
       do jao=1, nbasis 
          aodum = ao_mat(iao,jao)
          do kmo=1, nrorb
             tmp(kmo,jao) = tmp(kmo,jao) + rotmatAA(iao,kmo)*aodum
          end do
       end do
    end do
    
    do jao=1, nbasis
       do lmo=1, nrorb
          do kmo=1, nrorb
             mo_mat(kmo,lmo) = mo_mat(kmo,lmo) + tmp(kmo,jao)*rotmatBB(jao,lmo)
          end do
       end do
    end do

    
  end subroutine ao2mo_complex
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE FORM1H
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine form1det(noa,nva,nstates,i2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CIS states          !
    ! <C> Xmat(I,A,0,0) = X(I,A)Sqrt(2)  FOR UNRESTRICTED NO SQRT(2)!
    ! <C> Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)         !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
    
    real(8), intent(inout) :: i2a_mat(nstates*nstates)

    integer(8), optional, intent(in) :: nob,nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)


    real(8) :: rdum, sqrt2
    integer(8) :: ii, jj, aa, bb, i, j, a, b
    integer(8) :: k, ia2, ia, jb, kk
    

    sqrt2  = dsqrt(2.d0)
    if( present(nob) ) sqrt2 = 1.d0    
    

    i2a_mat    = 0.d0
    i2a_mat(1) = mo00    

    ! <C> < ia | h(1) | 0 > == 1/sqrt(2) { <ia|h(1)|0> + <IA|h(1)|0> } = 2/sqrt(2) <ia|h(1)|0>
    ! <C> < ia | h(1) | 0 > == < a | h(1) | i >
    ia0 : do ia = 2, nstates

       ia2 = (ia-1)*nstates + 1
       ii  = hole(ia)  ;  i = abs(ii)
       aa  = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       
       if ( ii.lt.0 ) then
          i2a_mat(ia)  = sqrt2 * mo_matAA(i,a)
          i2a_mat(ia2) = sqrt2 * mo_matAA(i,a)
       else
          i2a_mat(ia)  = sqrt2 * mo_matBB(i,a)
          i2a_mat(ia2) = sqrt2 * mo_matBB(i,a)
       end if

    end do ia0

    
    ia_jb : do ia = 2, nstates

       ii = hole(ia)  ;  i = abs(ii)
       aa = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       

       do jb = 2 , nstates

          jj = hole(jb)  ;  j = abs(jj)
          bb = part(jb)  ;  b = abs(bb) + noa
          if ( bb.gt.0 ) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          rdum = 0.d0

          if ( ii.eq.jj ) then
             if ( ii.lt.0 ) rdum = mo_matAA(a,b)
             if ( ii.gt.0 ) rdum = mo_matBB(a,b)
          end if

          if ( aa.eq.bb ) then
             if ( aa.lt.0 ) rdum = rdum - mo_matAA(i,j)
             if ( aa.gt.0 ) rdum = rdum - mo_matBB(i,j)
          end if

          if ( ia.eq.jb ) rdum = rdum + mo00
          i2a_mat(k) = rdum
          
       end do
    end do ia_jb



  end subroutine form1det
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ZFORM1H
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine Zform1det(noa,nva,nstates,Zi2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CIS states          !
    ! <C> Xmat(I,A,0,0) = X(I,A)Sqrt(2)  FOR UNRESTRICTED NO SQRT(2)!
    ! <C> Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)         !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00  

    integer(8), intent(in) :: nob,nvb
    real(8),    intent(in) :: mo_matBB(nob+nvb,nob+nvb)

    complex(8), intent(inout) :: Zi2a_mat(nstates*nstates)
    

    complex(8) :: cdum
    integer(8) :: ii, jj, aa, bb, i, j, a, b
    integer(8) :: k, ia2, ia, jb, kk

    
    Zi2a_mat    = dcmplx( 0.d0,0.d0 )
    Zi2a_mat(1) = dcmplx( mo00, 0.d0 )


    !: < ia | h(1) | 0 > == < a | h(1) | i >
    ia0 : do ia = 2, nstates

       ia2 = (ia-1)*nstates + 1
       ii  = hole(ia)  ;  i = abs(ii)
       aa  = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob
       
       if ( ii*aa.gt.0 ) then
          if ( ii.lt.0 ) then
             Zi2a_mat(ia)  = dcmplx( mo_matAA(i,a), 0.d0 )
             Zi2a_mat(ia2) = dcmplx( mo_matAA(i,a), 0.d0 )
          else
             Zi2a_mat(ia)  = dcmplx( mo_matBB(i,a), 0.d0 )
             Zi2a_mat(ia2) = dcmplx( mo_matBB(i,a), 0.d0 )
          end if
       end if
       
    end do ia0

    
    ia_jb : do ia = 2, nstates

       ii = hole(ia)  ;  i = abs(ii)
       aa = part(ia)  ;  a = abs(aa) + noa
       if ( aa.gt.0 ) a = abs(aa) + nob

       do jb = 2 , nstates
          
          jj = hole(jb)  ;  j = abs(jj)
          bb = part(jb)  ;  b = abs(bb) + noa
          if ( bb.gt.0 ) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          cdum = dcmplx( 0.d0,0.d0 )
          
          if ( ii.eq.jj ) then
             if ( aa*bb .gt. 0 ) then
                if ( aa.lt.0 ) cdum = dcmplx( mo_matAA(a,b),0.d0 )
                if ( aa.gt.0 ) cdum = dcmplx( mo_matBB(a,b),0.d0 )
             end if
          end if
          
          if ( aa.eq.bb ) then
             if ( ii*jj .gt. 0 ) then
                if ( ii.lt.0 ) cdum = cdum - dcmplx( mo_matAA(i,j),0.d0 )
                if ( ii.gt.0 ) cdum = cdum - dcmplx( mo_matBB(i,j),0.d0 )
             end if
          end if
          
          if ( ia.eq.jb ) cdum = cdum + dcmplx( mo00,0.d0 )
          zi2a_mat(k) = cdum
          
       end do
    end do ia_jb
    


  end subroutine Zform1det
  !==================================================================!
  !==================================================================!
  subroutine form1det_ip(noa,nva,nstates,i2a_mat,mo_matAA,mo00,hole,part,nob,nvb,mo_matBB)
    
    ! <C> ==========================================================!
    ! <C> Form one electron matrix elements for CISD-IP states      !
    ! <C> ==========================================================!

    implicit none
    integer(8), intent(in) :: noa, nva, nstates
    integer(8), intent(in) :: hole(nstates,2), part(nstates)
    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
    
    real(8), intent(inout) :: i2a_mat(nstates*nstates)

    integer(8), optional, intent(in) :: nob,nvb
    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)


    real(8) :: rdum
    integer(8) :: ii, jj, aa, bb, xx, yy, i, j, a, b, x, y
    integer(8) :: k, ia, jb


    i2a_mat = 0.d0
    
    do ia=1, nstates
       
       xx = hole(ia,1)   ;  x = abs(xx)
       ii = hole(ia,2)   ;  i = abs(ii)
       aa = part(ia)     ;  a = abs(aa) + noa 
       if (aa.gt.0) a = abs(aa) + nob


       do jb=1, nstates
          
          yy = hole(jb,1)  ;  y = abs(yy)
          jj = hole(jb,2)  ;  j = abs(jj)
          bb = part(jb)    ;  b = abs(bb) + noa 
          if (bb.gt.0) b = abs(bb) + nob
          
          k    = (ia-1)*nstates + jb
          rdum = 0.d0


          !: singles
          
          
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
             if ( xx*yy .gt. 0 ) then !: same sign
                !: < x |A| y >
                if ( xx.lt.0 ) rdum = - mo_matAA(y,x)
                !: < X |A| Y >
                if ( xx.gt.0 ) rdum = - mo_matBB(Y,X)
             end if
             go to 78
          end if SS


          !: singles doubles
          
          
          SD1: if ( yy.eq.xx .and. jj.eq.0 ) then             
             !: < i->a , x |A| x >  =  < a |A| i >
             if ( ii.lt.0 ) rdum = mo_matAA(a,i) 
             !: < I->A, x |A| x >  =  < A |A| I > 
             if ( ii.gt.0 ) rdum = mo_matBB(A,I) 
             go to 78
          end if SD1
          
          SD2: if ( yy.eq.xx .and. ii.eq.0 ) then
             !: < x |A| j->b, x >  =  < j |A| b >
             if ( jj.lt.0 ) rdum = mo_matAA(b,j)
             !: < x |A| J->B, x >  =  < J |A| B >
             if ( jj.gt.0 ) rdum = mo_matBB(b,j)
             go to 78
          end if SD2

          SD3 : if ( yy.eq.ii .and. jj.eq.0 ) then
             same_spin0 : if ( xx*aa .gt. 0 ) then
                !: < y->a, x |A| y >  =  - < a |A| x >
                if ( xx.lt.0 ) rdum = - mo_matAA(a,x)
                !: < Y->A, X |A| Y >  =  - < A |A| X >
                if ( xx.gt.0 ) rdum = - mo_matBB(A,X)
             end if same_spin0
             go to 78
          end if SD3
          
          SD4: if ( xx.eq.jj .and. ii.eq.0 ) then
             same_spin00 : if ( yy*bb .gt. 0 ) then
                !: < x |A| x->b, y >  =  - < y | A | b >
                if ( yy.lt.0 ) rdum = - mo_matAA(b,y)
                !: < X |A| X->B, Y >  =  - < Y |A| B >
                if ( yy.gt.0 ) rdum = - mo_matBB(B,Y)
             end if same_spin00
             go to 78
          end if SD4
          

          if ( ii*jj.eq.0 ) go to 78
          

          !: doubles
          
          
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             !: < i->a, x |A| j->a, x >  =  - < j |A| i > 
             if ( ii.lt.0 ) rdum = - mo_matAA(i,j)
             !: < I->A, x |A| J->A, x >  =  - < J |A| I > 
             if ( ii.gt.0 ) rdum = - mo_matBB(I,J)
          end if kdelta_xy_ab
          
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             !: < i->a, x |A| i->b, x >  =  < a |A| b > 
             if ( aa.lt.0 ) rdum = rdum + mo_matAA(a,b)
             !: < I->A, x |A| I->B, x >  =  < A |A| B > 
             if ( aa.gt.0 ) rdum = rdum + mo_matBB(A,B)
          end if kdelta_xy_ij
          
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             same_spin1 : if ( xx*yy .gt. 0 ) then
                !: < i->a, x |A| i->a, y >  =  - < y |A| x >
                if ( xx.lt.0 ) rdum = rdum - mo_matAA(x,y)
                !: < I->A, X |A| I->A, Y >  =  - < Y |A| X >
                if ( xx.gt.0 ) rdum = rdum - mo_matBB(X,Y)
             end if same_spin1
          end if kdelta_ij_ab
          
          kdelta_jx_ab : if ( jj.eq.xx .and. aa.eq.bb ) then
             same_spin2 : if ( ii*yy .gt. 0 ) then
                !: < i->a, x |A| x->a, y >  =  < y |A| i >
                if ( ii.lt.0 ) rdum = rdum + mo_matAA(i,y)
                !: < I->A, X |A| X->A, Y >  =  < Y |A| I >
                if ( ii.gt.0 ) rdum = rdum + mo_matBB(I,Y) 
             end if same_spin2
          end if kdelta_jx_ab
          
          kdelta_yi_ab : if ( yy.eq.ii .and. aa.eq.bb ) then
             same_spin3 : if ( xx*jj .gt. 0 ) then
                !: < y->a, x |A| j->a, y >  =  < j |A| x >
                if ( jj.lt.0 ) rdum = rdum + mo_matAA(j,x)
                !: < Y->A, X |A| J->A, Y >  =  < J |A| X >
                if ( jj.gt.0 ) rdum = rdum + mo_matBB(J,X)
             end if same_spin3
          end if kdelta_yi_ab
          
          kdelta_iy_xj : if ( ii.eq.yy .and. xx.eq.jj ) then
             same_spin4 : if ( aa*bb .gt. 0 ) then
                !: < y->a, x |A| x->b, y >  =  - < a |A| b >
                if ( aa.lt.0 ) rdum = rdum - mo_matAA(a,b)
                !: < Y->A, X |A| X->B, Y >  =  - < A |A| B >
                if ( aa.gt.0 ) rdum = rdum - mo_matBB(A,B)
             end if same_spin4
          end if kdelta_iy_xj
          
78        continue
          
          if( ia.eq.jb ) rdum = rdum + mo00
          
          i2a_mat(k) = rdum
          
       end do
    end do

    
  end subroutine form1det_ip
  !==================================================================!
  !==================================================================!
  subroutine form1cis(nstates,nstuse,mat,cis_vec)
    
    ! <C> ==========================================================!
    ! <C>  Transformation B^T A B, with result written back to A    !
    ! <C>  A(N*N), B(N*M), V(N), results written back to A(M*M)     !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in) :: nstates, nstuse
    real(8), intent(in)    :: cis_vec(nstates*nstates)
    real(8), intent(inout) :: mat(nstates*nstates)
    
    integer(8) :: istate, jstate, ii, jj, iuse, kuse, kk, m, n, juse, jjj
    real(8)    :: dum
    real(8)    :: tmpmat(nstates)
    
    
    ! <C>  A*B and store in A
    do istate=1, nstates       

       ! <C> store a row from A to tmpmat
       do jstate=1,nstates
          jj = (jstate-1)*nstates
          tmpmat(jstate) = mat(jj+istate)
       end do
       
       ! <C> mat*cis_vec
       do kuse=1, nstuse
          kk = (kuse-1)*nstates
          mat(istate+kk) = dot_product( tmpmat(:), cis_vec(kk+1:kk+nstates) )
       end do
       
    end do

    
    ! <C>  B^T*(A*B) and store in A
    do kuse=1, nstuse       

       kk = (kuse-1)*nstates
       tmpmat(:) = mat(kk+1:kk+nstates)

       do iuse=1, nstuse
          ii = (iuse-1)*nstates
          mat(iuse+kk) = dot_product( tmpmat(:), cis_vec(ii+1:ii+nstates) )
       end do

    end do
    

    ! <C> Pack down if nstuse < nstates
    if(nstuse.lt.nstates) then
       do juse=2, nstuse
          jj  = (juse-1)*nstuse
          jjj = (juse-1)*nstates          
          do iuse=1, nstuse
             mat(iuse+jj)=mat(iuse+jjj)
          end do
       end do
    end if


  end subroutine form1cis
  !==================================================================!
  !==================================================================!
  subroutine Zform1cis(nstates,nstuse,Zmat,Zcis_vec)
    
    ! <C> ==========================================================!
    ! <C>  Transformation B^T A B, with result written back to A    !
    ! <C>  A(N*N), B(N*M), V(N), results written back to A(M*M)     !
    ! <C> ==========================================================!
    
    implicit none
    
    integer(8), intent(in) :: nstates, nstuse
    complex(8), intent(in) :: Zcis_vec(nstates*nstates)
    complex(8), intent(inout) :: Zmat(nstates*nstates)
    
    integer(8) :: istate, jstate, ii, jj, iuse, kuse, kk, m, n, juse, jjj
    real(8)    :: rdum
    complex(8) :: cdum, Ztmpmat(nstates)
    
    
    ! <C>  A*B and store in A
    do istate=1, nstates       
       
       ! <C> store a row from A to tmpmat
       !: Zmat initially all real
       do jstate=1, nstates
          jj = (jstate-1)*nstates
          Ztmpmat(jstate) = Zmat(jj+istate)
       end do
       
       ! <C> mat*cis_vec
       do kuse=1, nstuse
          kk = (kuse-1)*nstates
          cdum = dcmplx( 0.d0,0.d0 )
          do jstate=1, nstates
             cdum = cdum + Ztmpmat(jstate)*Zcis_vec(kk+jstate)
          end do
          Zmat(kk+istate) = cdum
          !Zmat(istate+kk) = dot_product( Zcis_vec(kk+1:kk+nstates) , Ztmpmat(:) ) 
       end do
       
    end do

    
    ! <C>  B^T*(A*B) and store in A
    do kuse=1, nstuse       
       
       kk = (kuse-1)*nstates
       do jstate=1, nstates
          Ztmpmat(jstate) = Zmat(kk+jstate)
       end do
       
       !Ztmpmat(:) = Zmat(kk+1:kk+nstates)
       
       do iuse=1, nstuse
          ii = (iuse-1)*nstates
          cdum = dcmplx(0.d0,0.d0)
          do jstate=1, nstates
             cdum = cdum + dconjg( Zcis_vec(ii+jstate) ) * Ztmpmat(jstate)
          end do
          !Zmat(kk+iuse) = dot_product( Ztmpmat(:), Zcis_vec(ii+1:ii+nstates) ) 
          Zmat(kk+iuse) = cdum
       end do
       
    end do
    
    
    ! <C> Pack down if nstuse < nstates
    if(nstuse.lt.nstates) then
       do juse=2, nstuse
          jj  = (juse-1)*nstuse
          jjj = (juse-1)*nstates          
          do iuse=1, nstuse
             Zmat(iuse+jj) = Zmat(iuse+jjj)
          end do
       end do
    end if


  end subroutine Zform1cis
  !==================================================================!
  !==================================================================!
  subroutine exp_mat(nstuse,mat,expmat,hdt)

    ! <C> get exponential of Hermetian matrix
    ! <C> diag_mat = [R]^T [mat] [R]
    ! <C> Exp_mat  = [R] exp( [diag_mat] ) [R]^T 
    ! <C> Exp_mat  = [R] [exp(eig)] [R]^T
    
    
    implicit none
    
    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: mat(nstuse*nstuse), hdt
    real(8), intent(inout) :: expmat(nstuse*nstuse)
    
    integer(8) :: i,k,kk,j,jj
    real(8)    :: start, finish, mat_k, exp_eig
    
    integer(8) :: info, lwork
    real(8)    :: eig(nstuse), mat_cp(nstuse*nstuse)
    real(8), allocatable :: work(:)
    

    mat_cp = mat

    ! <C> diagonalize matrix
    ! <C> dysev(  'v'   ,  'u'  ,  N  , A , LDA , W ,WORK,LWORK,INFO)
    lwork= -1
    call dsyev('vectors','upper',nstuse,mat_cp,nstuse,eig,eig,lwork,info)
    lwork = max(3*nstuse,int(eig(1)))
    allocate( work(lwork) )
    call dsyev('vectors','upper',nstuse,mat_cp,nstuse,eig,work,lwork,info)
    deallocate( work )
    
    
    ! <C> Form the exponential of the matrix    
    expmat = 0.d0    
    do j=1, nstuse
       
       jj   = (j-1)*nstuse
       exp_eig = dexp( hdt*eig(j) )
       do k=1, nstuse
          kk = (k-1)*nstuse
          mat_k = mat_cp(jj+k) * exp_eig
          expmat(kk+1:kk+nstuse) = expmat(kk+1:kk+nstuse) + mat_k * mat_cp(jj+1:jj+nstuse)
       end do
       
    end do
    

  end subroutine exp_mat
  !==================================================================!
  !==================================================================!
  subroutine Zexp_mat(nstuse,Zmat,Zexpmat,hdt)

    ! <C> get exponential of Hermetian matrix
    ! <C> diag_mat = [R]^T [mat] [R]
    ! <C> Exp_mat  = [R] exp( [diag_mat] ) [R]^T 
    ! <C> Exp_mat  = [R] [exp(eig)] [R]^T
    
    
    implicit none
    
    integer(8), intent(in) :: nstuse
    real(8),    intent(in) :: hdt
    complex(8), intent(in) :: Zmat(nstuse*nstuse)
    complex(8), intent(inout) :: Zexpmat(nstuse*nstuse)
    
    integer(8) :: i,k,kk,j,jj
    real(8)    :: start, finish, mat_k, exp_eig
    
    integer(8) :: info, lwork
    real(8)    :: eig(nstuse)
    complex(8) :: mat_cp(nstuse*nstuse)
    real(8), allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    
    
    mat_cp = Zmat

    ! <C> diagonalize matrix
    lwork = 2*nstuse - 1 
    info = 10
    allocate( work(lwork) ) 
    allocate( rwork(3*nstuse-2) )
    !:   zheev(   'v'   ,  'u'  , N    ,   A  ,  LDA , W ,WORK,LWORK,RWORK,INFO)
    call zheev('vectors','upper',nstuse,mat_cp,nstuse,eig,work,lwork,rwork,info)
    deallocate( work, rwork )
    
    
    ! <C> Form the exponential of the matrix    
    Zexpmat = dcmplx( 0.d0,0.d0 )
    do j=1, nstuse
       
       jj   = (j-1)*nstuse
       exp_eig = dexp( hdt*eig(j) )
       do k=1, nstuse
          kk = (k-1)*nstuse
          !mat_k = mat_cp(jj+k) * exp_eig
          !Zexpmat(kk+1:kk+nstuse) = Zexpmat(kk+1:kk+nstuse) + mat_k * dconjg(mat_cp(jj+1:jj+nstuse))
          mat_k = dconjg(mat_cp(jj+k)) * exp_eig
          Zexpmat(kk+1:kk+nstuse) = Zexpmat(kk+1:kk+nstuse) + mat_k * mat_cp(jj+1:jj+nstuse)
       end do
       
    end do
    

  end subroutine Zexp_mat
  !==================================================================!
  !==================================================================!
  subroutine matmult(A,B,nsize)
   
    
    ! <C> commented section does
    ! <C> Matrix Multiplication A*B=C
    ! <C> C is written back to A       

    
    integer(8), intent(in) :: nsize
    real(8),    intent(in) :: B(nsize*nsize) !B(nsize,nsize)
    real(8), intent(inout) :: A(nsize*nsize) !A(nsize,nsize)

    integer(8) :: i,ii,j,jj,k,kk
    real(8) :: rdum, scratch(nsize)    !scratch(nsize,nsize)

    
    !scratch = matmul( A,B )
    !A = scratch

    do i=1, nsize
       ii=(i-1)*nsize
       do j=1, nsize
          scratch(j) = A(ii+j)
       end do
       do k=1, nsize
          kk = (k-1)*nsize
          rdum = 0.d0
          do j=1, nsize
             rdum = rdum+scratch(j)*B(j+kk)
          end do
          A(ii+k) = rdum
       end do
    end do


    do i=2, nsize
       ii=(i-1)*nsize
       do j=1, i
          jj = (j-1)*nsize
          rdum = A(i+jj)
          A(i+jj) = A(ii+j)
          A(ii+j) = rdum
       end do
    end do

    
  end subroutine matmult
  !==================================================================!
  !==================================================================!
  subroutine Zmatmult(A,B,nsize)
   
    
    ! <C> commented section does
    ! <C> Matrix Multiplication A*B=C
    ! <C> C is written back to A       

    
    integer(8), intent(in) :: nsize
    complex(8), intent(in) :: B(nsize*nsize) !B(nsize,nsize)
    complex(8), intent(inout) :: A(nsize*nsize) !A(nsize,nsize)

    integer(8) :: i,ii,j,jj,k,kk
    complex(8) :: cdum, scratch(nsize)    !scratch(nsize,nsize)

    
    !scratch = matmul( A,B )
    !A = scratch
    
    do i=1, nsize
       ii=(i-1)*nsize
       do j=1, nsize
          scratch(j) = A(ii+j)
       end do
       do k=1, nsize
          kk = (k-1)*nsize
          cdum = dcmplx( 0.d0,0.d0 )
          do j=1, nsize
             cdum = cdum+scratch(j)*B(j+kk)
          end do
          A(ii+k) = cdum
       end do
    end do


    do i=2, nsize
       ii=(i-1)*nsize
       do j=1, i
          jj = (j-1)*nsize
          cdum = A(i+jj)
          A(i+jj) = A(ii+j)
          A(ii+j) = cdum
       end do
    end do

    
  end subroutine Zmatmult
  !==================================================================!
  !==================================================================!
end module util
    
!  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!  ! SUBROUTINE FORM1H
!  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
!  subroutine form1det(noa,nva,nstates,i2a_mat,mo_matAA,mo00,nob,nvb,mo_matBB)
!    
!    ! <C> ==========================================================!
!    ! <C> Form one electron matrix elements for CIS states          !
!    ! <C> Xmat(I,A,0,0) = X(I,A)Sqrt(2)  FOR UNRESTRICTED NO SQRT(2)!
!    ! <C> Xmat(I,A,J,B) = delta(I,J)X(A,B)-delta(A,B)X(I,J)         !
!    ! <C> ==========================================================!
!
!    implicit none
!    integer(8), intent(in) :: noa, nva, nstates
!    real(8),    intent(in) :: mo_matAA(noa+nva,noa+nva), mo00
!    
!    real(8), intent(inout) :: i2a_mat(nstates*nstates)
!
!    integer(8), optional, intent(in) :: nob,nvb
!    real(8),    optional, intent(in) :: mo_matBB(nob+nvb,nob+nvb)
!
!
!    real(8) :: rdum, sqrt2
!    integer(8) :: ii, jj, aa, bb, i, j, a, b
!    integer(8) :: k, ia2, ia, jb, kk
!    
!
!    sqrt2  = dsqrt(2.d0)
!    if( present(nob) ) sqrt2 = 1.d0    
!    
!
!    i2a_mat    = 0.d0
!    i2a_mat(1) = mo00    
!
!    
!    !: <ia|h(1)|0> = 1/sqrt(2) { <ia|h(1)|0> + <IA|h(1)|0> } = 2/sqrt(2) <ia|h(1)|0> 
!    !: <ia|h(1)|0> = <a|h(1)|i> 
!    
!    ia0_alpha : do ia = 1, noanva
!       k  = ia + 1
!       kk = nstates * (k-1)
!       i  = (ia-1)/nva + 1
!       a  = ia - nva*(i-1) + noa
!       i2a_mat(kk+1) = sqrt2 * mo_matAA(i,a)
!       i2a_mat(k)    = i2a_mat(kk+1)
!    end do ia0_alpha
!
!    ! <C> < ia | h(1) | jb > alpha alpha block
!    iajb_alpha : do ia = 1, noanva
!       k  = ia + 1
!       kk = nstates*(k-1)
!       i  = (ia-1)/nva + 1
!       a  = ia - nva*(i-1) + noa
!       do jb = 1, noanva
!          l = jb + 1
!          j = (jb-1)/nva + 1
!          b = jb - nva*(j-1) + noa
!          rdum = 0.d0
!          if(i.eq.j) rdum = rdum + mo_matAA(a,b)
!          if(a.eq.b) rdum = rdum - mo_matAA(i,j)
!          i2a_mat(kk+l) = rdum
!       end do
!       i2a_mat(kk+k) = i2a_mat(kk+k) + mo00
!    end do iajb_alpha
!
!    unrestrictedpart : if( present(nob) ) then
!       
!       ia0_beta : do ia = 1, nobnvb
!          k  = ia + noanva + 1
!          kk = nstates * (k-1)
!          i  = (ia-1)/nvb + 1
!          a  = ia - nvb*(i-1) + nob
!          i2a_mat(kk+1) = mo_matBB(i,a)
!          i2a_mat(k)    = i2a_mat(kk+1)
!       end do ia0_beta
!       
!       iajb_beta : do ia = 1, nobnvb
!          k  = ia + noanva + 1
!          kk = nstates * (k-1)
!          i  = (ia-1)/nvb + 1
!          a  = ia - nvb*(i-1) + nob
!          do jb = 1, nobnvb
!             l = jb + noanva + 1
!             j = (jb-1)/nvb + 1
!             b = jb - nvb*(j-1) + nob
!             rdum = 0.d0
!             if(i.eq.j) rdum = rdum + mo_matBB(a,b)
!             if(a.eq.b) rdum = rdum - mo_matBB(i,j)
!             i2a_mat(kk+l) = rdum
!          end do
!          i2a_mat(kk+k) = i2a_mat(kk+k) + mo00
!       end do iajb_beta
!
!    end if unrestrictedpart
!
!
!
!  end subroutine form1det
!  !==================================================================!
!  !==================================================================!
