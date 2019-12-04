module get_ham0
  
  use global_variables
  
  implicit none

  integer(8) :: double_startAA, double_startBB, double_startAB, double_startBA


  !: contains subroutine get_cis
  !: contains subroutine get_cis_index
  !: contains subroutine get_ip_cisd
  !: contains subroutine get_ip_index
  !: contains subroutine diagonalize
  !: contains subroutine writeme_ham0( myroutine, option )
  !: contains subrotuine errors_getham0( myroutine, option, mystring )
  

contains
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_RCIS
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cis
    
    !: Form restricted or unrestricted  CIS Hamiltonian                                       
    !: ia = (i-1)*nva + a ; ib = (i-1)*nva + b
    !: jb = (j-1)*nvb + b ; ja = (j-1)*nvb + a     
    
    use read_integrals
    implicit none
    
    integer(8) :: istate, i, j, k, l, a, b
    integer(8) :: ia2, jb2, ia, ib, ja, jb
    real(8)    :: hmat(nstates,nstates)
    

    call write_header( 'get_cis','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'cis','form' ) 
    call get_cis_index
    
    !: temporarily store ham.  will transfer ham(rank2) to cis_vec(rank1)
    hmat   = 0.d0    
    
    if ( .not.unrestricted ) then    
       
       restricted : do ia = 2, nstates
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          ia2 = ia - 1 !: dijabAB index
          do jb = 2, nstates
             j = -hole_index(jb,1)
             b = -part_index(jb,1)   
             jb2 = jb - 1        !: dijabAB index
             ib  = (i-1)*nva + b !: diajbAA index
             ja  = (j-1)*nva + a !: diajbAA index
             hmat(jb,ia) = dijabAB(jb2,ia2) - diajbAA(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(noa+a) - orben(i)
       end do restricted
       
    else 
       
       aa : do ia = 2 , 1+noanva
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          do jb = 2 , 1+noanva
             j = -hole_index(jb,1)
             b = -part_index(jb,1)             
             ib = (i-1)*nva + b !: diajbAA index
             ja = (j-1)*nva + a !: diajbAA index
             hmat(jb,ia) = -diajbAA(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(noa+a) - orben(i)
       end do aa

       !: beta, beta block 
       bb : do ia = 2+noanva, nstates
          i = hole_index(ia,1)
          a = part_index(ia,1)
          do jb = 2+noanva, nstates
             j  = hole_index(jb,1)
             b  = part_index(jb,1)       
             ib = nvb*(i-1) + b !: diajbBB index
             ja = nvb*(j-1) + a !: diabBB index
             hmat(jb,ia) = -diajbBB(ib,ja)
          end do
          hmat(ia,ia) = hmat(ia,ia) + orben(a+nob+nrorb) - orben(i+nrorb)
       end do bb
       
       !: alpha,beta block ;; beta,alpha block
       ab : do ia = 2, 1+noanva
          i = -hole_index(ia,1)
          a = -part_index(ia,1)
          ia2 = (i-1)*nva + a
          do jb = 2+noanva, nstates
             j  = hole_index(jb,1)
             b  = part_index(jb,1)                          
             jb2 = (j-1)*nvb + b
             hmat(jb,ia) = dijabAB(jb2,ia2)
             hmat(ia,jb) = hmat(jb,ia)
          end do
       end do ab
       
    end if
    
    !: transfer to cis_vec
    istate=0
    do i=1, nstates
       do j=1, nstates
          istate = istate + 1 
          cis_vec( istate ) = hmat(j,i)
       end do
    end do

    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start
    
    
    call write_header( 'get_cis','get_ham0','leave' )

    
  end subroutine get_cis
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_RCIS_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_cis_index


    !: assign hole and particle excitation indices for each state
    implicit none
    
    integer(8) :: i, j, a, b, istate
    
    
    call writeme_ham0( 'cis_index', 'imhere' )

    !: HF reference state
    istate = 1
    hole_index(istate,1) = 0
    part_index(istate,1) = 0
    
    do i=1, noa
       do a=1, nva
          istate = istate + 1 
          hole_index(istate,1) = -i
          part_index(istate,1) = -a ! noa + a 
       end do
    end do    

    if ( unrestricted ) then
       do i=1, nob
          do a=1, nvb
             istate = istate + 1 
             hole_index(istate,1) = i ! noa+nva+i
             part_index(istate,1) = a ! noa+nva+nob+a
          end do
       end do
    end if
    
    call writeme_ham0( 'cis_index','imdone' )
    

  end subroutine get_cis_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_IP_CISD
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_ip_cisd

    !: get IP-CISD Hamiltonian

    use read_integrals
    implicit none

    integer(8) :: istate, xorb, iorb, aorb
    integer(8) :: xx, ii, aa, yy, jj, bb, x, y, i, a, j, b
    integer(8) :: ia, jb, ia2, jb2, yx, xb, yj, ji, ix, xi, jy, ab, ya, ib, ja, yb, xa
    real(8) :: sign, storeme
    
    
    call write_header( 'get_ip_cisd','get_ham0','enter' )
    
    call cpu_time(start)
    call writeme_ham0( 'ip', 'form' ) 
    call get_ip_index 
    
    ia : do ia=1, nstates
       
       xx = hole_index(ia,1) ; x = abs(xx)
       ii = hole_index(ia,2) ; i = abs(ii)
       aa = part_index(ia,1) ; a = abs(aa)
       
       !: orbital indices
       if ( xx.lt.0 ) xorb = -xx       ; if ( xx.gt.0 ) xorb = nrorb + xx
       if ( ii.lt.0 ) iorb = -ii       ; if ( ii.gt.0 ) iorb = nrorb + ii
       if ( aa.lt.0 ) aorb = -aa + noa ; if ( aa.gt.0 ) aorb = nrorb + nob + aa
       
       jb : do jb=ia, nstates
          
          yy = hole_index(jb,1) ; y = abs(yy)
          jj = hole_index(jb,2) ; j = abs(jj)
          bb = part_index(jb,1) ; b = abs(bb)
          
          storeme = 0.d0
          
          SS : if( ii.eq.0 .and. jj.eq.0 ) then
             if( xx.eq.yy ) then
                storeme = -orben(xorb)
             end if
             go to 78
          end if SS
          
          
          !: need to consider ii.ne.0 and jj.eq.0 if do loops go from jb=1
          SD_bb : if ( xx.gt.0 .and. yy.gt.0 ) then
             if( ii.eq.0 .and. jj.ne.0 ) then             
                a_or_b : if( jj.lt.0 ) then
                   !: <X|H|jb_Y> = -<Yj||Xb> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(jb2,Y<=X)
                   jb2 = (j-1)*nva + b
                   YX  = (Y-1)*nob - (Y-1)*(Y-2)/2 + (X-Y) + 1
                   if ( X.lt.Y ) YX = (X-1)*nob - (X-1)*(X-2)/2 + (Y-X) + 1
                   !:----------------------------!
                   storeme = - dijkaBA(jb2,YX)  
                else
                   !: <X|H|JB_Y> = -<YJ||XB> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(XB,Y<J)
                   XB = (X-1)*nvb + B
                   YJ = (Y-1)*nob - Y*(Y-1)/2 + (J-Y)    ; sign = 1.d0
                   if ( J.lt.Y ) then
                      YJ = (J-1)*nob - J*(J-1)/2 + (Y-J) ; sign = -sign
                   end if
                   !:-----------------------------------!
                   storeme = - sign * dijkaBB(XB,YJ) 
                end if a_or_b                
                go to 78                
             else if  ( jj.eq.0 .and. ii.ne.0 ) then
                a_or_b2 : if ( ii.lt.0 ) then
                   !: <Y|H|ia_X> = -<Xi||Ya> = -<Ij||Ka> : DijkaBA(ja,I<=K) = DijkaBA(ia, X<=Y)
                   ia2 = (i-1)*nva + a
                   YX  = (X-1)*nob - (X-1)*(X-2)/2 + (Y-X) + 1
                   if ( Y.lt.X ) YX = (Y-1)*nob - (Y-1)*(Y-2)/2 + (X-Y) + 1
                   !:----------------------------!
                   storeme = - dijkaBA(ia2,YX)
                else
                   !: <Y|H|IA_X> = -<XI||YA> = -<IJ||KA> : DijkaBB(KA,I<J) = DijkaBB(YA,X<I)
                   YA = (Y-1)*nvb + A
                   XI = (X-1)*nob - X*(X-1)/2 + (I-X)    ; sign=1.d0
                   if ( I.lt.X ) then
                      XI = (I-1)*nob - I*(I-1)/2 + (X-I) ; sign=-1.d0
                   end if
                   storeme = - sign * dijkaBB(YA,XI)
                end if a_or_b2                
                go to 78                
             end if
          end if SD_bb

          !SD_aa : if ( xx.lt.0 .and. yy.lt.0 ) then
             
          

          DD : if( ii.ne.0 .and. jj.ne.0 ) then
             
             kdelta_ab : if ( aa.eq.bb ) then
                if ( aa.lt.0 ) then
                   !: <jY||iX> = <iJ||kL> : DijklAB(J<=L,i<=k)=DijklAB(Y<=X,j<=i)
                   YX = (Y-1)*nob - (Y-1)*(Y-2)/2 + (X-Y) + 1
                   ji = (j-1)*noa - (j-1)*(j-2)/2 + (i-j) + 1
                   if ( X.lt.Y ) YX = (X-1)*nob - (X-1)*(X-2)/2 + (Y-X) + 1
                   if ( i.lt.j ) ji = (i-1)*noa - (i-1)*(i-2)/2 + (j-i) + 1
                   !:-----------------------------------!
                   storeme = storeme + dijklAB(YX,ji) 
                else
                   !: <JY||IX> = <IJ||KL> : DijklBB(k<l,i<j)=DijklBB(I<X,J<Y)
                   !: note, J<Y always in this setup 
                   IX = (I-1)*nob - I*(I-1)/2 + (X-I)
                   JY = (J-1)*nob - J*(J-1)/2 + (Y-J)    ; sign =  1.d0
                   if ( Y.lt.J ) then
                      JY = (Y-1)*nob - Y*(Y-1)/2 + (J-Y) ; sign = -1.d0
                   end if
                   if ( X.lt.I ) then
                      IX = (X-1)*nob - X*(X-1)/2 + (I-X) ; sign = -sign
                   end if
                   !:-----------------------------------------!
                   storeme = storeme + sign * dijklBB(IX,JY) 
                end if
             end if kdelta_ab

             kdelta_ij : if ( ii.eq.jj ) then
                if ( aa.lt.0 ) then
                   !: -<Ya||Xb> = <Ia||Jb> : DiajbBA(a<=b,I<=J)=DiajbBA(a<=b,Y<=X)
                   YX = (Y-1)*nob - (Y-1)*(Y-2)/2 + (X-Y) + 1
                   ab = (a-1)*nva - (a-1)*(a-2)/2 + (b-a) + 1
                   if ( X.lt.Y ) YX = (X-1)*nob - (X-1)*(X-2)/2 + (Y-X) + 1
                   if ( b.lt.a ) ab = (b-1)*nva - (b-1)*(b-2)/2 + (a-b) + 1
                   !:-----------------------------------!
                   storeme = storeme - diajbBA(ab,YX)  
                else
                   !: -<YA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,YA)
                   XB = (X-1)*nvb + B
                   YA = (Y-1)*nvb + A
                   !:-----------------------------------!
                   storeme = storeme - diajbBB(XB,YA) 
                end if
             end if kdelta_ij
             

             kdelta_xy : if ( XX.eq.YY ) then
                aa_block : if ( ii.lt.0 .and. jj.lt.0 ) then
                   !: -<ja||ib> = -<ia||jb> : DiajbAA(jb,ia) = DiajbAA(ib,ja)
                   ja = (j-1)*nva + a
                   ib = (i-1)*nva + b
                   !:-----------------------------------!
                   storeme = storeme - diajbAA(ib,ja)
                end if aa_block

                bb_block :if ( II.gt.0 .and. JJ.gt.0 ) then
                   !: -<JA||IB> = -<ia||jb> : DiajbBB(jb,ia) = DiajbBB(IB,JA)
                   JA = (J-1)*nvb + A
                   IB = (I-1)*nvb + B
                   !:----------------------------------!
                   storeme = storeme - diajbBB(IB,JA) 
                end if bb_block

                ab_block : if( ii.lt.0 .and. JJ.gt.0 ) then
                   !: <iJ||aB>  : DijabAB(JB,ia)
                   ia2 = (i-1)*nva + a
                   JB2 = (J-1)*nvb + B
                   !:--------------------------!
                   storeme = dijabAB(JB2,ia2) 
                end if ab_block

                ba_block : if( II.gt.0 .and. jj.lt.0 ) then
                   !: <Ij||Ab> = <jI||bA> : DijabAB(IA,jb)
                   IA2 = (I-1)*nvb + A
                   jb2 = (j-1)*nva + b
                   !:---------------------------!
                   storeme = dijabAB(IA2,jb2)  
                end if ba_block
             end if kdelta_xy

             
             kdelta_jx : if ( JJ.eq.XX ) then
                if ( ii.gt.0 ) then
                   !: <YA||IB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(IB,YA)
                   IB = (I-1)*nvb + B
                   YA = (Y-1)*nvb + A
                   !:----------------------------------!
                   storeme = storeme + diajbBB(IB,YA) 
                else
                   !: <Ya||iB> = -<iY||aB> : DijabAB(YB,ia)
                   YB  = (Y-1)*nvb + B
                   ia2 = (i-1)*nva + a
                   !:---------------------------!
                   storeme = -dijabAB(YB,ia2) 
                end if
             end if kdelta_jx
             
             kdelta_iy : if ( II.eq.YY ) then
                if ( JJ.gt.0 ) then
                   !: <JA||XB> = <IA||JB> : DiajbBB(JB,IA)=DiajbBB(XB,JA)
                   XB = (X-1)*nvb + B
                   JA = (J-1)*nvb + A
                   !:----------------------------------!
                   storeme = storeme + diajbBB(XB,JA) 
                else
                   !: <jX||bA> : DijabAB(XA,jb)
                   XA  = (X-1)*nvb + A
                   jb2 = (j-1)*nva + b
                   !:---------------------------!
                   storeme = -dijabAB(XA,jb2)  
                end if
             end if kdelta_iy

             diagonal : if ( jb.eq.ia ) then
                if ( aa.lt.0 ) storeme = storeme - orben(i) - orben(x) + orben(noa+a)
                if ( AA.gt.0 ) storeme = storeme - orben(nrorb+I) - orben(nrorb+X) + orben(nrorb+nob+A)
             end if diagonal
             
             !special_off : if ( ii.eq.yy .and. xx.eq.jj .and. aa.eq.bb ) then
             
          end if DD
          
78        continue

          istate = (ia-1)*nstates + jb
          cis_vec(istate) = storeme

          istate = (jb-1)*nstates + ia
          cis_vec(istate) = storeme
          
       end do jb
    end do ia
   
    call cpu_time(finish)
    write(iout,"(' Hamiltonian generation time:',f12.4,' seconds')") finish-start    

    call write_header( 'get_ip_cisd','get_ham0','leave' )
    

  end subroutine get_ip_cisd
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE GET_IP_INDEX
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine get_ip_index

    !: assign hole and particle excitation indices for each state
    implicit none

    integer(8) :: i, j, a, b, x, istate


    call writeme_ham0( 'ip_index','imhere' )
    
    part_index(:,:) = 0
    hole_index(:,:) = 0
    istate = 0

    !: only beta electrons will be ionized
    singles : do x=1, nob
       istate = istate + 1 
       hole_index(istate,1) = x
       hole_index(istate,2) = 0
       part_index(istate,1) = 0
    end do singles

    ab_doubles : do x=1, nob
       aa : do i=1, noa
          do a=1, nva
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = -i
             part_index(istate,1) = -a
          end do
       end do aa       
    end do ab_doubles

    bb_doubles : do x=1, nob
       bb : do i=x+1, nob
          do a=1, nvb
             istate = istate + 1 
             hole_index(istate,1) = x
             hole_index(istate,2) = i
             part_index(istate,1) = a
          end do
       end do bb
    end do bb_doubles
    
    call writeme_ham0( 'ip_index','imdone' )

    
  end subroutine get_ip_index
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE DIAGONALIZE
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine diagonalize

    !: call lapack dysev to get eigenstates and eigenvalues    
    
    implicit none
    
    integer(8) :: lscratch, info, lwork
    real(8), allocatable :: work(:)
    character(20) :: cinfo

    
    call write_header( 'diagonalize', 'get_ham0', 'enter' )    

    call cpu_time(start)


    !<C> if LWORK = -1, workspace query is assumed.  routine ONLY calculates the
    !<C> optimal size of WORK array, returns this value as 1rst entry of WORK array.
    lwork = -1
    !<C> dysev(   'v'   ,  'u'  , N     ,   A   ,  LDA  ,   W   , WORK  ,LWORK,INFO)
    call dsyev('vectors','upper',nstates,cis_vec,nstates,cis_eig,cis_eig,lwork,info)

    
    lwork = max( 3*nstates,int(cis_eig(1)) )!write(iout,"(' optimal WORK space required:  ',i0)") lwork
    info  = 10
    call writeme_ham0( 'diagonalize' , 'calling' )
    

    allocate( work(lwork) )
    !: dysev(  'v'   ,  'u'  , N     ,   A   ,  LDA  ,   W   ,WORK,LWORK,INFO)
    call dsyev('vectors','upper',nstates,cis_vec,nstates,cis_eig,work,lwork,info)
    deallocate( work )
    

    write( cinfo,'(i0)' ) info
    call cpu_time(finish)
    

    if ( info.eq.0 ) call errors_getham0( 'diagonalize','success', trim(cinfo) )
    if ( info.lt.0 ) call errors_getham0( 'diagonalize','illegal', trim(cinfo) )
    if ( info.gt.0 ) call errors_getham0( 'diagonalize','no_converge', trim(cinfo) )
    if ( info.eq.0 ) write(iout,"(' LAPACK dysev diagonalization time:',f12.4,' seconds')") finish-start

    call write_header( 'diagonalize','get_ham0','leave' )

    
  end subroutine diagonalize
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE WRITEME_HAM0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine writeme_ham0(myroutine,option)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option
    

    select case( trim(myroutine) )
    case( 'cis' ) 
       select case( trim(option) )
       case('form')    ; write(iout,"(' forming ',i0,' x ',i0,' CIS Hamiltonian')") nstates, nstates
       end select
    case( 'cis_index' )
       select case( trim(option) ) 
       case( 'imhere' ) ; write(iout,"(' assigning hole_index(nstates,1) and part_index(nstates,1) for CIS')")
       case( 'imdone' ) ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'ip' )
       select case( trim(option) )
       case('form')    ; write(iout,"(' forming ',i0,' x ',i0,' IP-CISD Hamiltonian')") nstates, nstates
       end select
    case( 'ip_index' )
       select case( trim(option) )
       case( 'imhere' ) ; write(iout,"(' assigning hole_index(nstates,2) and part_index(nstates,1) for IP-CISD')")
       case( 'imdone' ) ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'cisd_index' )
       select case( trim(option) )
       case('imhere') ; write(iout,"(' assigning hole_index(nstates,2) and part_index(nstates,2) for CISD')")
       case('imdone') ; write(iout,"(' finished assigning indices.  Negative indices for alpha electrons')")
       end select
    case( 'diagonalize' )
       select case( trim(option) )
       case('calling')  ; write(iout,'(A)') ' calling LAPACK dysev'
       end select
    end select

    
    flush(iout)


  end subroutine writeme_ham0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! SUBROUTINE ERRORS_GETHAM0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  subroutine errors_getham0(myroutine,option,mystring)

    implicit none
    character(*), intent(in) :: myroutine
    character(*), intent(in) :: option
    character(*), optional,intent(in) :: mystring

    
    select case( trim(myroutine) )
    case( 'diagonalize' )
       select case( trim(option) )
       case( 'illegal' ) 
          write(iout,'(A)') ' ERROR:  illegal value in '//trim(adjustl(mystring))//' ham0 element'
          go to 100
       case('no_converge')
          write(iout,'(A)') ' ERROR:  dysev failed to converge!  INFO = '//trim(adjustl(mystring))
          go to 100
       case('success')
          write(iout,'(A)') ' INFO = '//trim(adjustl(mystring))
          go to 200 
       end select
    end select
    
    
100 write(iout,'(A)') divide
    call dnt(iout)
    write(iout,'(A)') " SOMEDAY YOU'LL SEE THE REASON WHY THERE'S GOOD IN GOODBYE - unknown"
    stop


200 continue

    
  end subroutine errors_getham0
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!    
end module get_ham0
