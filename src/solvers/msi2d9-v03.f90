module msi2d9
  implicit none
  !
  ! Goals:
  !    Given the matrix of coefficients B of a linear system Bx = rhs 
  !    1) Find the lower L and upper U triangular matrix of the factorization
  !       B = LU - P. (P matrix is not important) 
  !       B must have only nine non-null diagonals.
  !    2) Solve the linear system Bx = rhs iteratively using the LU factorization (MSI method)
  !
  !    Latest version: 28 Apr 2011
  !
  ! Subroutines given by this module
  !
  !    1) fb2d9 
  !       (ASSUMES THE COEFFICIENTS OF B MATRIX CORRESPONDING TO
  !        BOUNDARY NEIGHBOORS OF FICTITIOUS VOLUMES ARE ZERO)
  !
  !    2) lu2d9 
  !       (SETS ALPHA=0 AT FICTITIOUS VOLUMES INDEPENDENTLY OF ALPHA ELSEWHERE)
  !       (the parameter alpha must be set inside this subroutine)
contains
  !
  ! fb2d9 ASSUMES THE COEFFICIENTS OF B MATRIX CORRESPONDING TO
  ! BOUNDARY NEIGHBOORS OF FICTITIOUS VOLUMES ARE ZERO
  !
  ! fb2d9 solves the linear system b.var = rhs iteratively
  !
  ! b     - matrix  of coefficients of the original system 
  ! dl    - lower triangular matrix 
  ! du    - upper triangular matrix
  ! nj    - number of volumes in the j direction
  ! ni    - number of volumes in the i direction
  ! nij   - nij=ni*nj
  ! var   - vector of unkowns
  ! gamma - L2( r_n ) < L2( r_0 ) * gamma
  ! L2    - L2 norm
  ! r_n   - residual vector: r_n = b var_n - rhs
  ! nitm  - maximum number of iteractions
  ! rhs   - vector with the source term of the linear system
  !
  subroutine fb2d9(b,dl,du,nj,ni,nij,var,gamma,nitm,rhs)
    implicit none
    integer, intent(in)    ::   ni, nj, nij, nitm
    real(8), intent(in)    ::   gamma    ! 
    real(8), intent(in)    ::   b(nij,9) ! matrix of coefficients
    real(8), intent(in)    ::  dl(nij,5) ! lower triangular matrix
    real(8), intent(in)    ::  du(nij,4) ! upper triangular matrix
    real(8), intent(in)    :: rhs(nij)   ! source vector of the linear system
    real(8), intent(inout) :: var(nij)   ! vector of unknowns
    !
    integer :: nit
    integer :: isw, is, ise, iw, ip, ie, inw, in, ine
    real(8) :: gammarp
    real(8) :: r(nij) ! residual vector
    real(8) :: w(nij) ! vector used to solve linear system Lw = -r
    real(8) :: z(nij) ! vector used to solve linear system Uz =  w
    !
    ! Calculation of the residual vector r = A.var-rhs
    !--------------------------------------------------------------------------------
    ip = 1
    ie  = ip + 1
    in  = ip + ni
    ine = in + 1
    r(ip) = -rhs(ip) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
         +             b(ip,8) * var(in ) + b(ip,9) * var(ine) 
    !--------------------------------------------------------------------------------
    do ip = 2, ni-1
       iw = ip - 1
       ie  = ip + 1
       in  = ip + ni
       inw = in - 1
       ine = in + 1
       r(ip) = -rhs(ip) + b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
            +             b(ip,7) * var(inw) + b(ip,8) * var(in ) + b(ip,9) * var(ine) 
    end do
    !--------------------------------------------------------------------------------
    ip = ni
    iw  = ip - 1
    in  = ip + ni
    inw = in - 1
    r(ip) = -rhs(ip) + b(ip,4) * var(iw ) + b(ip,5) * var(ip ) &
         +             b(ip,7) * var(inw) + b(ip,8) * var(in )                      
    !--------------------------------------------------------------------------------
    do ip = ni+1, nij-ni
       iw  = ip - 1
       ie  = ip + 1
       is  = ip - ni
       in  = ip + ni
       isw = is - 1
       ise = is + 1
       inw = in - 1
       ine = in + 1
       if ( mod(ip,ni) == 1 ) then
          r(ip) = -rhs(ip) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
               +             b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
               +             b(ip,8) * var(in ) + b(ip,9) * var(ine) 
       else if ( mod(ip,ni) == 0 ) then
          r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) &
               +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) &
               +             b(ip,7) * var(inw) + b(ip,8) * var(in ) 
       else
          r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
               +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
               +             b(ip,7) * var(inw) + b(ip,8) * var(in ) + b(ip,9) * var(ine) 
       end if
    end do
    !--------------------------------------------------------------------------------
    ip = nij-ni+1
    ie  = ip + 1
    is  = ip - ni
    ise = is + 1
    r(ip) = -rhs(ip) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
         +             b(ip,5) * var(ip ) + b(ip,6) * var(ie ) 
    !--------------------------------------------------------------------------------
    do ip = nij-ni+2, nij-1
       iw  = ip - 1
       ie  = ip + 1
       is  = ip - ni
       isw = is - 1
       ise = is + 1
       r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
            +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) 
    end do
    !--------------------------------------------------------------------------------
    ip = nij
    iw  = ip - 1
    is  = ip - ni
    isw = is - 1
    r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) &
         +             b(ip,4) * var(iw ) + b(ip,5) * var(ip )
    !--------------------------------------------------------------------------------
    gammarp = gamma * sqrt(dot_product(r,r))
    !--------------------------------------------------------------------------------
    nit = 0
    do
       nit = nit + 1
       !
       ! Solving linear system Lw = -r
       !------------------------------
       ip = 1
       w(ip) = -r(ip) / dl(ip,1)
       !------------------------------
       do ip = 2, ni-1
          w(ip) = -( r(ip) + dl(ip,2)*w(ip-1) ) / dl(ip,1)
       end do
       !--------------------------------------------------
       ip = ni
       w(ip) = -( r(ip) + dl(ip,2)*w(ip-1) + dl(ip,3)*w(ip-ni+1) ) / dl(ip,1)
       !---------------------------------------------------------------------
       ip = ni+1
       w(ip) = -( r(ip) + dl(ip,2)*w(ip-1) + dl(ip,3)*w(ip-ni+1) + dl(ip,4)*w(ip-ni) ) / dl(ip,1)
       !-----------------------------------------------------------------------------------------
       do ip = ni+2, nij
          w(ip) = -( r(ip) + dl(ip,2)*w(ip-1) + dl(ip,3)*w(ip-ni+1) + dl(ip,4)*w(ip-ni) + dl(ip,5)*w(ip-ni-1) ) / dl(ip,1)
       end do
       !-------------------------------------------------------------------------------------------------------------------
       !
       ! Solving linear system Uz = w
       !-----------------------------
       ip = nij
       z(ip) = w(ip)
       !------------
       do ip = nij-1, nij-ni+2, -1
          z(ip) = w(ip) - du(ip,1)*z(ip+1)
       end do
       !----------------------------------
       ip = nij-ni+1
       z(ip) = w(ip) - du(ip,1)*z(ip+1) - du(ip,2)*z(ip+ni-1)
       !-----------------------------------------------------
       ip = nij-ni
       z(ip) = w(ip) - du(ip,1)*z(ip+1) - du(ip,2)*z(ip+ni-1) - du(ip,3)*z(ip+ni)
       !-------------------------------------------------------------------------
       do ip = nij-ni-1, 1, -1
          z(ip) = w(ip) - du(ip,1)*z(ip+1) - du(ip,2)*z(ip+ni-1) - du(ip,3)*z(ip+ni) - du(ip,4)*z(ip+ni+1)
       end do
       !--------------------------------------------------------------------------------------------------
       !
       ! Updating solution
       !------------------
       var = z + var        
       !
       ! Calculation of the residual vector r = A.var-rhs
       !--------------------------------------------------------------------------------
       ip = 1
       ie  = ip + 1
       in  = ip + ni
       ine = in + 1
       r(ip) = -rhs(ip) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
            +             b(ip,8) * var(in ) + b(ip,9) * var(ine) 
       !--------------------------------------------------------------------------------
       do ip = 2, ni-1
          iw = ip - 1
          ie  = ip + 1
          in  = ip + ni
          inw = in - 1
          ine = in + 1
          r(ip) = -rhs(ip) + b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
               +             b(ip,7) * var(inw) + b(ip,8) * var(in ) + b(ip,9) * var(ine) 
       end do
       !--------------------------------------------------------------------------------
       ip = ni
       iw  = ip - 1
       in  = ip + ni
       inw = in - 1
       r(ip) = -rhs(ip) + b(ip,4) * var(iw ) + b(ip,5) * var(ip ) &
            +             b(ip,7) * var(inw) + b(ip,8) * var(in )                      
       !--------------------------------------------------------------------------------
       do ip = ni+1, nij-ni
          iw  = ip - 1
          ie  = ip + 1
          is  = ip - ni
          in  = ip + ni
          isw = is - 1
          ise = is + 1
          inw = in - 1
          ine = in + 1
          if ( mod(ip,ni) == 1 ) then
             r(ip) = -rhs(ip) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
                  +             b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
                  +             b(ip,8) * var(in ) + b(ip,9) * var(ine) 
          else if ( mod(ip,ni) == 0 ) then
             r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) &
                  +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) &
                  +             b(ip,7) * var(inw) + b(ip,8) * var(in ) 
          else
             r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
                  +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) &
                  +             b(ip,7) * var(inw) + b(ip,8) * var(in ) + b(ip,9) * var(ine) 
          end if
       end do
       !--------------------------------------------------------------------------------
       ip = nij-ni+1
       ie  = ip + 1
       is  = ip - ni
       ise = is + 1
       r(ip) = -rhs(ip) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
            +             b(ip,5) * var(ip ) + b(ip,6) * var(ie ) 
       !--------------------------------------------------------------------------------
       do ip = nij-ni+2, nij-1
          iw  = ip - 1
          ie  = ip + 1
          is  = ip - ni
          isw = is - 1
          ise = is + 1
          r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) + b(ip,3) * var(ise) &
               +             b(ip,4) * var(iw ) + b(ip,5) * var(ip ) + b(ip,6) * var(ie ) 
       end do
       !--------------------------------------------------------------------------------
       ip = nij
       iw  = ip - 1
       is  = ip - ni
       isw = is - 1
       r(ip) = -rhs(ip) + b(ip,1) * var(isw) + b(ip,2) * var(is ) &
            +             b(ip,4) * var(iw ) + b(ip,5) * var(ip )
       !--------------------------------------------------------------------------------
       !
       if ( sqrt(dot_product(r,r)) < gammarp ) exit
       if ( nit == nitm ) exit
       !
    end do
  end subroutine fb2d9
  ! 
  ! ld2d9 SETS ALPHA=0 AT FICTITIOUS VOLUMES INDEPENDENTLY OF ALPHA OTHERWISE
  !
  ! lu2d9 calculates the entries of the L and U matrix in the decomposition LU = B + P
  ! where B is the matrix of the coefficients of the linear system to be solved and P is a
  ! matrix whose entries are not necessary.
  !
  ! b     - matrix  of coefficients of the original system (has 9 diagonals)
  ! dl    - lower triangular matrix 
  ! du    - upper triangular matrix
  ! nj    - number of volumes in the j direction
  ! ni    - number of volumes in the i direction
  ! nij   - nij=ni*nj
  !
  subroutine lu2d9(b,dl,du,nij,ni,nj)
    implicit none
    integer, intent(in)  :: ni, nj, nij
    real(8), intent(in)  ::  b(nij,9)  ! matrix of coefficients
    real(8), intent(out) :: dl(nij,5)  ! lower diagonal matrix
    real(8), intent(out) :: du(nij,4)  ! upper diagonal matrix

    integer :: isw, is, ise, iw, ip
    real(8) :: alpha, phi1, phi2, phi3, phi4, alpha0

    alpha0= 0.5d0
    alpha = alpha0

    du = 0.d0
    dl = 0.d0

    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = 1
    !
    dl(ip,1) = b(ip,5)
    du(ip,1) = b(ip,6) / dl(ip,1)
    du(ip,2) = b(ip,7) / dl(ip,1)
    du(ip,3) = b(ip,8) / dl(ip,1)
    du(ip,4) = b(ip,9) / dl(ip,1)
    !------------------------------------------------------------------------------------------------------------------------
    !
    do ip = 2, ni-1
       iw  = ip - 1
       dl(ip,2) = b(ip,4)
       dl(ip,1) = b(ip,5) - dl(ip,2)*du(iw,1)
       du(ip,1) = b(ip,6) / dl(ip,1)
       du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
       du(ip,3) = ( b(ip,8) - dl(ip,2)*du(iw,4) ) / dl(ip,1)
       du(ip,4) = b(ip,9) / dl(ip,1)
    end do
    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = ni
    !
    iw  = ip - 1
    is  = ip - ni
    ise = is + 1
    dl(ip,3) = b(ip,3)
    dl(ip,2) = b(ip,4)
    dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)
    du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) ) / dl(ip,1)
    du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
    du(ip,3) = ( b(ip,8) - dl(ip,2)*du(iw,4) ) / dl(ip,1)
    du(ip,4) = b(ip,9) / dl(ip,1)
    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = ni+1
    !
    iw  = ip - 1
    is  = ip - ni
    ise = is + 1
    dl(ip,4) = b(ip,2)
    dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
    dl(ip,2) = b(ip,4) - dl(ip,4)*du(is,2) 
    dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)
    du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) - dl(ip,4)*du(is,4) ) / dl(ip,1)
    du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
    du(ip,3) = ( b(ip,8) - dl(ip,2)*du(iw,4) ) / dl(ip,1)
    du(ip,4) = b(ip,9) / dl(ip,1)
    !------------------------------------------------------------------------------------------------------------------------
    !
    alpha = alpha0
    do ip = ni+2, nij-ni-1
       !
       if ( mod(ip,ni) == 0 .or. mod(ip,ni) == 1 ) then
          alpha = 0.d0
       else
          alpha = alpha0
       end if
       iw  = ip - 1
       is  = ip - ni
       isw = is - 1
       ise = is + 1
       dl(ip,5) = b(ip,1)
       dl(ip,4) = ( b(ip,2) -dl(ip,5)*du(isw,1) -alpha*b(ip,3)*du(ise,1) ) / ( 1.d0 - alpha*du(is,1)*du(ise,1) )
       dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
       dl(ip,2) = (b(ip,4) - dl(ip,4)*du(is,2) -dl(ip,5)*du(isw,3) -2.d0*alpha*dl(ip,5)*du(isw,2) ) / ( 1 + 2.d0*alpha*du(iw,2) )
       phi1 = dl(ip,3) * du(ise,1)
       phi2 = dl(ip,5) * du(isw,2)
       phi3 = dl(ip,3) * du(ise,4)
       phi4 = dl(ip,2) * du(iw ,2)
       dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)-dl(ip,5)*du(isw,4) &
            +     alpha * ( 2.d0*phi1 + phi2 + phi3 + 2.d0*phi4 )
       du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) - dl(ip,4)*du(is,4) - 2.d0*alpha*(phi1+phi3) ) / dl(ip,1)
       du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
       du(ip,3) = ( b(ip,8) - dl(ip,2)*du(iw,4) - alpha*phi4 ) / dl(ip,1)
       du(ip,4) = b(ip,9) / dl(ip,1)
    end do
    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = nij-ni
    !
    iw  = ip - 1
    is  = ip - ni
    isw = is - 1
    ise = is + 1
    dl(ip,5) = b(ip,1)
    dl(ip,4) = b(ip,2) - dl(ip,5)*du(isw,1) 
    dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
    dl(ip,2) = b(ip,4) - dl(ip,4)*du(is,2) -dl(ip,5)*du(isw,3)
    dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)-dl(ip,5)*du(isw,4)
    du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) - dl(ip,4)*du(is,4) ) / dl(ip,1)
    du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
    du(ip,3) = ( b(ip,8) - dl(ip,2)*du(iw,4) ) / dl(ip,1)
    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = nij-ni+1
    !
    iw  = ip - 1
    is  = ip - ni
    isw = is - 1
    ise = is + 1
    dl(ip,5) = b(ip,1)
    dl(ip,4) = b(ip,2) -dl(ip,5)*du(isw,1)
    dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
    dl(ip,2) = b(ip,4) - dl(ip,4)*du(is,2) -dl(ip,5)*du(isw,3)
    dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)-dl(ip,5)*du(isw,4) 
    du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) - dl(ip,4)*du(is,4) ) / dl(ip,1)
    du(ip,2) = ( b(ip,7) - dl(ip,2)*du(iw,3) ) / dl(ip,1)
    !------------------------------------------------------------------------------------------------------------------------
    !
    do ip = nij-ni+2, nij-1
       !
       iw  = ip - 1
       is  = ip - ni
       isw = is - 1
       ise = is + 1
       dl(ip,5) = b(ip,1)
       dl(ip,4) = b(ip,2) -dl(ip,5)*du(isw,1)
       dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
       dl(ip,2) = b(ip,4) - dl(ip,4)*du(is,2) -dl(ip,5)*du(isw,3)
       dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)-dl(ip,5)*du(isw,4) 
       du(ip,1) = ( b(ip,6) - dl(ip,3)*du(ise,3) - dl(ip,4)*du(is,4) ) / dl(ip,1)
    end do
    !------------------------------------------------------------------------------------------------------------------------
    !
    ip = nij
    !
    iw  = ip - 1
    is  = ip - ni
    isw = is - 1
    ise = is + 1
    dl(ip,5) = b(ip,1)
    dl(ip,4) = b(ip,2) -dl(ip,5)*du(isw,1)
    dl(ip,3) = b(ip,3) - dl(ip,4) * du(is,1)
    dl(ip,2) = b(ip,4) - dl(ip,4)*du(is,2) -dl(ip,5)*du(isw,3)
    dl(ip,1) = b(ip,5)-dl(ip,2)*du(iw,1)-dl(ip,3)*du(ise,2)-dl(ip,4)*du(is,3)-dl(ip,5)*du(isw,4)
    !------------------------------------------------------------------------------------------------------------------------
    !
  end subroutine lu2d9

end module msi2d9
