#define ALP_PC

module local

  use elems
  implicit none
  ! Legacy; check ipoints(:) and weights(:)
  ! These will be removed
  real(8),allocatable :: ipoint(:,:),weight(:)
  
  ! Integation points for types
  ! ipoints(getElTypeNo(<TYPE_NAME>)%array) --> Retains integration points for the type
  ! weights(getElTypeNo(<TYPE_NAME>)%array) --> Retains weights for the type
  type(array2d), target :: ipoints(elTypeCount)
  type(array1d), target :: weights(elTypeCount)

contains

  !!$ ! --- Form Elastic stiffness matrix
  ! Edit: Removed estress  
  subroutine FormElKE(eltype, ecoords, E,nu,dt,k)
    implicit none
    character(*) :: eltype
    integer :: nodecount, typeIndex, nip
    real(8) :: ecoords(:,:) !---> coords of nodes
    real(8) :: E,nu,G, Kb, dt !--> material properties
    ! ---> E --> Youngs modulus
    ! ---> nu --> poisson ration
    ! ---> Rest dont know yet
    real(8) :: D(cpdim,cpdim)
    ! ---> D ---> Dmat ---> Matrix of elastic coefficients
    real(8) :: detj
    !          dN   ----> derivative of shape funtions
    !          detj ----> determinant of jacobian
    !          B    ----> B matrix
    ! Dimensions:
    !          dN(pdim, nodecount)
    !          B(cpdim, nodecount * pdim)
    !          N(1,nodecount)
    ! nodecount = getNodeCount(eltype)
    real(8), allocatable ::  dN(:, :), B(:, :), N(:,:)
    real(8), target :: k(:,:)
    integer :: i
    real(8) :: Q(pdim,pdim),pN,m(cpdim,1)
    real(8) :: s
    ! Removed variables:
    ! real(8), pointer :: kl(:,:)
    ! real(8) :: T(pdim, eldof), estress(:,:)
       
    nodeCount = getNodeCount(eltype)
    allocate(dN(pdim, nodecount), B(cpdim, nodecount*pdim), N(1, nodecount))
       
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    pN=f1/dble(nodecount)
    G=E/(f2*(f1+nu))
    Kb=E/(f3*(f1-f2*nu))
    k = f0; ! Initialize k to zero
    call DMat(D,E,nu) ! --- Get Dmat   
    ! ---- Loop over integration points and form stiffness matrix
    do i = 1, nip        
       call FormdNdetJ(eltype, ipoints(typeIndex)%array(i,:),ecoords,dN,detj)   
       call BMat(nodeCount, dN,B)   
       k = k+matmul(transpose(B),matmul(D,B))*weights(typeIndex)%array(i)*detj
    end do
    deallocate(dN, B, N)
  end subroutine FormElKE

  ! Form element index
  subroutine FormElIndx(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:),j,j1
    do j=1,size(enodes)
       do j1=1,pdim
          indx(pdim*j-j1+1)=pdim*enodes(j)-j1+1
       end do
    end do
  end subroutine FormElIndx

  ! Form element indexp
  subroutine FormElIndxp(enodes,indxp)
    implicit none
    integer :: enodes(nodesperel),indxp(eldof+eldofp),j,j1
    do j=1,nodesperel
       do j1=1,pdim
          indxp(pdim*j-j1+1)=(pdim+1)*enodes(j)-j1
       end do
       indxp(nodesperel*pdim+j)=(pdim+1)*enodes(j)
    end do
  end subroutine FormElIndxp

  ! Calculate element stress
  subroutine CalcElStress(eltype, ecoords,edisp,E,nu,estress)
    implicit none
    character(*) :: eltype
    ! Dimensions:
    !         ecoords(nodecount, pdim)
    !         edisp(nodecount * pdim)
    !         dN(pdim, nodecount)
    !         B(cpdim, nodecount * pdim)
    !         estress(getNip(eltype), cpdim)
    !         estrain(getNip(eltype), cpdim)
    real(8), allocatable :: dN(:, :), B(:, :), estrain(:, :)
    real(8) :: ecoords(:,:),edisp(:),E,nu,          &
       estress(:, :),D(cpdim,cpdim),detj
    integer :: i, nodecount, typeIndex, nip
    
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    typeIndex = getElTypeNo(eltype)
    allocate(dN(pdim, nodecount), B(cpdim, nodecount * pdim))
    allocate(estrain(nip, cpdim))
    call DMat(D,E,nu)     

    do i=1,nip
      call FormdNdetJ(eltype, ipoints(typeIndex)%array(i,:),ecoords,dN,detj)
      call BMat(nodecount, dN,B) 
       estrain(i,:)=matmul(B,edisp)
       estress(i,:)=matmul(D,estrain(i,:))
    end do
    
    deallocate(dN, B, estrain)
  end subroutine CalcElStress

  ! Reform element RHS
  ! NOT EDITED
  subroutine ReformElRHS(ecoords,estress,E,nu,visc,expn,dt,f)
    implicit none
    real(8) :: ecoords(nodesperel,pdim),estress(nip,cpdim),E,nu,visc,expn,dt,f(eldof), &
       alpha,D(cpdim,cpdim),S(cpdim,cpdim),dN(pdim,nodesperel),detj,B(cpdim,eldof),        &
       beta(cpdim),betad(cpdim,cpdim)
    integer :: i
    call DMat(D,E,nu); call MatInv(D,S)
    alpha=f1;f=f0
    do i=1,nip
       ! TODO DUMMY ARG1, EXCISE
       call FormdNdetJ("tri", ipoint(i,:),ecoords,dN,detj)
       ! TODO DUMMY ARG1, EXCISE
       call BMat(0, dN,B)
       call Matbetad(betad,estress(i,:),visc,expn)
       call MatInv(S+alpha*dt*betad,D)
       call Matbeta(beta,estress(i,:),visc,expn)
       f=f+matmul(transpose(B),matmul(D,dt*beta))*weight(i)*detj
    end do
  end subroutine ReformElRHS

  ! Computes the strain displacement matrix 'B'
  subroutine BMat(nodecount, dN,B)
    implicit none
    integer :: nodecount
    real(8) :: dN(pdim,nodecount),B(cpdim,nodecount*pdim)
    integer :: j
    B=f0
    select case(pdim)
    case(2)
       do j=1,nodecount
          B(1,2*j-1)=dN(1,j); B(1,2*j)=f0
          B(2,2*j-1)=f0     ; B(2,2*j)=dN(2,j)
          B(3,2*j-1)=dN(2,j); B(3,2*j)=dN(1,j)
       end do
    case(3)
       do j=1,nodecount
          B(1,3*j-2)=dN(1,j); B(1,3*j-1)=f0     ; B(1,3*j)=f0
          B(2,3*j-2)=f0     ; B(2,3*j-1)=dN(2,j); B(2,3*j)=f0
          B(3,3*j-2)=f0     ; B(3,3*j-1)=f0     ; B(3,3*j)=dN(3,j)
          B(4,3*j-2)=dN(2,j); B(4,3*j-1)=dN(1,j); B(4,3*j)=f0
          B(5,3*j-2)=f0     ; B(5,3*j-1)=dN(3,j); B(5,3*j)=dN(2,j)
          B(6,3*j-2)=dN(3,j); B(6,3*j-1)=f0     ; B(6,3*j)=dN(1,j)
       end do
    end select
  end subroutine BMat
  
 
  ! Computes dN/dx(yz) and the determinant of Jacobian 'detj'
  subroutine FormdNdetJ(eltype, ipcoord,ecoords,dN,detj)
    implicit none
    character(*) :: eltype
    ! Dimensions:
    !             ecoords(nodecount, pdim)
    !             dN(pdim, nodecount)
    ! Note: nodecount = getNodeCount(eltype)
    real(8) :: ipcoord(pdim),detj,ecoords(:, :),dN(:, :),               &
       jacob(pdim,pdim),invj(pdim,pdim)
    call ShapeFuncd(eltype, dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
    call MatInv(jacob,invj)
    dN=matmul(invj,dN) ! Form dN/dx, dN/dy, (dN/dz)
  end subroutine FormdNdetJ

  ! Computes only the determinant of Jacobian 'detj' for use in Hs
  ! NOT EDITED
  subroutine FormdetJ(ipcoord,ecoords,detj)
    implicit none
    real(8) :: ipcoord(pdim),detj,ecoords(nodesperel,pdim),dN(pdim,nodesperel),               &
       jacob(pdim,pdim)
    ! TODO DUMMY ARG1, EXCISE
    call ShapeFuncd("tri", dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
  end subroutine FormdetJ

  ! Computes D, the matrix of elastic properties
  subroutine DMat(D,E,nu)
    implicit none
    real(8) :: D(:,:),E,nu
    if (size(D,1)==3) call DMat2d(D,E,nu)
    if (size(D,1)==6) call DMat3d(D,E,nu)
  end subroutine DMat

  ! Computes D, the matrix of elastic properties
  subroutine DMat2d(D,E,nu)
    implicit none
    real(8) :: D(3,3),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,f0,nu,(f1-nu),f0,f0,f0,(f1-f2*nu)/f2/),           &
       shape=(/3,3/))
  end subroutine DMat2d

  ! Computes D, the matrix of elastic properties
  subroutine DMat3d(D,E,nu)
    implicit none
    real(8) :: D(6,6),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,nu,f0,f0,f0,nu,(f1-nu),nu,f0,f0,f0,nu,nu,(f1-nu), &
       f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,  &
       f0,f0,f0,f0,(f1-f2*nu)/f2/),shape=(/6,6/))
  end subroutine DMat3d

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(:),estress(:),visc,expn
    if (size(estress,1)==3) call Matbeta2d(beta,estress,visc,expn)
    if (size(estress,1)==6) call Matbeta3d(beta,estress,visc,expn)
  end subroutine Matbeta

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta2d(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(3),estress(3),kappa,visc,expn,cMat(3,3),s1,s2,s3
    s1=estress(1); s2=estress(2); s3=estress(3)
    kappa=sqrt(((s1-s2)/f2)**2+s3**2)
    cMat=reshape((/f1,-f1,f0,-f1,f1,f0,f0,f0,f4/),shape=(/3,3/))
    beta=((kappa**(expn-f1))/(f4*visc))*matmul(cMat,estress)
  end subroutine Matbeta2d

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta3d(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(6),estress(6),kappa,visc,expn,cMat(6,6),s1,s2,s3,s4,s5,s6
    s1=estress(1); s2=estress(2); s3=estress(3)
    s4=estress(4); s5=estress(5); s6=estress(6)
    kappa=sqrt(((s1-s2)**2+(s2-s3)**2+(s1-s3)**2)/f6+s4**2+s5**2+s6**2)
    cMat=reshape((/                                                            &
        f4/f3,-f2/f3,-f2/f3, f0 , f0 , f0 ,                                    &
       -f2/f3, f4/f3,-f2/f3, f0 , f0 , f0 ,                                    &
       -f2/f3,-f2/f3, f4/f3, f0 , f0 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f4 , f0 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f0 , f4 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f0 , f0 , f4 /),shape=(/6,6/))
    beta=((kappa**(expn-f1))/(f4*visc))*matmul(cMat,estress)
  end subroutine Matbeta3d

  ! Computes \beta', Jacobian Matrix (differentiate \beta w.r.t. components of
  ! stress)
  subroutine Matbetad(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(:,:),estress(:),visc,expn
    if (size(estress,1)==3) call Matbetad2d(betad,estress,visc,expn)
    if (size(estress,1)==6) call Matbetad3d(betad,estress,visc,expn)
  end subroutine Matbetad

  ! Computes \beta', Jacobian Matrix (differentiate \beta w.r.t. components of
  ! stress)
  subroutine Matbetad2d(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(3,3),estress(3),kappa,c1,c2,c3,visc,expn,s1,s2,s3
    s1=estress(1); s2=estress(2); s3=estress(3)
    kappa=sqrt(((s1-s2)/f2)**2+s3**2)
    if (kappa==f0) betad=f0
    if (kappa==f0) return
    c1=f1+(expn-f1)*((s1-s2)/(f2*kappa))**2
    c2=f1+(expn-f1)*(s3/kappa)**2
    c3=(expn-f1)*(s1*s3-s2*s3)/(kappa)**2
    betad=((kappa**(expn-f1))/(f4*visc))*reshape((/c1,-c1,c3,-c1,c1,-c3,c3,    &
       -c3,f4*c2/),shape=(/3,3/))
  end subroutine Matbetad2d

  ! Computes \beta', Jacobian Matrix (differentiate \beta w.r.t. components of
  ! stress)
  subroutine Matbetad3d(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(6,6),estress(6),kappa,visc,expn,c,Sx,Sy,Sz,T1,T2,T3,s1,   &
       s2,s3,s4,s5,s6
    s1=estress(1); s2=estress(2); s3=estress(3)
    s4=estress(4); s5=estress(5); s6=estress(6)
    kappa=sqrt(((s1-s2)**2+(s2-s3)**2+(s1-s3)**2)/f6+s4**2+s5**2+s6**2)
    if (kappa==f0) betad=f0
    if (kappa==f0) return
    c=sqrt(expn-f1)
    Sx=c*(f2*s1-s2-s3)/(f3*kappa)
    Sy=c*(f2*s2-s3-s1)/(f3*kappa)
    Sz=c*(f2*s3-s1-s2)/(f3*kappa)
    T1=c*f2*s4/kappa; T2=c*f2*s5/kappa; T3=c*f2*s6/kappa
    betad=((kappa**(expn-f1))/(f4*visc))*reshape((/                            &
          f4/f3+Sx**2,-f2/f3+Sx*Sy,-f2/f3+Sx*Sz,  Sx*T1 ,  Sx*T2 ,  Sx*T3 ,    &
         -f2/f3+Sx*Sy, f4/f3+Sy**2,-f2/f3+Sy*Sz,  Sy*T1 ,  Sy*T2 ,  Sy*T3 ,    &
         -f2/f3+Sx*Sz,-f2/f3+Sy*Sz, f4/f3+Sz**2,  Sz*T1 ,  Sz*T2 ,  Sz*T3 ,    &
             Sx*T1   ,    Sy*T1   ,    Sz*T1   ,f4+T1**2,  T1*T2 ,  T1*T3 ,    &
             Sx*T2   ,    Sy*T2   ,    Sz*T2   ,  T2*T1 ,f4+T2**2,  T2*T3 ,    &
             Sx*T3   ,    Sy*T3   ,    Sz*T3   ,  T3*T1 ,  T3*T2 ,f4+T3**2/),  &
          shape=(/6,6/))
  end subroutine Matbetad3d

end module local
