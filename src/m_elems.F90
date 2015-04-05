module elems

  use utils
  implicit none
  character(3) :: eltype
  integer :: pdim,nodesperel,nps,eldof,eldofp,cpdim, nip
  
  
  ! Constant element info
  integer, parameter :: elTypeCount = 5
  character(4) :: elTypes(elTypeCount) = (/"tri", "qua", "tet", "hex", "coh"/)
  integer :: elNodeCount(elTypeCount) = (/3, 4, 4, 6, 4/)
  
  ! shapeFuncMem: ShapeFunc Memoization 
  ! shapeFuncMem(typeIndex)%array(nip) -> ShapeFunc(elType, nip)
  type(array1dWrapper), target :: shapeFuncMem(elTypeCount)
  ! nodalStressInvMem: Memoization of the N2Inv array during the stress calculation
  ! nodalStressInvMem(typeIndex)%array -> N2Inv of given element type
  type(array2d) :: nodalStressInvMem(elTypeCount)

contains

  ! Return the pertaining dimension of this element
  function getDim(eltype) result(pdim)
    implicit none
    character(*) :: eltype
    integer :: pdim
    select case(eltype)
    case("tri"); pdim = 2
    case("qua"); pdim = 2
    case("tet"); pdim = 3
    case("hex"); pdim = 3
    case("coh"); pdim = 2
    end select
  end function getDim
  ! Return the numerical value of the eltype
  ! NOTE: Exhaustive search can be altered into a logn search if
  ! eltypes proliferate
  function getElTypeNo(eltype) result(elTypeNo)
    implicit none
    character(*) :: eltype
    integer :: elTypeNo, i
    do i = 1, elTypeCount
      if(elTypes(i) == eltype) then
        elTypeNo = i
        exit
      end if
    end do
  end function getElTypeNo
  
  function getElName(elTypeNo) result(elName)
    implicit none
    integer :: elTypeNo
    character(4) :: elname
    
    elname = eltypes(elTypeNo)
  end function getElName

  ! Return the element type nodecount
  function getNodeCount(eltype) result(nodeCount)
    character(*) :: eltype
    integer :: nodeCount, typeIndex
    
    typeIndex = getElTypeNo(eltype)
    nodeCount = elNodeCount(typeIndex)   
  end function getNodeCount
  
  ! Return the element's number of integration points
  function getNip(eltype) result(nip)
    implicit none
    character(*) :: eltype
    integer :: nip
    
    select case(eltype)
    case("tri"); nip = 1
    case("qua"); nip = 4
    case("tet"); nip = 1
    case("hex"); nip = 8
    case("coh"); nip = 2
    end select
  end function getNip
  
  ! Return the element's vtkid
  function getVtkid(eltype) result(vtkid)
    implicit none
    character(*) :: eltype
    integer :: vtkid
    
    select case(eltype)
    case("tri"); vtkid = 5
    case("qua"); vtkid = 9
    case("tet"); vtkid = 10
    case("hex"); vtkid = 12
    case("coh"); vtkid = 9
    end select
  end function getVtkid

  ! Sets element specific constants
  ! Initializes caches
  subroutine Initialize(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    if (pdim==2) cpdim=3
    if (pdim==3) cpdim=6
    call SamPts(ipoints, weights)  
    call ShapeFuncPrecomp(ipoints)
    call NodalStressInv
  end subroutine Initialize
  
  ! Get the shape function
  subroutine ShapeFunc(eltype, nip, N)
    character(*), intent(in) :: eltype
    integer, intent(in) :: nip
    real(8), pointer, intent(out) :: N(:)
    N=>ShapeFuncMem(getElTypeNo(eltype))%array(nip)%array
  end subroutine ShapeFunc
  
  ! Calculate quadrature ipoints and weights
  subroutine SamPts(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    
    if(pdim == 2) then
      call SamPtsTri(ipoints, weights)
      call SamPtsQua(ipoints, weights)
      call SamPtsCoh(ipoints, weights)
    end if
    
    if(pdim == 3) then
      call SamPtsTet(ipoints, weights)
      call SamPtsHex(ipoints, weights)
    end if
  end subroutine SamPts

  ! Precomputes 'N', the shape functions
  subroutine ShapeFuncPrecomp(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    
    if(pdim == 2) then
      call ShapeFuncPrecompTri(ipoints)
      call ShapeFuncPrecompQua(ipoints)
      call ShapeFuncPrecompCoh(ipoints)
    else if (pdim == 3) then
      call ShapeFuncPrecompTet(ipoints)
      call ShapeFuncPrecompHex(ipoints)
    end if
  end subroutine ShapeFuncPrecomp
  
  ! Computes 'dN', derivative of shape functions
  subroutine ShapeFuncd(eltype, dN,coord)
    implicit none
    character(*) :: eltype
    real(8) :: dN(:,:),coord(:)
    if (eltype=="tri") call ShapeFuncdTri(dN,coord)
    if (eltype=="qua") call ShapeFuncdQua(dN,coord)
    if (eltype=="tet") call ShapeFuncdTet(dN,coord)
    if (eltype=="hex") call ShapeFuncdHex(dN,coord)
  end subroutine ShapeFuncd

  ! Computes the element volume
  subroutine ElVol(ecoords,vol)
    implicit none
    real(8) :: ecoords(:,:),vol
    if (eltype=="tri") call ElVolTri(ecoords,vol)
    if (eltype=="qua") call ElVolQua(ecoords,vol)
    if (eltype=="tet") call ElVolTet(ecoords,vol)
    if (eltype=="hex") call ElVolHex(ecoords,vol)
  end subroutine ElVol

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(:),side,snodes(:)
    real(8) :: area, ecoords(:,:)
    if (eltype=="tri") call EdgeAreaNodesTri(enodes,ecoords,side,area,snodes)
    if (eltype=="qua") call EdgeAreaNodesQua(enodes,ecoords,side,area,snodes)
    if (eltype=="tet") call EdgeAreaNodesTet(enodes,ecoords,side,area,snodes)
    if (eltype=="hex") call EdgeAreaNodesHex(enodes,ecoords,side,area,snodes)
  end subroutine EdgeAreaNodes

! Linear Tri ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsTri(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    integer :: nip, typeIndex    
    nip = getNip("tri")
    typeIndex = getElTypeNo("tri")
    
    if(allocated(ipoints(typeIndex)%array)) deallocate(ipoints(typeIndex)%array)
    if(allocated(weights(typeIndex)%array)) deallocate(weights(typeIndex)%array)
    allocate(ipoints(typeIndex)%array(nip, pdim))
    allocate(weights(typeIndex)%array(nip))
    
    if (nip==1) then
       ipoints(typeIndex)%array(1,:)=(/f1/f3,f1/f3/)
       weights(typeIndex)%array=0.5d0 ! Weight = Wi = wi*wj
    end if
    if (nip==3) then
       ipoints(typeIndex)%array(1,:)=(/f1/f6,f1/f6/)
       ipoints(typeIndex)%array(2,:)=(/f4/f6,f1/f6/)
       ipoints(typeIndex)%array(3,:)=(/f1/f6,f4/f6/)
       weights(typeIndex)%array=f1/f6 ! Weight = Wi = wi*wj
    end if
  end subroutine SamPtsTri

  ! Computes the element volume
  subroutine ElVolTri(ecoords,vol)
    implicit none
    real(8) :: ecoords(3,2),vol,s,l1,l2,l3
    l1=sqrt((ecoords(1,1)-ecoords(2,1))**2+(ecoords(1,2)-ecoords(2,2))**2)
    l2=sqrt((ecoords(2,1)-ecoords(3,1))**2+(ecoords(2,2)-ecoords(3,2))**2)
    l3=sqrt((ecoords(3,1)-ecoords(1,1))**2+(ecoords(3,2)-ecoords(1,2))**2)
    s=0.5d0*(l1+l2+l3)
    vol=sqrt(s*(s-l1)*(s-l2)*(s-l3))
  end subroutine ElVolTri

  ! Computes 'N', the shape functions
  subroutine ShapeFuncPrecompTri(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    character(*), parameter :: eltype = "tri"
    integer :: nip, typeIndex, nodecount, i
    real(8) :: eta,nu   
    
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    
    if(allocated(shapeFuncMem(typeIndex)%array)) deallocate(shapeFuncMem(typeIndex)%array)
    allocate(shapeFuncMem(typeIndex)%array(nip))
    do i = 1, nip
      allocate(shapeFuncMem(typeIndex)%array(i)%array(nodecount))
      eta = ipoints(typeIndex)%array(i, 1)
      nu = ipoints(typeIndex)%array(i, 2)
      shapeFuncMem(typeIndex)%array(i)%array(1) = f1 - eta - nu
      shapeFuncMem(typeIndex)%array(i)%array(2) = eta
      shapeFuncMem(typeIndex)%array(i)%array(3) = nu
    end do
  end subroutine ShapeFuncPrecompTri

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta' and 'nu'
  subroutine ShapeFuncdTri(dN,coord)
    implicit none
    real(8) :: dN(2,3),coord(2),e,n
    e=coord(1); n=coord(2)
    dN(1,:)=(/-f1,f1,f0/) ! dN/de
    dN(2,:)=(/-f1,f0,f1/) ! dN/dn
  end subroutine ShapeFuncdTri

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesTri(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(3),side,snodes(2),j(2)
    real(8) :: ecoords(3,2),area
    if (side==1) then; j(1)=1; j(2)=2
    else if (side==2) then; j(1)=2; j(2)=3
    else if (side==3) then; j(1)=3; j(2)=1
    end if
    snodes=enodes(j)
    area=sqrt((ecoords(j(1),1)-ecoords(j(2),1))**2+                            &
              (ecoords(j(1),2)-ecoords(j(2),2))**2)
  end subroutine EdgeAreaNodesTri

! Linear Quad ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsQua(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    integer :: nip, typeIndex    
    nip = getNip("qua")
    typeIndex = getElTypeNo("qua")
    
    if(allocated(ipoints(typeIndex)%array)) deallocate(ipoints(typeIndex)%array)
    if(allocated(weights(typeIndex)%array)) deallocate(weights(typeIndex)%array)
    allocate(ipoints(typeIndex)%array(nip, pdim))
    allocate(weights(typeIndex)%array(nip))
    
    ipoints(typeIndex)%array(1,:)=(/-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(2,:)=(/-sqrt(f1/f3), sqrt(f1/f3)/)
    ipoints(typeIndex)%array(3,:)=(/ sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(4,:)=(/ sqrt(f1/f3), sqrt(f1/f3)/)
    weights(typeIndex)%array=f1 ! Weight = Wi = wi*wj
  end subroutine SamPtsQua

  ! Computes the element volume
  subroutine ElVolQua(ecoords,vol)
    implicit none
    real(8) :: ecoords(4,2),vol,l1,l2,l3,l4,d1,d2
    l1=sqrt((ecoords(1,1)-ecoords(2,1))**2+(ecoords(1,2)-ecoords(2,2))**2)
    l2=sqrt((ecoords(2,1)-ecoords(3,1))**2+(ecoords(2,2)-ecoords(3,2))**2)
    l3=sqrt((ecoords(3,1)-ecoords(4,1))**2+(ecoords(3,2)-ecoords(4,2))**2)
    l4=sqrt((ecoords(4,1)-ecoords(1,1))**2+(ecoords(4,2)-ecoords(1,2))**2)
    d1=sqrt((ecoords(1,1)-ecoords(3,1))**2+(ecoords(1,2)-ecoords(3,2))**2)
    d2=sqrt((ecoords(2,1)-ecoords(4,1))**2+(ecoords(2,2)-ecoords(4,2))**2)
    vol=0.25d0*sqrt(f4*(d1**2)*(d2**2)-(l1**2+l3**2-l2**2-l4**2)**2)
  end subroutine ElVolQua

  ! Computes 'N', the shape functions
  subroutine ShapeFuncPrecompQua(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    character(*), parameter :: eltype = "qua"
    integer :: nip, typeIndex, nodecount, i
    real(8) :: eta,nu   
    
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    
    if(allocated(shapeFuncMem(typeIndex)%array)) deallocate(shapeFuncMem(typeIndex)%array)
    allocate(shapeFuncMem(typeIndex)%array(nip))
    do i = 1, nip
      allocate(shapeFuncMem(typeIndex)%array(i)%array(nodecount))
      eta = ipoints(typeIndex)%array(i, 1)
      nu = ipoints(typeIndex)%array(i, 2)
      shapeFuncMem(typeIndex)%array(i)%array(1) = 0.25d0*(f1-eta)*(f1-nu)
      shapeFuncMem(typeIndex)%array(i)%array(2) = 0.25d0*(f1+eta)*(f1-nu)
      shapeFuncMem(typeIndex)%array(i)%array(3) = 0.25d0*(f1+eta)*(f1+nu)
      shapeFuncMem(typeIndex)%array(i)%array(4) = 0.25d0*(f1-eta)*(f1+nu)
    end do
  end subroutine ShapeFuncPrecompQua

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta' and 'nu'
  subroutine ShapeFuncdQua(dN,coord)
    implicit none
    real(8) :: dN(2,4),coord(2),e,n
    e=coord(1); n=coord(2)
    dN(1,:)=0.25d0*(/-(f1-n), (f1-n),(f1+n),-(f1+n)/) ! dN/de
    dN(2,:)=0.25d0*(/-(f1-e),-(f1+e),(f1+e), (f1-e)/) ! dN/dn
  end subroutine ShapeFuncdQua

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesQua(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(4),side,snodes(2),j(2)
    real(8) :: ecoords(4,2),area
    if (side==1) then; j(1)=1; j(2)=2
    else if (side==2) then; j(1)=2; j(2)=3;
    else if (side==3) then; j(1)=3; j(2)=4;
    else if (side==4) then; j(1)=4; j(2)=1;
    end if
    snodes=enodes(j)
    area=sqrt((ecoords(j(1),1)-ecoords(j(2),1))**2+                            &
              (ecoords(j(1),2)-ecoords(j(2),2))**2)
  end subroutine EdgeAreaNodesQua

! Linear Tet ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsTet(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    real(8),parameter :: f24=24.0d0,a=0.58541020d0,b=0.13819660d0   
    integer :: nip, typeIndex    

    nip = getNip("tet")
    typeIndex = getElTypeNo("tet")
    
    if(allocated(ipoints(typeIndex)%array)) deallocate(ipoints(typeIndex)%array)
    if(allocated(weights(typeIndex)%array)) deallocate(weights(typeIndex)%array)
    allocate(ipoints(typeIndex)%array(nip, pdim))
    allocate(weights(typeIndex)%array(nip))
 
    if (nip==1) then
       ipoints(typeIndex)%array(1,:)=(/f1/f4,f1/f4,f1/f4/)
       weights(typeIndex)%array=f1/f6 ! Weight = Wi = wi*wj*wk
    end if
    if (nip==4) then
       ipoints(typeIndex)%array(1,:)=(/a,b,b/) ! a=(5.0+3.0*sqrt(5.0))/20.0
       ipoints(typeIndex)%array(2,:)=(/b,a,b/) ! b=(5.0-sqrt(5.0))/20.0
       ipoints(typeIndex)%array(3,:)=(/b,b,a/)
       ipoints(typeIndex)%array(4,:)=(/b,b,b/)
       weights(typeIndex)%array=f1/f24 ! Weight = Wi = wi*wj*wk
    end if
  end subroutine SamPtsTet

  ! Computes the element volume
  subroutine ElVolTet(ecoords,vol)
    implicit none
    real(8) :: ecoords(4,3),vol,M(4,4)
    M=reshape((/ f1,ecoords(1,:),f1,ecoords(2,:),f1,ecoords(3,:),f1,           &
       ecoords(4,:) /),shape=(/4,4/))
    vol=M(1,4)*M(2,3)*M(3,2)*M(4,1)-M(1,3)*M(2,4)*M(3,2)*M(4,1)-               &
        M(1,4)*M(2,2)*M(3,3)*M(4,1)+M(1,2)*M(2,4)*M(3,3)*M(4,1)+               &
        M(1,3)*M(2,2)*M(3,4)*M(4,1)-M(1,2)*M(2,3)*M(3,4)*M(4,1)-               &
        M(1,4)*M(2,3)*M(3,1)*M(4,2)+M(1,3)*M(2,4)*M(3,1)*M(4,2)+               &
        M(1,4)*M(2,1)*M(3,3)*M(4,2)-M(1,1)*M(2,4)*M(3,3)*M(4,2)-               &
        M(1,3)*M(2,1)*M(3,4)*M(4,2)+M(1,1)*M(2,3)*M(3,4)*M(4,2)+               &
        M(1,4)*M(2,2)*M(3,1)*M(4,3)-M(1,2)*M(2,4)*M(3,1)*M(4,3)-               &
        M(1,4)*M(2,1)*M(3,2)*M(4,3)+M(1,1)*M(2,4)*M(3,2)*M(4,3)+               &
        M(1,2)*M(2,1)*M(3,4)*M(4,3)-M(1,1)*M(2,2)*M(3,4)*M(4,3)-               &
        M(1,3)*M(2,2)*M(3,1)*M(4,4)+M(1,2)*M(2,3)*M(3,1)*M(4,4)+               &
        M(1,3)*M(2,1)*M(3,2)*M(4,4)-M(1,1)*M(2,3)*M(3,2)*M(4,4)-               &
        M(1,2)*M(2,1)*M(3,3)*M(4,4)+M(1,1)*M(2,2)*M(3,3)*M(4,4)
    vol=vol/f6
  end subroutine ElVolTet

  ! Computes 'N', the shape functions
  subroutine ShapeFuncPrecompTet(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    character(*), parameter :: eltype = "tet"
    integer :: nip, typeIndex, nodecount, i
    real(8) :: eta, nu, psi 
    
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    
    if(allocated(shapeFuncMem(typeIndex)%array)) deallocate(shapeFuncMem(typeIndex)%array)
    allocate(shapeFuncMem(typeIndex)%array(nip))
    do i = 1, nip
      allocate(shapeFuncMem(typeIndex)%array(i)%array(nodecount))
      eta = ipoints(typeIndex)%array(i, 1)
      nu = ipoints(typeIndex)%array(i, 2)
      psi = ipoints(typeIndex)%array(i, 3)
      shapeFuncMem(typeIndex)%array(i)%array(1) = f1 - eta - nu - psi
      shapeFuncMem(typeIndex)%array(i)%array(2) = eta
      shapeFuncMem(typeIndex)%array(i)%array(3) = nu
      shapeFuncMem(typeIndex)%array(i)%array(4) = psi
    end do
  end subroutine ShapeFuncPrecompTet

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta','nu' and 'psi'
  subroutine ShapeFuncdTet(dN,coord)
    implicit none
    real(8) :: dN(3,4),coord(3),e,n,s
    e=coord(1); n=coord(2); s=coord(3)
    dN(1,:)=(/-f1,f1,f0,f0/) ! dN/de
    dN(2,:)=(/-f1,f0,f1,f0/) ! dN/dn
    dN(3,:)=(/-f1,f0,f0,f1/) ! dN/ds
  end subroutine ShapeFuncdTet

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesTet(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(4),side,snodes(3),j(3)
    real(8) :: ecoords(4,3),area
    if (side==1) then; j(1)=1; j(2)=2; j(3)=4
    else if (side==2) then; j(1)=2; j(2)=3; j(3)=4
    else if (side==3) then; j(1)=1; j(2)=3; j(3)=4
    else if (side==4) then; j(1)=1; j(2)=2; j(3)=3
    end if
    snodes=enodes(j)
    call TriArea(ecoords(j(1),1),ecoords(j(1),2),ecoords(j(1),3),              &
                 ecoords(j(2),1),ecoords(j(2),2),ecoords(j(2),3),              &
                 ecoords(j(3),1),ecoords(j(3),2),ecoords(j(3),3),area)
  end subroutine EdgeAreaNodesTet

! Linear Hex ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsHex(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    integer :: nip, typeIndex
    nip = getNip("hex")
    typeIndex = getElTypeNo("hex")
    
    if(allocated(ipoints(typeIndex)%array)) deallocate(ipoints(typeIndex)%array)
    if(allocated(weights(typeIndex)%array)) deallocate(weights(typeIndex)%array)
    allocate(ipoints(typeIndex)%array(nip, pdim))
    allocate(weights(typeIndex)%array(nip))
    
    ipoints(typeIndex)%array(1,:)=(/-sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(2,:)=(/ sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(3,:)=(/ sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(4,:)=(/-sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoints(typeIndex)%array(5,:)=(/-sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    ipoints(typeIndex)%array(6,:)=(/ sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    ipoints(typeIndex)%array(7,:)=(/ sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    ipoints(typeIndex)%array(8,:)=(/-sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    weights(typeIndex)%array=f1 ! Weight = Wi = wi*wj*wk
  end subroutine SamPtsHex

  ! Computes the element volume
  ! A hex can be sub-divided into 6 tets without adding new nodes
  subroutine ElVolHex(ecoords,vol)
    implicit none
    real(8) :: ecoords(8,3),tetcoords(4,3),vol,tetvol
    vol=f0
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(3,:)
    tetcoords(3,:)=ecoords(8,:); tetcoords(4,:)=ecoords(7,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(3,:)
    tetcoords(3,:)=ecoords(7,:); tetcoords(4,:)=ecoords(4,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(5,:)
    tetcoords(3,:)=ecoords(6,:); tetcoords(4,:)=ecoords(7,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(5,:)
    tetcoords(3,:)=ecoords(7,:); tetcoords(4,:)=ecoords(8,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(2,:)
    tetcoords(3,:)=ecoords(5,:); tetcoords(4,:)=ecoords(8,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
    tetcoords(1,:)=ecoords(1,:); tetcoords(2,:)=ecoords(2,:)
    tetcoords(3,:)=ecoords(8,:); tetcoords(4,:)=ecoords(3,:)
    call ElVolTet(tetcoords,tetvol); vol=vol+abs(tetvol)
  end subroutine ElVolHex

  ! Computes 'N', the shape functions
  subroutine ShapeFuncPrecompHex(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    character(*), parameter :: eltype = "hex"
    integer :: nip, typeIndex, nodecount, i
    real(8) :: eta, nu, psi 
    real(8), parameter :: c = 0.125d0
    
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    
    if(allocated(shapeFuncMem(typeIndex)%array)) deallocate(shapeFuncMem(typeIndex)%array)
    allocate(shapeFuncMem(typeIndex)%array(nip))
    do i = 1, nip
      allocate(shapeFuncMem(typeIndex)%array(i)%array(nodecount))
      eta = ipoints(typeIndex)%array(i, 1)
      nu = ipoints(typeIndex)%array(i, 2)
      psi = ipoints(typeIndex)%array(i, 3)
      shapeFuncMem(typeIndex)%array(i)%array(1) = c*(f1-eta)*(f1-nu)*(f1-psi)
      shapeFuncMem(typeIndex)%array(i)%array(2) = c*(f1+eta)*(f1-nu)*(f1-psi)
      shapeFuncMem(typeIndex)%array(i)%array(3) = c*(f1+eta)*(f1+nu)*(f1-psi)
      shapeFuncMem(typeIndex)%array(i)%array(4) = c*(f1-eta)*(f1+nu)*(f1-psi)
      shapeFuncMem(typeIndex)%array(i)%array(5) = c*(f1-eta)*(f1-nu)*(f1+psi)
      shapeFuncMem(typeIndex)%array(i)%array(6) = c*(f1+eta)*(f1-nu)*(f1+psi)
      shapeFuncMem(typeIndex)%array(i)%array(7) = c*(f1+eta)*(f1+nu)*(f1+psi)
      shapeFuncMem(typeIndex)%array(i)%array(8) = c*(f1-eta)*(f1+nu)*(f1+psi)
    end do
  end subroutine ShapeFuncPrecompHex

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta','nu' and 'psi'
  subroutine ShapeFuncdHex(dN,coord)
    implicit none
    real(8) :: dN(3,8),coord(3),e,n,s
    real(8),parameter :: c=0.125d0
    e=coord(1); n=coord(2); s=coord(3)
    dN(1,:)=c*(/-(f1-n)*(f1-s), (f1-n)*(f1-s), (f1+n)*(f1-s),-(f1+n)*(f1-s),   &
           -(f1-n)*(f1+s), (f1-n)*(f1+s), (f1+n)*(f1+s),-(f1+n)*(f1+s)/) ! dN/de
    dN(2,:)=c*(/-(f1-e)*(f1-s),-(f1+e)*(f1-s), (f1+e)*(f1-s), (f1-e)*(f1-s),   &
           -(f1-e)*(f1+s),-(f1+e)*(f1+s), (f1+e)*(f1+s), (f1-e)*(f1+s)/) ! dN/dn
    dN(3,:)=c*(/-(f1-e)*(f1-n),-(f1+e)*(f1-n),-(f1+e)*(f1+n),-(f1-e)*(f1+n),   &
            (f1-e)*(f1-n), (f1+e)*(f1-n), (f1+e)*(f1+n), (f1-e)*(f1+n)/) ! dN/ds
  end subroutine ShapeFuncdHex

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesHex(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(8),side,snodes(4),j(4)
    real(8) :: ecoords(8,3),area
    if (side==1) then; j(1)=1; j(2)=2; j(3)=6; j(4)=5
    else if (side==2) then; j(1)=2; j(2)=3; j(3)=7; j(4)=6
    else if (side==3) then; j(1)=3; j(2)=4; j(3)=8; j(4)=7
    else if (side==4) then; j(1)=4; j(2)=1; j(3)=5; j(4)=8
    else if (side==5) then; j(1)=1; j(2)=2; j(3)=3; j(4)=4
    else if (side==6) then; j(1)=5; j(2)=6; j(3)=7; j(4)=8
    end if
    snodes=enodes(j)
    call QuadArea(ecoords(j(1),1),ecoords(j(1),2),ecoords(j(1),3),             &
                  ecoords(j(2),1),ecoords(j(2),2),ecoords(j(2),3),             &
                  ecoords(j(3),1),ecoords(j(3),2),ecoords(j(3),3),             &
                  ecoords(j(4),1),ecoords(j(4),2),ecoords(j(4),3),area)
  end subroutine EdgeAreaNodesHex    

  ! Cohesive element
    subroutine SamPtsCoh(ipoints, weights)
    implicit none
    type(array2d) :: ipoints(:)
    type(array1d) :: weights(:)
    integer :: nip, typeIndex    
    ! CN: Precomputed value
    real(8), parameter :: CN = 0.5773502691896260D0 
    
    nip = getNip("coh")
    typeIndex = getElTypeNo("coh")
    
    if(allocated(ipoints(typeIndex)%array)) deallocate(ipoints(typeIndex)%array)
    if(allocated(weights(typeIndex)%array)) deallocate(weights(typeIndex)%array)
    allocate(ipoints(typeIndex)%array(nip, pdim))
    allocate(weights(typeIndex)%array(nip))
    
    if(nip==2) then
      ipoints(typeIndex)%array(1, :) = (/-CN, f0/)
      ipoints(typeIndex)%array(2, :) = (/CN, f0/)
      weights(typeIndex)%array = f1
    end if
    
  end subroutine SamPtsCoh  
  
  subroutine ShapeFuncPrecompCoh(ipoints)
    implicit none
    type(array2d), intent(in) :: ipoints(:)
    character(*), parameter :: eltype = "coh"
    integer :: nip, typeIndex, nodecount, i
    real(8) :: eta  
    
    typeIndex = getElTypeNo(eltype)
    nip = getNip(eltype)
    nodecount = getNodeCount(eltype)
    
    if(allocated(shapeFuncMem(typeIndex)%array)) deallocate(shapeFuncMem(typeIndex)%array)
    allocate(shapeFuncMem(typeIndex)%array(nip))
    do i = 1, nip
      allocate(shapeFuncMem(typeIndex)%array(i)%array(nodecount))
      eta = ipoints(typeIndex)%array(i, 1)
      shapeFuncMem(typeIndex)%array(i)%array(1) = 0.5*(1-eta)
      shapeFuncMem(typeIndex)%array(i)%array(2) = 0.5*(1+eta)
      shapeFuncMem(typeIndex)%array(i)%array(3) = shapeFuncMem(typeIndex)%array(i)%array(2)
      shapeFuncMem(typeIndex)%array(i)%array(4) = shapeFuncMem(typeIndex)%array(i)%array(1)
    end do
  end subroutine ShapeFuncPrecompCoh

  subroutine ShapeFuncdCoh(dN,coord)
    implicit none
    real(8) :: dN(2,4),coord(2)
    real(8), parameter :: half=0.5
    dN(1,:)=(/-half, half, half, -half/) ! dN/de
    dN(2,:)=dN(1,:) ! dN/dn
  end subroutine ShapeFuncdCoh
  
  ! Get cohesive element's tangent, det and normal vector
  subroutine getCohValues(ecoords, tangent, normal, det)
    implicit none
    real(8), intent(in) :: ecoords(:, :)
    real(8), intent(out) :: tangent(pdim), normal(pdim), det
    real(8), parameter :: HALF=0.5d0
    
    tangent(1)=HALF*(ecoords(2,1)-ecoords(1,1)+   &
                     ecoords(3,1)-ecoords(4,1))
    tangent(2)=HALF*(ecoords(2,2)-ecoords(1,2)+   &
                     ecoords(3,2)-ecoords(4,2))
    call magnitude(tangent, det)
    tangent=tangent/det
    det=det*HALF
    normal(1)=-tangent(2)
    normal(2)=tangent(1)        
  end subroutine getCohValues
  
  ! Get cohesive's relative displacement and velocity
  subroutine getCohRels(nodes, local_u, nip, dt, urel, vrel)
    implicit none
    integer, intent(in) :: nodes(:), nip
    ! Note at local_u, node n's dofs: 2n-1, 2n
    real(8), intent(in) :: local_u(:), dt
    real(8), intent(out) :: urel(pdim), vrel(pdim)
    
    ! Variables
    ! N: ShapeFunc
    integer :: i, node1, node2
    real(8), pointer :: N(:)
    integer :: nodecount
    
    ! Constants
    real(8), parameter :: ZERO=0.0d0
    
    call ShapeFunc("coh", nip, N)
    nodecount = getNodeCount("coh")
    urel = ZERO
    vrel = ZERO
    
    do node1=1,nodecount/2
      node2=node1+nodecount/2
      do i=1,pdim
        urel(i) =   urel(i)                           &
                  + N(node2) * local_u(2*(node2-1)+i) &
                  - N(node1) * local_u(2*(node1-1)+i)
      end do
    end do
    vrel = urel / dt
  end subroutine getCohRels
  
  ! Form Gap and VGap for cohesive
  subroutine getCohGaps(tangent, normal, urel, vrel, gap, vgap)
    implicit none
    real(8), intent(in) :: tangent(pdim), normal(pdim), urel(pdim), vrel(pdim)
    real(8), intent(out) :: gap(pdim), vgap(pdim)
    
    gap(1) = normal(1) * urel(1) + normal(2) * urel(2)
    gap(2) = tangent(1) * urel(1) + tangent(2) * urel(2)
    
    vgap(1) = normal(1) * vrel(1) + normal(2) * vrel(2)
    vgap(2) = tangent(1) * vrel(1) + tangent(2) * vrel(2)
    
  end subroutine getCohGaps

  ! Memoize the nodal stress inverse matrix
  ! ShapeFuncMem must be initialized
  ! TODO: Handle coh
  subroutine NodalStressInv
    implicit none
    character(4) :: eltype
    integer :: nip, nodecount, i, j
    real(8), allocatable :: N2(:, :)
    
    do i = 1, elTypeCount
      eltype = eltypes(i)
      nip = getNip(eltype)
      nodeCount = getNodeCount(eltype)
      if (getDim(eltype) == pdim .and. nip == nodecount) then
        if(allocated(nodalStressInvMem(i)%array)) deallocate(nodalStressInvMem(i)%array)
        allocate(N2(nodecount, nip), nodalStressInvMem(i)%array(nip, nodecount))
        do j = 1, nip
          N2(j, :) = ShapeFuncMem(i)%array(j)%array
        end do
        call MatInv(N2, nodalStressInvMem(i)%array)
        deallocate(N2)
      end if
    end do
  end subroutine NodalStressInv

end module elems
