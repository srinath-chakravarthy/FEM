#define ALP_PC

module global
  use local
  use seplaw
  use elems
  implicit none
#if defined ALP_PC
#include <finclude/petscdef.h>
#else
#include <petsc-finclude/petscdef.h>
#endif
  ! Global constants
  integer, parameter :: BC_PRESENT=0
  real(8), parameter :: PENALTY_PARAM=1.0e30
  ! Global variables
  integer :: nnds,nels,nmts,nceqs,nfrcs,ntrcs,nbcs,steps,maxnodesperel
  real(8) :: alpha,beta,t,dt,val,scale
  integer,allocatable :: nodes(:, :), work(:)
  real(8),allocatable :: vvec(:),mat(:,:)
  type(array2d), allocatable :: stress(:)
  real(8),allocatable,target :: uu(:),tot_uu(:),uup(:),coords(:,:)
  character(12) :: stype
  Vec :: Vec_F,Vec_U,Vec_Um,Vec_Up,Vec_lambda,Vec_I
  Vec,pointer :: Vec_W(:),Vec_Wlm(:)
  Mat :: Mat_K,Mat_M,Mat_Minv,Mat_Gt,Mat_G,Mat_GMinvGt,Mat_Kc, PC_MAT
  KSP :: Krylov
  SNES :: Solver
  PC  :: PreCon
  PetscInt :: ival, icntl, infog34
  PetscReal :: cntl, rinfo12, rinfo13, valmum, tol
  ! Variables to store loading information
  integer,allocatable :: fnode(:),telsd(:,:), bcnode(:, :)
  real(8),allocatable :: cval(:,:),fval(:,:),tval(:,:), bcval(:,:)
  ! Local element/side/node variables
  integer :: el,side,node
  real(8) :: E,nu,dns,visc,expn,H,B,phi,Kf,area,vol
  integer,allocatable :: indx(:),indxp(:)
  real(8),allocatable :: k(:,:),kc(:,:),Hs(:)
  type(element), allocatable, target :: global_elements(:), local_elements(:)
  ! Variables for parallel code
  integer :: nprcs,rank,ierr
  integer,allocatable :: epart(:),npart(:) ! Partitioning
  integer,allocatable :: nmap(:),emap(:),nl2g(:),indxmap(:,:) ! Mappings
  Vec :: Seq_U
  IS  :: From,To,RI
  VecScatter :: Scatter
  real(8),pointer :: pntr(:)
  
  ! New variables for node stress
  
  ! Dimensions:
  !         count_node(nlnds)
  !         stress_node(nlnds, cpdim)
  !         stress_at_el(local_elements(i)%nodecount, cpdim)
  
  integer, allocatable :: count_node(:) ! Contribution count to stress at given node
  integer :: nlnds ! Number of local nodes 
  real(8), allocatable, target :: stress_node(:, :) ! Aggregate stress at local nodes
  real(8), allocatable :: stress_at_el(:, :) ! Stress at the element's nodes
  Mat :: Mat_Stress, Seq_Stress ! Matrices of node stress
  integer, allocatable :: nodeIndex(:) ! Local to global node index mapping
  integer, allocatable :: global_count_node(:) ! Global contribution count
  IS :: Stress_row, Cpdim_Iteration ! Index sets for sequential Mat_Stress extraction
  real(8), allocatable, target :: aggregate_u(:), aggregate_stress(:, :) ! Running totals
  
  ! Iterative soluton
  real(8) :: t_init ! Start of the current interval
  integer :: dtNo ! The interval no
  
  ! Nonlinear solution
  Vec :: Residual ! SNES residual
  Mat :: Jacobian ! Jacobian
  integer :: iterationCount ! SNES KSP iteration count
  integer, allocatable :: nonlins(:) ! Local indices of nonlinear elements
  integer :: l_nonlin_ec ! Local nonlinear el count
  integer, allocatable :: nonlin_els(:) ! Indices of linear elements
  
  ! Extractions
  real(8), allocatable :: Array_K(:, :) ! Local Mat_K
  
  ! Bandwidth guess
  integer :: nodal_bw

  ! Cohesive materials
  integer :: ncohmats
  type(cohMat), allocatable :: cohmats(:)

  interface ApplyNodalForce
    module procedure ApplyNodalForce_All, ApplyNodalForce_BC
  end interface ApplyNodalForce

contains

! Test props: 100.0_8, 0.01_8, 0.01_8, 1.0_8, 0.0_8, 0.0_8

  ! Calculate Jacobian for SNES
  subroutine CalcJacobian(Solver, du, Jacobian, PCon, args, ierr) 
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    ! Input
    SNES :: Solver
    Vec :: du
    Mat :: Jacobian, PCon
    integer(8), optional :: args(*)
    PetscErrorCode :: ierr
    
    integer :: i, j, k
    integer :: node1, node2, dof1, dof2
    real(8), allocatable :: appliedStiff(:, :)

    ! Convenience
    type(element), pointer :: el

    
    ! Constants
    real(8), parameter :: HALF = 0.5
    
    call MatCopy(Mat_K, Jacobian,SAME_NONZERO_PATTERN,ierr)
!    call PrintMsg("Elastic Matrix")
!    call MatView(Jacobian, Petsc_Viewer_Stdout_World, ierr)
    call GetVec_U(du, uu)
    ! uu is the total displacement at the given time
    uu = uu + aggregate_u
    
    do i = 1, nels
      el => local_elements(i)
      if(allocated(appliedStiff)) deallocate(appliedStiff)
      allocate(appliedStiff(pdim * el%nodeCount, pdim * el%nodeCount))

      if(el%nlMat == 0) then
        cycle
      else if (local_elements(i)%nlMat == 1) then
        call applyStiff_1(local_elements(i), uu, dt, appliedStiff)
      else
        call PrintMsg("Invalid material")
      endif

      do node1 = 1, el%nodeCount; do node2 = 1, el%nodeCount; do dof1 = 1, pdim; do dof2 = 1, pdim
      call MatSetValue(Jacobian, &
                       pdim*(nl2g(el%nodes(node1))-1)+dof1-1, &
                       pdim*(nl2g(el%nodes(node2))-1)+dof2-1, & 
                       appliedStiff((node1 - 1) * pdim + dof1, (node2 - 1) * pdim + dof2), &
                       ADD_VALUES, ierr)
      enddo; enddo; enddo; enddo
    end do
    
    call MatAssemblyBegin(Jacobian,Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(Jacobian, Mat_Final_Assembly,ierr)
    call ApplyKBC(Jacobian)
    call MatAssemblyBegin(Jacobian,Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(Jacobian, Mat_Final_Assembly,ierr)
    !!$call MatView(Jacobian, Petsc_Viewer_Stdout_World, ierr)
  end subroutine

  ! Calculate Residual for SNES
  subroutine CalcResidual(Solver, du, Residual, args, ierr)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    SNES :: Solver
    Vec :: du, Residual
    integer :: args(*)
    PetscErrorCode :: ierr
    integer :: i, j
    real(8) :: current_u(pdim*nlnds)
    Vec :: AppliedForce
    real(8), allocatable :: nodeTraction(:)

    ! Convenience
    type(element), pointer :: el

    ! Constants
    real(8), parameter :: MINUS = -1.0d0 
    
    call GetVec_U(du, current_u)
    current_u = current_u + aggregate_u
    
    call VecCopy(Vec_F, Residual, ierr)
    ! For each element
    do i = 1, nels
      el => local_elements(i)
      if(allocated(nodeTraction)) deallocate(nodeTraction)
      allocate(nodeTraction(pdim * el%nodeCount))

      ! Get force
      if(el%nlMat == 0) then
        cycle
      else if (el%nlMat == 1) then 
        call applyTract_1(el, current_u, dt, nodeTraction)
      else 
        call PrintMsg("Invalid nonlinear material")
      endif
      
      ! Apply force
      do j = 1, local_elements(i)%nodeCount
        call ApplyNodalForce(Residual, nl2g(el%nodes(j)), nodeTraction((j - 1) * pdim + 1:(j - 1) * pdim + 2), .true.)
      enddo
    end do
    
    call VecAssemblyBegin(Residual, ierr)
    call VecAssemblyEnd(Residual, ierr)
    call EnforceBCForce(Residual, dt)
    call VecAssemblyBegin(Residual, ierr)
    call VecAssemblyEnd(Residual, ierr)
    
    ! --- Apply Displacement BC's to Residual here
    call EnforceBCForce(Residual,dt)

    call VecAssemblyBegin(Residual, ierr)
    call VecAssemblyEnd(Residual, ierr)
    
    
    !!$call Printmsg("Viewing the Force Vector")
    !!$call VecView(Residual, PETSC_VIEWER_STDOUT_WORLD, ierr)
    
    
    call VecScale(Residual, MINUS, ierr)
    call MatMultAdd(Jacobian, du, Residual, Residual, ierr)
    
    call VecAssemblyBegin(Residual, ierr)
    call VecAssemblyEnd(Residual, ierr)

    !!$call Printmsg("Viewing the Residual Vector")
    !!$call VecView(Residual, PETSC_VIEWER_STDOUT_WORLD, ierr)
    
    
  end subroutine CalcResidual

  ! Form local [K]

  subroutine FormLocalK(el,k,indx)
    implicit none
    integer :: el,indx(:)
    integer, pointer :: enodes(:)
    real(8) :: k(:,:), E, nu
    real(8), pointer :: ecoords(:, :)
    character(2) :: strng    

    enodes=>local_elements(el)%nodes
    ecoords=>local_elements(el)%ecoords

    if (local_elements(el)%mat == 0) then
      E = 0.0
      nu = 0.0
    else 
      E=mat(local_elements(el)%mat,1)
      nu=mat(local_elements(el)%mat,2)
    end if

    call FormElKE(local_elements(el)%eltype, ecoords, E, nu, dt, k)
    call FormLocalIndx(enodes,indx)    

  end subroutine FormLocalK

  ! Apply BCs to given matrix
  subroutine ApplyKBC(Matrix)
    implicit none
    Mat :: Matrix
    integer :: i, j, node
    real(8) :: vvec(pdim)
    !!$ ! Apply displacement boundary conditions 
    do i=1,nbcs
       vvec = 0.0
       node=bcnode(i, 1)
       do j = 1, pdim
          ! TODO: Redundant vvec composition; vvec can be equal to penalty_param
          !       Clean this with force and bc parallelization
          if(bcnode(i, j + 1) == BC_PRESENT) vvec(j) = penalty_param
       end do
       call UpdateKMatBcs(Matrix, node,bcnode(i,2:pdim+1),vvec)
    end do
  end subroutine ApplyKBC
!!$ ! Apply displacement bc's to correct row,column of stiffness matrix
  subroutine UpdateKMatBcs(Matrix, node,bc_presence,vvec)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    Mat :: Matrix
    integer :: node, i,j, bc_presence(pdim+1)
    real(8) :: vvec(:)
    do i = 1, pdim
       if(bc_presence(i) == BC_PRESENT) then
         j=(pdim)*node - (pdim)+i-1
         !!$print *, 'Applying BCs at dof ', node, j
         call MatSetValue(Matrix,j,j,vvec(i),Insert_values,ierr)
       end if
    end do
  end subroutine UpdateKMatBcs
 
  ! Apply nodal force
  subroutine ApplyNodalForce_All(Vec_F, node,vvec,iadd)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    Vec, intent(out) :: Vec_F
    integer :: node,i,j
    real(8) :: vvec(:)
    logical :: iadd
    do i=1,pdim
       j=(pdim)*node-(pdim)+i-1
       ! if (i==pdim+1) vvec(i)=scale*vvec(i)
       if (iadd) then 
          call VecSetValue(Vec_F,j,vvec(i),Add_Values,ierr)
       else
          call VecSetValue(Vec_F,j,vvec(i),Insert_Values,ierr)
       end if
    end do
  end subroutine ApplyNodalForce_All
  ! Apply BC's nodal force
  subroutine ApplyNodalForce_BC(Vec_F, node, bcs, vvec)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    Vec, intent(out) :: Vec_F
    integer :: node,bcs(pdim),i,j
    real(8) :: vvec(:)
    do i = 1, pdim
      if(bcs(i) /= BC_PRESENT) cycle
      j = pdim*node - pdim + i - 1
      call VecSetValue(Vec_F, j, vvec(i), Insert_Values, ierr)
    end do
  end subroutine ApplyNodalForce_BC
  
  
  ! Apply traction (EbEAve)
  subroutine ApplyTraction(el,side,vvec,iadd)
    implicit none
    type(element), target :: el
    integer :: side,i, nps
    integer, pointer :: enodes(:)
    real(8) :: vvec(:),area
    logical :: iadd
    real(8), pointer :: ecoords(:, :)
    integer, allocatable :: snodes(:)


    enodes=>el%nodes
    ecoords=>el%ecoords

    nps = getNps(el%eltype)
    allocate(snodes(nps))

    call EdgeAreaNodes(el%eltype, enodes,ecoords,side,area,snodes)
    vvec=vvec*area/dble(nps)
    snodes=nl2g(snodes)
    do i=1,nps
       call ApplyNodalForce(Vec_F, snodes(i),vvec,iadd)
    end do

    deallocate(snodes)
  end subroutine ApplyTraction

  ! Form RHS
  ! Edit: Removed Mat_K modification. FormRHS accesses only RHS
  ! Permit not used yet
  subroutine FormRHS(t_init, dt)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    real(8), intent(in) :: t_init, dt
    real(8) :: t_end, applied_interval
    integer :: i,j
    real(8) :: t1,t2
    logical :: iadd
    real(8) :: vvec(pdim)
    call VecZeroEntries(Vec_F, ierr)
    t_end = t_init + dt
    
    ! TODO: Modify nceqs for dt
    do i=1,nceqs
       t1=cval(i,2)/dt; t2=cval(i,3)/dt
       if(.true.) then! if (tstep>=nint(t1) .and. tstep<=nint(t2)) then
          j=(pdim)*nnds+i-1
          if (stype/="explicit" .and. rank==0)                                 &
             call VecSetValue(Vec_F,j,cval(i,1),Add_Values,ierr)
       end if
    end do

    ! Modified for time
    do i=1,nfrcs
       !t1=fval(i,pdim+1)/dt; t2=fval(i,pdim+2)/dt
       t1 = fval(i, pdim + 1)
       t2 = fval(i, pdim + 2)
       if (t_end < t1 .or. t_init > t2) cycle
       applied_interval = min(t2, t_end) - max(t1, t_init)
       node=fnode(i); vvec=fval(i,1:pdim)
       vvec = vvec * applied_interval / (t2 - t1)
       iadd = .true.
       call ApplyNodalForce(Vec_F, node,vvec,iadd)

    end do
    
    do i=1,ntrcs
      t1=tval(i,pdim+1)/dt
      t2=tval(i,pdim+2)/dt

      if(t_end < t1 .or. t_init > t2) cycle
      
      applied_interval = min(t2, t_end) - max(t1, t_init)
      el=telsd(i,1)
      side=telsd(i,2)
      vvec=tval(i,1:pdim) * (applied_interval / (t2 - t1))
      iadd = .true.
      iadd = .true.
      if (el/=0) call ApplyTraction(local_elements(el),side,vvec, iadd)
    end do

    call VecAssemblyBegin(Vec_F, ierr)
    call VecAssemblyEnd(Vec_F, ierr)

!!$ ! Apply displacement boundary conditions
    call EnforceBCForce(Vec_F, dt)

    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine FormRHS

  ! Override BC's forces onto the Vec_F
  subroutine EnforceBCForce(Vec_F, dt)
    implicit none
    Vec, intent(out) :: Vec_F
    real(8), intent(in) :: dt
    
    real(8) :: vvec(pdim)
    integer :: i, j
    
    do i=1,nbcs
      vvec = 0.0
      node=bcnode(i, 1)
      do j=1,pdim
        if(bcnode(i,j+1)==0) vvec(j)= PENALTY_PARAM * bcval(i,j) * dt / t
      end do
      call ApplyNodalForce(Vec_F, node,bcnode(i,2:pdim+1),vvec)
    end do
  end subroutine EnforceBCForce

  ! Form local index
  subroutine FormLocalIndx(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalIndx

  ! Recover stress
  ! Note: removed parameter stress
  subroutine RecoverStress(el)
    implicit none
    integer :: el, nodecount
    ! Dimensions:
    !         indx(nodecount * pdim)
    integer, allocatable :: indx(:)
    integer, pointer :: enodes(:)
    real(8), pointer :: ecoords(:, :)
    
    nodecount = local_elements(el)%nodecount
    allocate(indx(nodecount * pdim))
    
    enodes=>local_elements(el)%nodes
    ecoords=>local_elements(el)%ecoords
    E=mat(local_elements(el)%mat,1); nu=mat(local_elements(el)%mat,2)
    call FormLocalIndx(enodes,indx)
     
    call CalcElStress(local_elements(el)%eltype, ecoords,uu(indx),E,nu,stress(el)%array)

    deallocate(indx)
  end subroutine RecoverStress

  subroutine RecoverNodalStress(el, stress_at_el)
    implicit none
    ! Dimensions
    !         stress_at_el(nodecount, cpdim)
    !         N(nodecount)
    !         N2(nip, nodecount)    ; nip == nodecount
    !         N2Inv(nodeCount, nip) ; nip == nodecount
    integer, intent(in) :: el
    real(8), intent(out) :: stress_at_el(:, :)
    
    character(4) :: eltype
    integer :: nip, nodecount, typeIndex, i
    
    nodecount = local_elements(el)%nodecount
    eltype = local_elements(el)%eltype
    nip = getNip(eltype)
    typeIndex = getElTypeNo(eltype)
    
    if (nip == nodecount) then
      do i = 1, cpdim
        stress_at_el(:, i) = matmul(nodalStressInvMem(typeIndex)%array, stress(el)%array(:, i))
      end do
    else if (nip == 1) then
      do i = 1, nodecount
        stress_at_el(i, :) = stress(el)%array(1, :)
      end do     
    end if
  end subroutine RecoverNodalStress  

  
  ! Reform RHS
  ! UNUSED, WAS NOT EDITED
  subroutine ReformLocalRHS(el,f,indx)
    implicit none
    integer :: el,indx(:),j,j1
    integer, pointer :: enodes(:)
    real(8) :: f(:)
    real(8), pointer :: ecoords(:, :)
    enodes=>local_elements(el)%nodes
    ecoords=>local_elements(el)%ecoords
    E=mat(local_elements(el)%mat,1); nu=mat(local_elements(el)%mat,2)
    visc=mat(local_elements(el)%mat,3); expn=mat(local_elements(el)%mat,4)
    f=f0
    call ReformElRHS(ecoords,stress(el)%array,E,nu,visc,expn,dt,f(1:eldof))
    call FormLocalIndx(enodes,indx)
    ! Fix BCs (i.e., zero out entries)
    do j=1,nodesperel
       do j1=1,pdim
          ! if (bc(nodes(el,j),j1)==0) f(pdim*j-pdim+j1)=f0
       end do
    end do
  end subroutine ReformLocalRHS

  ! Print message
  subroutine PrintMsg(msg)
    implicit none
    character(*) :: msg
    if (rank==0) print*,msg
  end subroutine PrintMsg

  ! Scatter given U and get all local values
  subroutine GetVec_U(Vec_U, u)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    Vec :: Vec_U
    real(8) :: u(nlnds * pdim)
    real(8), pointer :: pntr(:)
    call VecScatterBegin(Scatter,Vec_U,Seq_U,Insert_Values,Scatter_Forward,ierr)
    call VecScatterEnd(Scatter,Vec_U,Seq_U,Insert_Values,Scatter_Forward,ierr)
    call VecGetArrayF90(Seq_U,pntr,ierr)
    u=pntr
    call VecRestoreArrayF90(Seq_U,pntr,ierr)
  end subroutine GetVec_U

  subroutine applyTract_1_Object(el, du, dt, Vec_F)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    real(8) :: du(:), dt
    Vec :: Vec_F
    type(element), target :: el

    real(8) :: coh_tangent(pdim), coh_norm(pdim), coh_gap, det
    real(8) :: urel(pdim), vrel(pdim), gap(pdim), vgap(pdim)
    real(8) :: force_coh(pdim)
    real(8), pointer :: N1(:)
    
    character(4) :: eltype
    integer :: nip, nodecount, typeIndex
    
    real(8) sig
    integer :: node1, dof
    real(8), pointer :: ecoords(:,:)
    integer, pointer :: nodes(:)
    real(8) :: added_values(pdim) 
    integer :: intp, i, j, k

    eltype = el%eltype  
    typeIndex = getElTypeNo(eltype) 
    nip = getNip(eltype)
    nodecount = el%nodecount
    
    nodes=>el%nodes
    ecoords=>el%ecoords
    call getCohValues(ecoords, coh_tangent, coh_norm, det)
    do intp = 1, nip
      call getCohRels(nodes, du, intp, dt, urel, vrel)
      call getCohGaps(coh_tangent, coh_norm, urel, vrel, gap, vgap)
      call Seplaw_1_Tract(cohMats(el%nlMat)%props, &
                          gap, vgap, dt, force_coh)
      call ShapeFunc(eltype, intp, N1)
      
      sig = 1.0d0
      do node = 1, nodecount
        if(node .gt. nodecount/2) sig = -1.0d0
        added_values = 0.0
        do dof = 1, pdim
          added_values(dof) = sig*N1(node)* &
                              (force_coh(1)*coh_norm(dof)+force_coh(2)*coh_tangent(dof))  &
                              *weights(typeIndex)%array(intp)*det
        end do
        call ApplyNodalForce(Vec_F, nl2g(el%nodes(node)), added_values, .true.)
      end do
    end do

  end subroutine 

  subroutine applyTract_1(el, du, dt, result)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    real(8) :: du(:), dt
    type(element), target :: el
    !Vec :: Vec_F
    real(8) :: result(pdim * el%nodeCount)

    real(8) :: coh_tangent(pdim), coh_norm(pdim), coh_gap, det
    real(8) :: urel(pdim), vrel(pdim), gap(pdim), vgap(pdim)
    real(8) :: force_coh(pdim)
    real(8), pointer :: N1(:)
    
    character(4) :: eltype
    integer :: nip, nodecount, typeIndex
    
    real(8) sig
    integer :: node1, dof
    real(8), pointer :: ecoords(:,:)
    integer, pointer :: nodes(:)
    !real(8) :: added_values(pdim) 
    integer :: intp, i, j, k

    eltype = el%eltype  
    typeIndex = getElTypeNo(eltype) 
    nip = getNip(eltype)
    nodecount = el%nodecount
    
    nodes=>el%nodes
    ecoords=>el%ecoords
    call getCohValues(ecoords, coh_tangent, coh_norm, det)
    do intp = 1, nip
      call getCohRels(nodes, du, intp, dt, urel, vrel)
      call getCohGaps(coh_tangent, coh_norm, urel, vrel, gap, vgap)
      call Seplaw_1_Tract(cohMats(el%nlMat)%props, &
                          gap, vgap, dt, force_coh)
      call ShapeFunc(eltype, intp, N1)
      
      sig = 1.0d0
      do node = 1, nodecount
        if(node .gt. nodecount/2) sig = -1.0d0
        !added_values = 0.0
        do dof = 1, pdim
          result((node - 1) * pdim + dof) = sig*N1(node)* &
                              (force_coh(1)*coh_norm(dof)+force_coh(2)*coh_tangent(dof))  &
                              *weights(typeIndex)%array(intp)*det 
!          added_values(dof) = sig*N1(node)* &
 !                             (force_coh(1)*coh_norm(dof)+force_coh(2)*coh_tangent(dof))  &
  !                            *weights(typeIndex)%array(intp)*det
        end do
        !result((node - 1) * pdim + 1:(node - 1) * pdim + 2) = added_values
        !print*, result
        !call ApplyNodalForce(Vec_F, nl2g(el%nodes(node)), added_values, .true.)
      end do
    end do

  end subroutine !applyTract_1

  subroutine applyStiff_1_Obj(el, du, dt, Stiff_Mat)
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    real(8), intent(in) :: du(:), dt
    Mat, intent(inout) :: Stiff_Mat
    type(element), target :: el

    ! Variables
    ! det: Tangent's magnitude / 2
    ! uu: current total magnitude
    ! coh_tangent: tanget of the cohesive
    ! coh_norm: normal to the cohesive
    ! coh_gap: 
    ! urel: relative displacement of the cohesive
    ! vrel: relative velocity
    ! stiff_u: single integration point's stiffness contribution
    ! N1, N2: Shape funcs for cohesive
    integer :: i, j, k
    integer :: dof1, dof2
    real(8) :: det 
    real(8) :: coh_tangent(pdim), coh_norm(pdim), coh_gap
    real(8) :: urel(pdim), vrel(pdim), gap(pdim), vgap(pdim)
    real(8) :: stiff_coh(pdim, pdim)
    real(8), pointer :: N1(:)
    
    integer :: node1, node2
    real(8), pointer :: ecoords(:,:)
    integer, pointer :: nodes(:)
    integer :: nodecount, typeIndex, nip
    real(8) :: sig1, sig2 ! Sign of the operation
    
    character(4) :: eltype

    ! Jacobian resetter
    real(8) :: elasticVals(pdim, pdim)
    
    ! Constants
    real(8), parameter :: HALF = 0.5

      eltype = el%eltype

      nip=getNip(eltype)
      nodecount = el%nodecount
      ecoords=>el%ecoords
      nodes=>el%nodes
      call getCohValues(ecoords, coh_tangent, coh_norm, det)

      do j=1,nip
        call getCohRels(nodes, du, j, dt, urel, vrel)
        call getCohGaps(coh_tangent, coh_norm, urel, vrel, gap, vgap)
        call Seplaw_1_Stiff(cohMats(el%nlMat)%props, &
                            gap, vgap, dt, stiff_coh)
        sig1 = 1.0d0
        do node1=1, nodecount
          if(node1 .gt. nodecount/2) sig1 = -1.0d0
          sig2 = 1.0d0
          do node2=1, nodecount
            if(node2 .gt. nodecount/2) sig2 = -1.0d0
            call ShapeFunc(el%eltype,j,N1)          
            do dof1=1,PDIM
              do dof2=1,PDIM
                call MatSetValue(Stiff_Mat, &
                                 pdim*(nl2g(nodes(node1))-1)+dof1-1, &
                                 pdim*(nl2g(nodes(node2))-1)+dof2-1, & 
                                 ((stiff_coh(1,1)*coh_norm(dof1)+stiff_coh(2,1)*coh_tangent(dof1))*coh_norm(dof2)) + &
                                 ((stiff_coh(1,2)*coh_norm(dof1)+stiff_coh(2,2)*coh_tangent(dof1))*coh_tangent(dof2))*sig1*sig2*N1(node1)*N1(node2)*weights(getElTypeNo(eltype))%array(j)*det, &
                                 ADD_VALUES, ierr)
              end do
            end do
          end do
        end do
      end do

  end subroutine applyStiff_1_Obj

  subroutine applyStiff_1(el, du, dt, result)
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    real(8), intent(in) :: du(:), dt
    !Mat, intent(inout) :: Stiff_Mat
    type(element), target :: el
    real(8) :: result(pdim * el%nodeCount, pdim * el%nodeCount)

    ! Variables
    ! det: Tangent's magnitude / 2
    ! uu: current total magnitude
    ! coh_tangent: tanget of the cohesive
    ! coh_norm: normal to the cohesive
    ! coh_gap: 
    ! urel: relative displacement of the cohesive
    ! vrel: relative velocity
    ! stiff_u: single integration point's stiffness contribution
    ! N1, N2: Shape funcs for cohesive
    integer :: i, j, k
    integer :: dof1, dof2
    real(8) :: det 
    real(8) :: coh_tangent(pdim), coh_norm(pdim), coh_gap
    real(8) :: urel(pdim), vrel(pdim), gap(pdim), vgap(pdim)
    real(8) :: stiff_coh(pdim, pdim)
    real(8), pointer :: N1(:)
    
    integer :: node1, node2
    real(8), pointer :: ecoords(:,:)
    integer, pointer :: nodes(:)
    integer :: nodecount, typeIndex, nip
    real(8) :: sig1, sig2 ! Sign of the operation
    
    character(4) :: eltype

    ! Jacobian resetter
    real(8) :: elasticVals(pdim, pdim)
    
    ! Constants
    real(8), parameter :: HALF = 0.5

      eltype = el%eltype

      nip=getNip(eltype)
      nodecount = el%nodecount
      ecoords=>el%ecoords
      nodes=>el%nodes
      call getCohValues(ecoords, coh_tangent, coh_norm, det)

      result = 0

      do j=1,nip
        call getCohRels(nodes, du, j, dt, urel, vrel)
        call getCohGaps(coh_tangent, coh_norm, urel, vrel, gap, vgap)
        call Seplaw_1_Stiff(cohMats(el%nlMat)%props, &
                            gap, vgap, dt, stiff_coh)
        sig1 = 1.0d0
        do node1=1, nodecount
          if(node1 .gt. nodecount/2) sig1 = -1.0d0
          sig2 = 1.0d0
          do node2=1, nodecount
            if(node2 .gt. nodecount/2) sig2 = -1.0d0
            call ShapeFunc(el%eltype,j,N1)          
            do dof1=1,PDIM
              do dof2=1,PDIM
                result((node1 - 1 ) * pdim + dof1, (node2 - 1) * pdim + dof2) = &
                    result((node1 - 1 ) * pdim + dof1, (node2 - 1) * pdim + dof2) + ((stiff_coh(1,1)*coh_norm(dof1)+stiff_coh(2,1)*coh_tangent(dof1))*coh_norm(dof2)) + &
                    ((stiff_coh(1,2)*coh_norm(dof1)+stiff_coh(2,2)*coh_tangent(dof1))*coh_tangent(dof2))*sig1*sig2*N1(node1)*N1(node2)*weights(getElTypeNo(eltype))%array(j)*det
!                call MatSetValue(Stiff_Mat, &
 !                                pdim*(nl2g(nodes(node1))-1)+dof1-1, &
  !                               pdim*(nl2g(nodes(node2))-1)+dof2-1, & 
   !                              ((stiff_coh(1,1)*coh_norm(dof1)+stiff_coh(2,1)*coh_tangent(dof1))*coh_norm(dof2)) + &
    !                             ((stiff_coh(1,2)*coh_norm(dof1)+stiff_coh(2,2)*coh_tangent(dof1))*coh_tangent(dof2))*sig1*sig2*N1(node1)*N1(node2)*weights(getElTypeNo(eltype))%array(j)*det, &
     !                            ADD_VALUES, ierr)
              end do
            end do
          end do
        end do
      end do


  end subroutine ! applyStiff_1

end module global
