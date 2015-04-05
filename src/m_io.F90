#define ALP_PC
#define destroy(x) if(allocated(x)) deallocate(x)

module io
  use global
  implicit none 
  
  contains


  ! Read simulation parameters
  subroutine ReadParameters(file)
    implicit none
    integer, intent(in) :: file

    read(file,*)stype,pdim,nodal_bw
    read(file,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nbcs
    read(file,*)t,dt
    
  end subroutine ReadParameters

  subroutine BroadCastParameters(rank, broadcaster)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: rank, broadcaster

    ! Constants
    integer, parameter :: INT_BUFFER_SIZE = 10, REAL_BUFFER_SIZE = 2

    ! Buffers
    integer :: int_buffer(INT_BUFFER_SIZE)
    real(8) :: real_buffer(REAL_BUFFER_SIZE)
    integer :: ierr

    if(rank==broadcaster) then
      if(stype == "implicit") then
        int_buffer(1) = 0
      else
        int_buffer(1) = 1
      end if
      int_buffer(2:) = (/pdim, nodal_bw, nels, nnds, nmts, nceqs, nfrcs, ntrcs, nbcs/)
      real_buffer = (/t, dt/)
    end if

    call MPI_Bcast(int_buffer, INT_BUFFER_SIZE, &
                   MPI_INTEGER, broadcaster, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(real_buffer, REAL_BUFFER_SIZE,&
                   MPI_REAL8, broadcaster, MPI_COMM_WORLD, ierr)
    if (rank /= broadcaster) then
        if(int_buffer(1) == 0) stype = "implicit"
        if(int_buffer(1) == 1) stype = "explicit"
        pdim = int_buffer(2)
        nodal_bw = int_buffer(3)
        nels = int_buffer(4)
        nnds = int_buffer(5)
        nmts = int_buffer(6)
        nceqs = int_buffer(7)
        nfrcs = int_buffer(8)
        ntrcs = int_buffer(9)
        nbcs = int_buffer(10)
        t = real_buffer(1)
        dt = real_buffer(2)
    end if 
  end subroutine BroadCastParameters

  subroutine ReadElementsCoords(file, global_elements, aggregatenode, nonlinear)
    implicit none
    integer, intent(in) :: file
    type(element), allocatable, intent(out) :: global_elements(:)
    integer, intent(out) :: aggregatenode
    logical, intent(out) :: nonlinear

    integer :: i
    character(4) :: dummy
    real(8), allocatable :: coords(:, :)

    destroy(global_elements)
    allocate(global_elements(nels))

    aggregatenode = 0
    do i=1,nels  
      ! aggregatenode accrues total number of nodes per element for
      ! all elements
      read(file, *)global_elements(i)%eltype
      if(global_elements(i)%eltype == "coh") nonlinear = .true.
      backspace(file)
      global_elements(i)%nodecount = getNodeCount(global_elements(i)%eltype)
      allocate(global_elements(i)%nodes(global_elements(i)%nodecount))
      allocate(global_elements(i)%ecoords(global_elements(i)%nodecount, pdim))
      aggregatenode= aggregatenode+(global_elements(i)%nodecount)
      read(file, *)dummy, global_elements(i)%nodes, global_elements(i)%mat
    end do

    allocate(coords(nnds, pdim))
    do i=1,nnds
      read(file, *)coords(i, :)
    end do
    do i=1,nels
      global_elements(i)%ecoords=coords(global_elements(i)%nodes, :)
    end do
    deallocate(coords)
  end subroutine ReadElementsCoords

  subroutine PartitionBroadcast(rank, broadcaster, global_elements, aggregatenode, epart)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
  integer, intent(in) :: rank, broadcaster
  type(element), intent(in) :: global_elements(:)
  integer, intent(in) :: aggregatenode
  integer, allocatable, intent(out) :: epart(:)

  integer, allocatable :: nodes(:), npart(:), work(:)
  integer :: n, i, j, ierr
  integer, pointer :: null_i => null()
  real(8), pointer :: null_r => null()

  destroy(epart)
  allocate(epart(nels))

  if (rank==broadcaster) then
    allocate(nodes(aggregatenode),work(nels+1), npart(nnds)); work(1)=0
    n=0
    do i=1,nels
      j=n+1; n=n+global_elements(i)%nodecount
      nodes(j:n)=global_elements(i)%nodes
      work(i+1)=n
    end do
    nodes=nodes-1
    
    call METIS_PartMeshNodal(nels,nnds,work,nodes(:),null_i,null_i,nprcs,     &
           null_r,null_i,n,epart,npart)
    deallocate(npart, nodes,work)
  end if
  call MPI_Bcast(epart,nels,MPI_Integer,broadcaster,MPI_Comm_World,ierr)

  end subroutine PartitionBroadcast

  subroutine SendElement(dst, el)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
  integer, intent(in) :: dst
  type(element), intent(in) :: el

  integer :: ierr
  integer, parameter :: BUFFER_SIZE = 10
  integer :: i_buffer(BUFFER_SIZE)
  real(8) :: r_buffer(BUFFER_SIZE, BUFFER_SIZE)
  
  i_buffer(1:3) = (/getElTypeNo(el%eltype),&
                    el%nodecount,&
                    el%mat & 
                  /)
  call MPI_Send(i_buffer, BUFFER_SIZE, MPI_INT, dst, 0, MPI_COMM_WORLD, ierr)
  i_buffer(1:el%nodecount) = el%nodes
  call MPI_Send(i_buffer, BUFFER_SIZE, MPI_INT, dst, 0, MPI_Comm_World, ierr)
  r_buffer(:el%nodecount, :pdim) = el%ecoords
  call MPI_Send(r_buffer, BUFFER_SIZE*BUFFER_SIZE, MPI_REAL8, dst, 0, MPI_COMM_WORLD, ierr)
  
  end subroutine SendElement

  subroutine RecvElement(src, el)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
  integer, intent(in) :: src
  type(element), intent(out) :: el

  integer :: ierr
  integer, parameter :: BUFFER_SIZE = 10
  integer :: i_buffer(BUFFER_SIZE)
  real(8) :: r_buffer(BUFFER_SIZE, BUFFER_SIZE)

  call MPI_Recv(i_buffer, BUFFER_SIZE, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  el%eltype = elTypes(i_buffer(1))
  el%nodecount = i_buffer(2)
  el%mat = i_buffer(3)
  call MPI_Recv(i_buffer, BUFFER_SIZE, MPI_INT, src, 0, MPI_Comm_World, MPI_STATUS_IGNORE, ierr)
  el%nodes = i_buffer(:el%nodecount)
  call MPI_Recv(r_buffer, BUFFER_SIZE*BUFFER_SIZE, MPI_REAL8, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  el%ecoords = r_buffer(:el%nodecount, :pdim)
  
  end subroutine RecvElement

  subroutine DistributeElements(rank, distributor, nels, epart, global_elements, local_elements)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: rank, distributor, nels, epart(:)
    type(element), intent(in) :: global_elements(:)
    type(element), allocatable, intent(out) :: local_elements(:)

    integer :: i, j

    destroy(local_elements)
    j = 0
    do i = 1, nels
      if(epart(i)==rank) j = j + 1
    end do
    allocate(local_elements(j))

    j = 1
    do i = 1, nels

      if(rank==distributor) then
                                if(epart(i)==rank) then
                                  local_elements(j) = global_elements(i)
                                  j = j + 1
                                else
                                  call SendElement(epart(i), global_elements(i))
                                end if
      else
                                if(epart(i)==rank) then
                                  call RecvElement(distributor, local_elements(j))
                                  j = j + 1
                                end if 
      end if
    end do

  end subroutine DistributeElements

  subroutine SetNonLinEls(nlels, local_elements, l_nonlin_ec, nonlinels)
    implicit none
    integer, intent(in) :: nlels
    type(element), intent(in) :: local_elements(:)
    integer, intent(out) :: l_nonlin_ec
    integer, allocatable, intent(out) :: nonlinels(:)

    integer i, j

    l_nonlin_ec = 0
    destroy(nonlinels)

    do i = 1, nlels
      if(local_elements(i)%eltype == "coh") l_nonlin_ec = l_nonlin_ec + 1
    end do
    allocate(nonlinels(l_nonlin_ec))
    j = 1
    do i = 1, nlels
      if(local_elements(i)%eltype == "coh") then
        nonlinels(j) = i
        j = j + 1
      end if
    end do
    
  end subroutine SetNonlinEls

  subroutine SetCoords(nlels, local_elements, coords)
    implicit none
    integer, intent(in) :: nlels
    type(element), intent(in) :: local_elements(:)
    real(8), allocatable, intent(out) :: coords(:, :)

    integer :: i, j

    destroy(coords)
    allocate(coords(nlnds, pdim))
    do i = 1, nels
      do j = 1, local_elements(i)%nodecount
        coords(local_elements(i)%nodes(j), :) = local_elements(i)%ecoords(j, :)
      end do
    end do
  end subroutine SetCoords

  subroutine ReadDistMaterials(file, distributor, nmat, mats)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: nmat, file, distributor
    real(8), allocatable, intent(out) :: mats(:,:)

    integer, parameter :: ELASTIC_MAT_SIZE = 5

    integer :: i, ierr

    destroy(mats)
    allocate(mats(nmat, ELASTIC_MAT_SIZE))

    if(rank == distributor) then
      do i = 1, nmat
        read(file,*), mats(i, :)
      end do
    end if
    call MPI_Bcast(mats, nmat, MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)

  end subroutine ReadDistMaterials

  subroutine ReadDistForces(file, nfrcs, nprcs, rank, distributor, fnode, fval)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: file, nprcs, rank, distributor
    integer :: nfrcs
    integer, allocatable, intent(out) :: fnode(:)
    real(8), allocatable, intent(out) :: fval(:,:)
    
    integer :: ierr
    integer :: i, currentProcessor, localForceCount
    integer :: i_buffer(1), tempFNode(nfrcs)
    real(8) :: r_buffer(pdim + 2), tempFVal(nfrcs, pdim + 2)

    destroy(fnode)
    destroy(fval)
    localForceCount = nfrcs / nprcs
    if (mod(nfrcs, nprcs) > rank) localForceCount = localForceCount + 1
    ! allocate(fnode(localForceCount), fval(localForceCount, pdim + 2))
    allocate(fnode(localForceCount), fval(localForceCount, pdim + 2))

    if(rank == distributor) then
      do i = 1, nfrcs
      read(file,*)tempFNode(i),tempFVal(i,:)
       tempFVal(i, pdim + 1) = min(tempFVal(i, pdim + 1), t)
       tempFVal(i, pdim + 2) = min(tempFVal(i, pdim + 2), t)
      end do
    end if

    call MPI_Bcast(tempFVal, nfrcs * (pdim + 2), MPI_REAL8, distributor, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tempFNode, nfrcs, MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)

    do i = 1, localForceCount
      fnode(i) = tempFNode((i - 1) * nprcs + rank + 1)
      fval(i,:) = tempFVal((i - 1) * nprcs + rank + 1, :)
    end do

    nfrcs = localForceCount
  end subroutine ReadDistForces

  subroutine ReadDistBcs(file, nbcs, nprcs, rank, distributor, bcnode, bcval)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: file, nprcs, rank, distributor
    integer :: nbcs
    integer, allocatable, intent(out) :: bcnode(:, :)
    real(8), allocatable, intent(out) :: bcval(:,:)
    
    integer :: ierr
    integer :: i, currentProcessor, localBcCount
    integer :: i_buffer(1), tempBcNode(nbcs, 1 + pdim)
    real(8) :: r_buffer(pdim + 2), tempBcVal(nbcs, pdim)

    destroy(bcval)
    destroy(bcnode)
    localBcCount = nbcs / nprcs
    if (mod(nbcs, nprcs) > rank) localBcCount = localBcCount + 1
    allocate(bcnode(localBcCount, 1 + pdim), bcval(localBcCount, pdim))

    if(rank == distributor) then
      do i = 1, nbcs
      read(file,*)tempBcNode(i, :),tempBcVal(i,:)
      end do
    end if

    call MPI_Bcast(tempBcVal, nbcs * (pdim), MPI_REAL8, distributor, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tempBcNode, nbcs * (pdim + 1), MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)

    do i = 1, localBcCount
      bcnode(i, :) = tempBcNode((i - 1) * nprcs + rank + 1, :)
      bcval(i,:) = tempBcVal((i - 1) * nprcs + rank + 1, :)
    end do

    nbcs = localBcCount
  end subroutine ReadDistBcs
end module io
