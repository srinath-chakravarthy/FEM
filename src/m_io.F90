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
    read(file,*)nels,nnds,nmts,ncohmats,nceqs,nfrcs,ntrcs,nbcs
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
    integer, parameter :: INT_BUFFER_SIZE = 11, REAL_BUFFER_SIZE = 2

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
      int_buffer(2:) = (/pdim, nodal_bw, nels, nnds, nmts, nceqs, nfrcs, ntrcs, nbcs, ncohmats/)
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
        ncohmats = int_buffer(11)
        t = real_buffer(1)
        dt = real_buffer(2)
    end if 
  end subroutine BroadCastParameters

  subroutine ReadElementsCoords(file, global_elements, aggregatenode)
    implicit none
    integer, intent(in) :: file
    type(element), allocatable, intent(out) :: global_elements(:)
    integer, intent(out) :: aggregatenode

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
      backspace(file)
      global_elements(i)%nodecount = getNodeCount(global_elements(i)%eltype)
      allocate(global_elements(i)%nodes(global_elements(i)%nodecount))
      allocate(global_elements(i)%ecoords(global_elements(i)%nodecount, pdim))
      aggregatenode= aggregatenode+(global_elements(i)%nodecount)
      read(file, *)dummy, global_elements(i)%nodes, global_elements(i)%mat, global_elements(i)%nlMat
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
  
  i_buffer(1:4) = (/getElTypeNo(el%eltype),&
                    el%nodecount,&
                    el%mat, &
                    el%nlMat & 
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
  el%nlMat = i_buffer(4)
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

  subroutine ReadDistMaterials(file, distributor, nmat, mats, ncohmats, cohmats)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: nmat, ncohmats, file, distributor
    real(8), allocatable, intent(out) :: mats(:,:)
    type(cohMat), allocatable, intent(out) :: cohmats(:)

    integer :: i, ierr, temp

    destroy(mats)
    destroy(cohmats)
    allocate(mats(nmat, ELASTIC_MAT_SIZE))
    allocate(cohmats(ncohmats))

    if(rank == distributor) then
      do i = 1, nmat
        read(file,*), mats(i, :)
      end do
    end if
    call MPI_Bcast(mats, nmat, MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)

    if(rank == distributor) then
      do i = 1, ncohmats
        read(file,*) cohmats(i)%seplaw
        backspace(file)
        cohmats(i)%propCount = PROP_COUNTS(cohmats(i)%seplaw)
        allocate(cohmats(i)%props(cohmats(i)%propCount))
        read(file,*) temp, cohmats(i)%props
      end do
    end if

    do i = 1, ncohmats
      if (rank == distributor) temp = cohmats(i)%seplaw
      call MPI_Bcast(temp, 1, MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)
      if (rank /= distributor) then
        cohmats(i)%seplaw = temp
        cohmats(i)%propCount = PROP_COUNTS(cohmats(i)%seplaw)
        allocate(cohmats(i)%props(cohmats(i)%propCount))
      end if
      call MPI_Bcast(cohmats(i)%props, cohmats(i)%propCount, MPI_REAL8, distributor, MPI_COMM_WORLD, ierr)
    end do 

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

  ! epart must be set
  subroutine ReadDistTractions(file, ntrcs, nprcs, rank, distributor, telsd, tval)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    integer, intent(in) :: file, nprcs, rank, distributor
    integer :: ntrcs
    integer, allocatable, intent(out) :: telsd(:, :)
    real(8), allocatable, intent(out) :: tval(:,:)
    
    integer :: ierr, i, i2, j
    integer :: localTractCount
    integer :: i_buffer(2), tempTelsd(ntrcs, 2)
    real(8) :: r_buffer(pdim + 2), tempTval(ntrcs, pdim + 2)
    integer, allocatable :: emap(:)

    destroy(telsd)
    destroy(tval)
    localTractCount = 0

    if(rank == distributor) then
      do i = 1, ntrcs
      read(file,*)tempTelsd(i, :),tempTval(i,:)
      end do
    end if

    call MPI_Bcast(tempTelsd, ntrcs * 2, MPI_INTEGER, distributor, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tempTval, ntrcs * (pdim + 2), MPI_REAL8, distributor, MPI_COMM_WORLD, ierr)

    ! Calculate local traction count
    do i = 1, size(tempTelsd)
      if(epart(tempTelsd(i, 1)) == rank) localTractCount = localTractCount + 1
    end do

    allocate(telsd(localTractCount, 2), tval(localTractCount, pdim + 2))
    
    j = 1
    do i = 1, size(tempTelsd)
      if(epart(tempTelsd(i, 1)) == rank) then
        telsd(j, :) = tempTelsd(i, :)
        tval(j, :) = tempTval(i, :)
        j = j + 1
      end if  
    end do

    ! Localize element numbers
    allocate(emap(size(epart)))
    emap = 0
    j = 1
    do i = 1, size(epart)
      if(epart(i) == rank) then
        emap(i) = j
        j = j + 1
      end if
    end do
    telsd(:, 1) = emap(telsd(:, 1))
    deallocate(emap)

    ntrcs = localTractCount
  end subroutine ReadDistTractions


  
  ! Write results in ASCII VTK (legacy) format
  subroutine WriteOutput
    implicit none
#if defined ALP_PC
#include <finclude/petscmat.h90>
#include <finclude/petscsys.h>
#else
#include <petsc-finclude/petscmat.h90>
#include <petsc-finclude/petscsys.h>
#endif
    character(64) :: name,fmt
    character(32) :: buffer
    integer,save :: k=0
    integer :: i,j,j1,lnnds,lnels
    real(8),pointer :: field_val(:), stress_val(:, :)
    field_val => aggregate_u
    stress_val => aggregate_stress
    write(name,'(I0,A,I0.6,A)')rank,"_output_",k,".vtk"
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(local_elements)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Defmod"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    select case(pdim)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:))/)
       end do
    end select
    ! j stores the total number of columns
    j = 0
    do i=1, lnels
      j = j + local_elements(i)%nodecount+1
    end do
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels, j
    do i=1,lnels
       write(buffer, '(I0)') local_elements(i)%nodecount
       fmt = "(I0," // trim(buffer) // "(1X,I0))"
       write(10,fmt)local_elements(i)%nodecount,local_elements(i)%nodes-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)') getVtkid(local_elements(i)%eltype)
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    write(10,'(A,I0)') "SCALARS STRESS FLOAT ", cpdim
    write(10,'(A)') "LOOKUP_TABLE DEFAULT"
    write(buffer, '(I0)') cpdim
    fmt = "(" // trim(buffer) // "(F0.6,1X))"
    do i = 1, lnnds
      write(10, fmt) stress_val(i, :)
    end do
    j=pdim
    write(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    select case(pdim)
    case(2)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       end do
    end select
    close(10); k=k+1
  end subroutine WriteOutput
end module io
