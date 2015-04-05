! --------------- Srinath Chakravarthy --------------------
! --------- 2D simple Elastic finite element solver -------
! ------ copied from defmod google code ---------
! -----------------------------------------

#define ALP_PC

program main

  use global
  use io
  implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
#define destroy(x) if(allocated(x)) deallocate(x)
  character(256) :: input_file, dummy, buffer
  integer :: i,j,j1,j2,n, ef_eldof, aggregatenode
  integer,pointer :: null_i=>null(), el_nodes(:)
  real(8),pointer :: null_r=>null()

  type(element) :: dummy_el
  
  call PetscInitialize(Petsc_Null_Character,ierr)

  call MPI_Comm_Rank(MPI_Comm_World,rank,ierr)
  call MPI_Comm_Size(MPI_Comm_World,nprcs,ierr)

  call PetscOptionsGetString(Petsc_Null_Character,'-f',input_file,n,ierr)
  if (n==0) call PrintMsg("Usage: mpiexec -n <cores> defmod -f <filename>")
  if (n==0) go to 9

  ! Read input file parameters
  if(rank == 0) open(10,file=input_file,status='old')
  call PrintMsg("Reading input ...")
  
  if(rank == 0) call ReadParameters(10) !; if (stype /= "implicit") go to 9
  call BroadCastParameters(rank, 0)
  call Initialize(ipoints, weights)
  if(rank == 0) call ReadElementsCoords(10, global_elements, aggregatenode, nonlinear)
  call PartitionBroadcast(rank, 0, global_elements, aggregatenode, epart)
  call DistributeElements(rank, 0, nels, epart, global_elements, local_elements, nonlinear)
  if(rank == 0) deallocate(global_elements)
  nels=size(local_elements)
  call SetNonlinEls(nels, local_elements, l_nonlin_ec, nonlin_els)

  allocate(npart(nnds))
  allocate(work(nnds))
  do i=1,nels
     do j=1,local_elements(i)%nodecount
        npart(local_elements(i)%nodes(j)) = 1
     end do
  end do
  
  
  ! work is global to local node numbers
  !     work(i) gives global node i's local number
  ! j1 is the nonlinear element iteration
  j=1; work=0
  do i=1,nnds
    if (npart(i)==1) then
       work(i)=j; j=j+1
    end if
  end do
  nlnds = j - 1 ! nlnds set
  do i=1,nels
    do j=1,local_elements(i)%nodecount
       ! local_elements' nodes contains local node numbers
       local_elements(i)%nodes(j)=work(local_elements(i)%nodes(j))
    end do
  end do
  deallocate(npart)

  call SetCoords(nels, local_elements, coords)
  
  allocate(nl2g(nlnds),indxmap((pdim)*nlnds,2))
  j=1
  do i=1,nnds
     if (work(i)==j) then
        nl2g(j)=i; j=j+1
     end if
  end do
  do i=1,nlnds
     do j=1,pdim
        ! indxmap(i, j):
        !     local ith dof's global number.
        !     A dof's no is attained from its node's number.
        ! example: dof per el = 3
        !     node #1: 0, 1, 2
        !     node #3: 6, 7, 8
        indxmap((pdim)*i-j+1,1)=(pdim)*i-j
        indxmap((pdim)*i-j+1,2)=(pdim)*nl2g(i)-j
     end do
  end do
  deallocate(work)
  
  call ReadDistMaterials(10, 0, nmts, mat)
  
  call ReadDistBcs(10, nbcs, nprcs, rank, 0, bcnode, bcval)
  call ReadDistForces(10, nfrcs, nprcs, rank, 0, fnode, fval)
  deallocate(epart)
  if(rank == 0) close(10)

! TODO: nceqs
!  allocate(cval(nceqs,3)); cval=f0
! TODO: tractions
! allocate(telsd(ntrcs,2),tval(ntrcs,pdim+2)); tval=f0
!  allocate(emap(nels)); emap=0
!  j=1
!  do i=1,nels
!     if (epart(i)==rank) then
!        epart(i)=1; emap(i)=j; j=j+1
!     else
!        epart(i)=0
!     end if
!  end do

!  do i=1,ntrcs
!     read(10,*)telsd(i,:),tval(i,:)
!  end do
!  telsd(:,1)=emap(telsd(:,1)) ! Remap nodes/els
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT END !!!!!!!!!!!!!!!!!!!

  ! Initialize aggregate result values
  allocate(aggregate_u(nlnds * pdim), aggregate_stress(nlnds, cpdim))
  aggregate_u = 0
  aggregate_stress = 0

  ! Initialize local element variables and global U
  allocate(vvec(pdim))  
  ! n: total dof cont
  n=(pdim)*nnds  
  call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
  allocate(stress(nels))
  do i = 1, nels
    allocate(stress(i)%array(getNip(local_elements(i)%eltype), cpdim))
  end do
  ! Set cpdim index set, to be used later
  allocate(work(cpdim))
  do i = 1, cpdim
    work(i) = i - 1
  end do
  call ISCreateGeneral(Petsc_Comm_Self, cpdim, work, Petsc_Copy_Values, Cpdim_Iteration, ierr)
  deallocate(work)
  
  ! Form non-cohesive stiffness matrix
  call PrintMsg("Forming [K] ...")
  nodal_bw=pdim*(nodal_bw+1)
  call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,   &
     Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)
  call MatSetOption(Mat_K,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
  do i=1,nels
     if(local_elements(i)%eltype == "coh") cycle
     allocate(k(local_elements(i)%nodecount*pdim,local_elements(i)%nodecount*pdim))
     allocate(indx(local_elements(i)%nodecount*pdim))
     call FormLocalK(i,k,indx)
     indx=indxmap(indx,2) 
     call MatSetValues(Mat_K,local_elements(i)%nodeCount*pdim,indx,  &
                      local_elements(i)%nodeCount*pdim, indx,k,Add_Values,ierr)
     deallocate(k)
     deallocate(indx)
  end do
  call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
  call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
  call ApplyKBC(Mat_K)
  call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
  call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
!!$  call PrintMsg("Elastic Stiffness Only")
!!$  call MatView(Mat_K, PETSC_VIEWER_STDOUT_WORLD,ierr)
  allocate(Array_K(nlnds*pdim, nlnds*pdim))
  call GetStiffnessArray(Mat_K, Array_K)


  ! Initialize arrays to communicate ghost node values
  call PrintMsg("Setting up environment ...")
  ! Create nodal stress matrix
  call MatCreateDense(Petsc_Comm_World, Petsc_Decide, Petsc_Decide, nnds, cpdim, Petsc_Null_Scalar, Mat_Stress, ierr)
  call MatSetOption(Mat_Stress, Mat_Row_Oriented, Petsc_False, ierr)
  call VecDuplicate(Vec_U,Vec_F,ierr)
  j=pdim*nlnds
  call VecCreateSeq(Petsc_Comm_Self,j,Seq_U,ierr)
  call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,2),Petsc_Copy_Values,From,  &
     ierr)
  call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,1),Petsc_Copy_Values,To,    &
     ierr)
  call VecScatterCreate(Vec_U,From,Seq_U,To,Scatter,ierr)
  allocate(uu(j))  ! Array for Seq_U extraction 
  allocate(count_node(nlnds)) ! Stress contribution count
  allocate(global_count_node(nnds)) ! Global stress contribution count
  allocate(stress_node(nlnds, cpdim)) ! Local stress array
  call MatDuplicate(Mat_K, MAT_COPY_VALUES, Jacobian, ierr)
  call VecDuplicate(Vec_U, Residual, ierr)

  ! Implicit Solver
  if (stype/="explicit") then
     call SNESCreate(PETSC_COMM_WORLD, Solver, ierr)
     call SNESGetKSP(Solver, Krylov, ierr)
     call SetupKSPSolver(Krylov)
     call SNESSetType(Solver, SNESNEWTONLS, ierr)
     call SNESSetFunction(Solver, Residual, CalcResidual, PETSC_NULL_OBJECT, ierr)
     call SNESSetJacobian(Solver, Jacobian, Jacobian, CalcJacobian, PETSC_NULL_OBJECT, ierr)
     call SNESSetFromOptions(Solver, ierr)
     
     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#endif
     call SetupKSPSolver(Krylov)
  end if  
  
  dtNo = 1
  do
    t_init = dt * (dtNo - 1)
    if(t_init >= t) exit
    
    write(buffer, '(I0)'), dtNo 
    call PrintMsg("Interval: " // trim(buffer))
    
    call PrintMsg("    Resetting Environment ...")
    uu = 0
    global_count_node = 0
    count_node = 0  
    stress_node = 0
    call VecZeroEntries(Seq_U, ierr)
    ! call VecZeroEntries(Vec_U, ierr)
    call VecZeroEntries(Vec_F, ierr)
    call MatZeroEntries(Mat_Stress, ierr)
    
    ! Form RHS
    call PrintMsg("    Forming RHS ...")
    call FormRHS(t_init, dt)
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
    ! call VecView(Vec_F,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    if (stype/="explicit") then
      call PrintMsg("    Solving ...")
      if(nonlinear) then
        call CalcJacobian(PETSC_NULL_OBJECT, Vec_U, Jacobian, Jacobian, ierr=ierr)
        call SNESSolve(Solver, PETSC_NULL_OBJECT, Vec_U, ierr)
        call SNESGetIterationNumber(Solver, iterationCount, ierr)
        write(buffer, '(I0)'), iterationCount
        call PrintMsg("    SNES Iteration Count: " // trim(buffer))
      else
        call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
      end if
      call VecView(Vec_U,PETSC_VIEWER_STDOUT_WORLD,ierr)
      call GetVec_U(Vec_U, uu); aggregate_u=aggregate_u+uu
     
      call PrintMsg("    Getting stress ....")
      do i = 1, nels
         call RecoverStress(i)
      end do
      !!! Nodal Stress Calculations !!!
      ! Recover the nodal stress
      do i = 1, nels
        allocate(stress_at_el(local_elements(i)%nodecount, cpdim))
        call RecoverNodalStress(i, stress_at_el)
        el_nodes => local_elements(i)%nodes
        do j = 1, local_elements(i)%nodecount
          count_node(el_nodes(j)) = count_node(el_nodes(j)) + 1
          stress_node(el_nodes(j), :) = stress_node(el_nodes(j), :) + &
                                      stress_at_el(j, :)
        end do
        deallocate(stress_at_el)
      end do
      ! Work utilized for forming local contribution array to global_count_node
      allocate(work(nnds)); work = 0; global_count_node = 0
      do i = 1, nlnds
        work(nl2g(i)) = count_node(i)
      end do
      call MPI_AllReduce(work, global_count_node, nnds, MPI_Integer, MPI_Sum, MPI_Comm_World, ierr)
      deallocate(work)
      ! Adjust local stress contribution according to the global contribution
      do i = 1, nlnds
        stress_node(i, :) = stress_node(i, :) / global_count_node(nl2g(i))
      end do
      ! Work utilized for column-indexing of Mat_Stress
      allocate(work(cpdim))
      do i = 1, cpdim
        work(i) = i - 1
      end do
      ! Set the stress values
      call MatSetValues(Mat_Stress, nlnds, nl2g - 1, cpdim, work, stress_node, Add_Values, ierr)
      deallocate(work)
      call MatAssemblyBegin(Mat_Stress, Mat_Final_Assembly, ierr)
      call MatAssemblyEnd(Mat_Stress, Mat_Final_Assembly, ierr)
      ! call MatView(Mat_Stress, Petsc_Viewer_Stdout_World, ierr)
      call GetMat_Stress(stress_node); aggregate_stress = aggregate_stress + stress_node
      !!! Nodal Stress Calculations, End !!!
    end if
      
    dtNo = dtNo + 1
  end do
  
  ! Write output     
  call WriteOutput
  if(stype /= "explicit")  call KSPDestroy(Krylov,ierr)
    
  ! Cleanup
  call PrintMsg("Cleaning up ...")
  call VecScatterDestroy(Scatter,ierr)
  call ISDestroy(To,ierr)
  call ISDestroy(From,ierr)
  
  call VecDestroy(Seq_U,ierr)
  call VecDestroy(Vec_F,ierr)
  call VecDestroy(Vec_U,ierr)
  call MatDestroy(Mat_K,ierr)
  call MatDestroy(Mat_Stress,ierr)
  ! call SNESDestroy(Solver, ierr)
  
  deallocate(stress)
  deallocate(local_elements, count_node, stress_node)
  deallocate(aggregate_u, aggregate_stress)
  deallocate(coords,mat, vvec,indxmap,uu,fnode,fval,nl2g)
  destroy(cval)
  destroy(tval)
  destroy(telsd)
  deallocate(nonlin_els)
  do i = 1, elTypeCount
    destroy(ipoints(i)%array)
    destroy(weights(i)%array)
    destroy(shapeFuncMem(i)%array)
    destroy(nodalStressInvMem(i)%array)
  end do
  destroy(Array_K)
  
  call PrintMsg("Finished")

9 call PetscFinalize(ierr)

contains

  ! Setup solver
  subroutine SetupKSPSolver(Krylov)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
    KSP :: Krylov
    
       call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
      tol = 1.e-9
      call KSPSetTolerances(Krylov,tol,PETSC_DEFAULT_REAL,                       &
     &     PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)

!  Test MUMPS
#if defined(PETSC_HAVE_MUMPS)
      call KSPSetType(Krylov,KSPPREONLY,ierr)
      call KSPGetPC(Krylov,PreCon,ierr)
      call PCSetType(PreCon,PCLU,ierr)
      call PCFactorSetMatSolverPackage(PreCon,MATSOLVERMUMPS,ierr)
      call PCFactorSetUpMatSolverPackage(PreCon,ierr)
      call PCFactorGetMatrix(PreCon,PC_MAT,ierr)

!     sequential ordering
      icntl = 7 
      ival  = 2
      call MatMumpsSetIcntl(PC_MAT,icntl,ival,ierr)

!     threshhold for row pivot detection
      call MatMumpsSetIcntl(PC_MAT,24,1,ierr)
      icntl = 3
      valmum = 1.e-6
      call MatMumpsSetCntl(PC_MAT,icntl,valmum,ierr)

!     compute determinant of A
      call MatMumpsSetIcntl(PC_MAT,33,1,ierr)
#endif

      call KSPSetFromOptions(Krylov,ierr)
      call KSPSetUp(Krylov,ierr)
#if defined(PETSC_HAVE_MUMPS)
      icntl = 3;
      call MatMumpsGetCntl(PC_MAT,icntl,cntl,ierr)
      call MatMumpsGetInfog(PC_MAT,34,infog34,ierr)
      call MatMumpsGetRinfog(PC_MAT,12,rinfo12,ierr)
      call MatMumpsGetRinfog(PC_MAT,13,rinfo13,ierr)
      if (rank .eq. 0) then
         write(6,98) cntl
         write(6,99) rinfo12,rinfo13,infog34
      endif
 98   format('Mumps row pivot threshhold = ',1pe11.2)
 99   format('Mumps determinant=(',1pe11.2,1pe11.2,')*2^',i3)
#endif
!!$    call KSPSetType(Krylov,"gmres",ierr)
!!$    call KSPGetPC(Krylov,PreCon,ierr)
!!$    call PCSetType(PreCon,"asm",ierr)
!!$#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
!!$    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Real,                    &
!!$       Petsc_Default_Real,Petsc_Default_Integer,ierr)
!!$#else
!!$    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Double_Precision,        &
!!$       Petsc_Default_Double_Precision,Petsc_Default_Integer,ierr)
!!$#endif
!!$    call KSPSetFromOptions(Krylov,ierr)
  end subroutine SetupKSPSolver
  
  ! Extract Mat_Stress
  ! Note: IS Cpdim_Iteration must be set
  subroutine GetMat_Stress(local_stress)
    implicit none
    real(8) :: local_stress(nlnds, cpdim)
    real(8), pointer :: pntr(:, :)
    integer, allocatable :: work(:)
    integer :: i
    Mat :: Extraction(1)
    IS :: IS_row
    
    call ISGetSize(Cpdim_Iteration, i, ierr)
    if (i /= cpdim) then
      call PrintMsg("Cpdim_Iteration is not set correctly")
      return
    end if
    
    allocate(work(nlnds))
    do i = 1, nlnds
      work(i) = nl2g(i) - 1
    end do
    call ISCreateGeneral(Petsc_Comm_Self, nlnds, work, Petsc_Copy_Values, IS_row, ierr)
    deallocate(work)
    call MatGetSubMatrices(Mat_Stress, 1, IS_row, Cpdim_Iteration, Mat_Initial_Matrix, Extraction, ierr)
    call MatDenseGetArrayF90(Extraction(1), pntr, ierr)
    local_stress = pntr
    call MatDenseRestoreArrayF90(Extraction(1), pntr, ierr)
    call ISDestroy(IS_row, ierr)
    call MatDestroy(Extraction(1), ierr)
  end subroutine GetMat_Stress

  ! Get the stiffness matrix of local nodes
  subroutine GetStiffnessArray(src, Array_K)
    implicit none
    Mat, intent(in) :: src
    real(8), intent(out) :: Array_K(nlnds*pdim, nlnds*pdim)
    
    Mat :: Extraction(1) 
    IS :: Indices
    integer :: i, j
    real(8), pointer :: pntr(:, :)
    
    call ISCreateGeneral(Petsc_Comm_Self, nlnds*pdim, indxmap(:,2), PETSC_COPY_VALUES, Indices, ierr)
    call MatGetSubMatrices(src, 1, Indices, Indices, MAT_INITIAL_MATRIX, Extraction, ierr)
    call MatSeqAIJGetArrayF90(Extraction(1), pntr, ierr)
    
    Array_K = pntr

    call MatSeqAIJRestoreArrayF90(Extraction(1), pntr, ierr)
    call MatDestroy(Extraction(1), ierr)
    call ISDestroy(Indices, ierr)
  end subroutine GetStiffnessArray

end program main
