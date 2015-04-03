program main
  implicit none
  ! Output parameters
  character(*), parameter :: solver = "implicit"
  integer :: pdim, guess, nels, nnds, nmats, nceqs, nfrcs
  integer :: ntrcts, nbcs
  real(8) :: total_time, dt
  integer :: obsolete_output_freq, writeType
  integer, allocatable :: elements(:, :)
  real(8), allocatable :: coords(:, :)
  real(8), allocatable :: mats(:, :)
  integer, allocatable :: bc_nodes(:, :)
  real(8), allocatable :: bc_vals(:, :)
  integer, allocatable :: force_nodes(:)
  real(8), allocatable :: forces(:, :)
  
  ! Calculation parameters
  integer :: x_nels, y_nels, x_nnds, y_nnds
  
  ! Iterators
  integer :: i, j
  real(8) :: current_x
  real(8) :: current_y
  
  ! Convenience identifiers
  integer :: bc_0_node, bc_1_node, force_0_node, force_1_node
  
  ! Input buffers
  character(32) :: buffer
  

  
  pdim = 2
  x_nels = 10
  y_nels = 1
  
  if(iargc() == 2) then
    call getarg(1, buffer)
    read(buffer, '(I10)'), x_nels
    call getarg(2, buffer)
    read(buffer, '(I10)'), y_nels 
  end if
  
  x_nnds = x_nels + 1
  y_nnds = y_nels + 1
  guess = x_nels * y_nels
  
  nels = x_nels * y_nels
  nnds = x_nnds * y_nnds
  nmats = 1
  nceqs = 0 ! Not supported yet
  nfrcs = 2
  ntrcts = 0 ! Not supported yet
  nbcs = 2
  
  total_time = 0.01
  dt = 0.01
  obsolete_output_freq = 1
  writeType = 1
  
  !! Set the elements
  allocate(elements(nels, 5))
  elements(:, 5) = 1 ! Set the element material
  do i = 1, y_nels
    do j = 1, x_nels
      elements((i-1)*x_nels+j,1:4)=(/j + (i-1)*x_nnds,     &
                                     j + 1 +(i-1)*x_nnds,  &
                                     j + 1 + i*x_nnds,     &
                                     j + i*x_nnds          &
                                   /)
    end do
  end do
  !! Elements set
  
  !! Set the coordinates
  bc_0_node = 1
  bc_1_node = 1 + (y_nnds - 1) * x_nnds
  allocate(coords(nnds, pdim))  
  do i = 1, y_nnds
    current_y = (i - 1) * 1.0
    do j = 1, x_nnds
      current_x = (j - 1) * 1.0
      coords(j + (i - 1)*x_nnds, :) = (/current_x, current_y/)
    end do
  end do
  !! Coordinates set
  
  !! Set materials
  allocate(mats(nmats, 5))
  do i = 1, nmats
    mats(i, :) = (/3.0E10, 0.25, 1.0E18, 1.0, 3000.0/)
  end do
  !! Materials set
  
  !! Set bcs
  allocate(bc_nodes(2, pdim + 1), bc_vals(2, 2))
  bc_nodes(1, :) = (/bc_0_node, 0, 0/)
  bc_nodes(2, :) = (/bc_1_node, 0, 0/)
  bc_vals = 0.0
  !! Bcs set
  
  !! Set forces
  force_0_node = x_nnds
  force_1_node = x_nnds * y_nnds
  allocate(force_nodes(2), forces(2, 4))
  force_nodes = (/force_0_node, force_1_node/)
  forces(:, 1) = -10.0E10
  forces(:, 2) = 0.0
  forces(:, 3) = 0.0
  forces(:, 4) = total_time
  !! Forces set
  
  open(10, file="examples/generated_example.inp", status='replace')
  write(10, '(A,I0,1X,I0)'), solver // " ", pdim, guess
  write(10, '(7(I0,1X))'), nels, nnds, nmats, nceqs, nfrcs, ntrcts, nbcs
  write(10, '((2(F0.6,1X))(2(I0,1X)))'), total_time, dt, obsolete_output_freq, writeType
  write(10, '()')
  do i = 1, nels
    write(10, '((A)5(I0,1X))'), "qua ", elements(i, :)
  end do
  write(10, '()')
  do i = 1, nnds
    write(10, '(2(F0.6,1X))'), coords(i, :)
  end do
  write(10, '()')
  do i = 1, nmats
    write(10, '(5(F0.6,1X))'), mats(i, :)
  end do
  
  write(10, '()')
  do i = 1, nbcs
    write(10, '(3(I0,1X)2(F0.6,1X))'), bc_nodes(i,:), bc_vals(i, :)
  end do
  
  write(10, '()')
  do i = 1, nfrcs
    write(10, '((I0,1X)4(F0.6,1X))'), force_nodes(i), forces(i, :)
  end do
  
  close(10)
  deallocate(elements, coords, mats, bc_nodes, bc_vals)
  deallocate(force_nodes, forces)
  
end program main
