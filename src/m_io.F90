module io
  use elems
  implicit none 
  
  contains
  
  subroutine broadcastElement(rank, broadcaster, el)
    implicit none
#if defined ALP_PC
#include <finclude/petsc.h90>
#else
#include <petsc-finclude/petsc.h90>
#endif
#define destroy(x) if(allocated(x)) deallocate(x)
    ! Input
    integer, intent(in) :: rank, broadcaster
    type(element) :: el
    
    ! Constants
    integer, parameter :: CONST_BUFFER_COUNT = 3
    
    ! Variables
    integer, allocatable :: int_buffer(:)
    integer :: const_int_buffer(CONST_BUFFER_COUNT)
    real(8), allocatable :: real_buffer(:, :)
    integer :: pdim, ierr
    
    if(rank==broadcaster) then
      const_int_buffer = (/getElTypeNo(el%eltype), el%nodecount, el%mat/)  
    end if
    
    call MPI_Bcast(const_int_buffer, CONST_BUFFER_COUNT, &
                   MPI_INTEGER, broadcaster, MPI_COMM_WORLD, ierr) 
    
    if(rank/=broadcaster) then
      el%eltype = getElName(const_int_buffer(1))
      el%nodecount = const_int_buffer(2)
      el%mat = const_int_buffer(3)
      destroy(el%nodes)
      destroy(el%ecoords)
    end if
    pdim = getDim(el%eltype)
    if(rank /= broadcaster) then
      allocate(el%nodes(el%nodecount))
      allocate(el%ecoords(el%nodecount, pdim))
    end if
    

    
    allocate(int_buffer(el%nodecount))
    allocate(real_buffer(el%nodecount, pdim))
    if(rank==broadcaster) then
      int_buffer = el%nodes
      real_buffer = el%ecoords
    end if
    
    call MPI_Bcast(int_buffer, el%nodecount, MPI_INTEGER, &
                   broadcaster, MPI_COMM_WORLD, ierr)
    
  end subroutine broadcastElement
  
end module io
