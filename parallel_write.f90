!
! Use MPI-IO to write a diagonal matrix distributed block-cyclically,
!
program darraympiio
  use mpi
  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer :: n, nb ! problem size and block size
  real(dp), dimension(:), allocatable :: buff ! buffer
  real(dp) :: rTmp

  integer, external :: numroc ! blacs routine
  integer :: procs, prow, pcol ! blacs data

  character(12) :: filename = "testdat.bin"
  integer :: mpirank
  integer :: ierr
  integer, dimension(2) :: pdims, dims, distribs, dargs
  integer :: fh
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus
  integer :: newtype
  integer :: locsize, nelements
  integer(kind=MPI_ADDRESS_KIND) :: lb, locextent
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer :: ii

  n = 27

  ! Initialize MPI (for MPI-IO)
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, procs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierr)


  ! May as well get the process grid from MPI_Dims_create
  pdims = 0
  call MPI_Dims_create(procs, 2, pdims, ierr)
  prow = pdims(1)
  pcol = pdims(2)

  ! read size of array from first file entry
  call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL,&
      & fh, ierr)
  call MPI_File_write_all(fh, n, 1, MPI_INTEGER, mpistatus, ierr)
  call MPI_File_close(fh, ierr)

  ! create the darray that will be passed to write the data.
  nb = 4 ! small
  if (nb > (N/prow)) nb = N/prow
  if (nb > (N/pcol)) nb = N/pcol
  dims = [n,n]
  distribs(:) = MPI_DISTRIBUTE_CYCLIC
  dargs = [nb, nb]
  call MPI_Type_create_darray(procs, mpirank, 2, dims, distribs, dargs, pdims, MPI_ORDER_FORTRAN,&
      & MPI_DOUBLE_PRECISION, newtype, ierr)
  call MPI_Type_commit(newtype,ierr)

  call MPI_Type_size(newtype, locsize, ierr)
  call MPI_Type_get_extent(newtype, lb, locextent, ierr)

  ! Initialize local arrays
  nelements = locsize / sizeof(rTmp)
  allocate(buff(nelements))
  do ii = 1, nElements
    buff(ii) = ii + mpirank * 10000
  end do

  ! write the data
  call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_APPEND + MPI_MODE_WRONLY, MPI_INFO_NULL,&
      & fh, ierr)
  disp = 4 ! skip 4 bytes for first integer
  call MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, newtype, "native", MPI_INFO_NULL, ierr)
  call MPI_File_write_all(fh, buff, nelements, MPI_DOUBLE_PRECISION, mpistatus, ierr)
  call MPI_File_close(fh,ierr)

  ! Deallocate the local arrays
  deallocate(buff)

  ! End blacs for processors that are used
  call MPI_Finalize(ierr)

end program darraympiio
