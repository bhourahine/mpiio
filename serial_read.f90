program read
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer n, fh, ii
  real(dp), allocatable :: matrix(:,:)

  OPEN(newunit=fh, file="testdat.bin", form='unformatted', access="stream")
  read(fh) n
  write(*,*) n

  allocate(matrix(n,n))
  matrix(:,:) = 0.0
  read(fh)matrix

  do ii = 1, n
    write(*,"(40I6)")int(matrix(:, ii))
  end do

  close(fh)

end program read
