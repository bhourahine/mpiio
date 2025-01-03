program write
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer n, fh, ii, jj, kk
  real(dp), allocatable :: matrix(:,:)

  n = 27

  OPEN(newunit=fh, file="testdat.bin", form='unformatted', access="stream")
  write(fh)n

  allocate(matrix(n,n), source=0.0_dp)
  kk = 0
  do ii = 1, n
    do jj = 1, n
      kk = kk + 1
      matrix(jj, ii) = kk
    end do
  end do

  write(fh)matrix

  close(fh)

end program write
