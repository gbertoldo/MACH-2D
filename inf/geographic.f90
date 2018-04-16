program main
  implicit none
  integer :: i, j, nx, ny
  integer :: np, npsw, nps, npse, npw, npe, npnw, npn, npne
  nx = 4
  ny = 3
  do j = 1, ny
     do i = 1, nx
        np   = nx * (j-1) + i
        nps  = np - nx
        npn  = np + nx
        npw  = np - 1
        npe  = np + 1
        npsw = nps - 1
        npse = nps + 1
        npnw = npn - 1
        npne = npn + 1
        write(*,*) npsw, nps, npse, npw, np, npe, npnw, npn, npne
     end do
  end do
end program main
