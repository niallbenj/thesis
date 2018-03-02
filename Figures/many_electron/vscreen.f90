program vscreen

integer :: i,j,k
real :: r(40001), V1(40000), V2(40000), step
!
step = 0.02
r = 0.0
!
open(10,file='toplot')
do i=1,40000
r(i+1) = r(i) + step
V1(i) = -10.0/r(i)
V2(i) = -1.0/r(i)
write(10,*) r(i+1), V1(i), V2(i)
enddo


endprogram
