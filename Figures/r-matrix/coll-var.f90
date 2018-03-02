program coll_var

integer :: i,j,k
real :: E(40001), COLL1(40000), COLL2(40000), COLL3(40000), &
      & step, r
!
step = 0.02
r = 0.0
!
open(10,file='toplot')

do i=1,40000
  r = r + step
COLL1(i) = 5.0*log(r)
COLL2(i) = 5.0
COLL3(i) = 5.0/(r**(2.0))
write(10,*) r, COLL1(i), COLL2(i), COLL3(i)
enddo


endprogram
