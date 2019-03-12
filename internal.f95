subroutine inter
use vars
integer nx,ny,c1
real*8::linmax(2),linmin(2)
real*8::lx,ly
nx=30
ny=30
ni=nx*ny
allocate(cordi(ni,2))


linmax=0
linmin=0

linmax(1)=maxval(coord(:,1))
linmin(1)=minval(coord(:,1))
linmax(2)=maxval(coord(:,2))
linmin(2)=minval(coord(:,2))



lx=linmax(1)-linmin(1)
ly=linmax(2)-linmin(2)
lx=lx/(nx+1)
ly=ly/(ny+1)


c1=0
do i=1,nx
  do j=1,ny
    c1=c1+1
    cordi(c1,1)=linmin(1)+lx*i     
    cordi(c1,2)=linmin(2)+ly*j
  enddo
enddo

do i=1,ni
open(7,file='Internal nodes.txt',status='unknown')
	write(7,*)i,cordi(i,1),cordi(i,2)
enddo
 close(7)
 
end