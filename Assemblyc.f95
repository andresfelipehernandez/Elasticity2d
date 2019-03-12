subroutine assembly

use vars

allocate (ma(2*nn,2*nn))
allocate (b(2*nn))
allocate (bcc(2*nn))
allocate (bc(2*nn))


open(20,file='BoundCond.txt',status='unknown')


do j=1,nn
  if(coord(j,1)==minval(coord(:,1)))then
    bcc(2*j-1)=1
    bc(2*j-1)=0
    bcc(2*j)=1
    bc(2*j)=0
    write(20,*)'Node number',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Left'
  else if (coord(j,1)==maxval(coord(:,1))) then
    bcc(2*j-1)=2
    bc(2*j-1)=0
    bcc(2*j)=2
    bc(2*j)=0
    write(20,*)'Node number',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Right'
  else if ((coord(j,2)==minval(coord(:,2))).and.(coord(j,1)>minval(coord(:,1))).and.(coord(j,1)<maxval(coord(:,1)))) then
   bcc(2*j-1)=2
    bc(2*j-1)=0
    bcc(2*j)=2
    bc(2*j)=0
    write(20,*)'Node number',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Down'
  else if ((coord(j,2)==maxval(coord(:,2))).and.(coord(j,1)>minval(coord(:,1))).and.(coord(j,1)<maxval(coord(:,1)))) then
    bcc(2*j-1)=2
    bc(2*j-1)=0
    bcc(2*j)=2
    bc(2*j)=-0.01
    write(20,*)'Node number',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Up'
  endif
  
enddo

close(20)
open(15,file='matrizA.txt',status='unknown')

!Create matrix A, called ma, and vector b. To create Ax=b system

b=0.0d0
ma=0.0d0

!$$$$$$ 
 do i=1,2*nn
 if (bcc(i).eq.1) then
    ma(:,i)=-G(:,i)
    b=b-H(:,i)*bc(i)
  elseif (bcc(i).eq.2) then
    ma(:,i)=H(:,i)
    b=b+G(:,i)*bc(i)
   endif
enddo

 do i=1,2*nn
  open(5,file='B.txt',status='unknown')
  write(5,*)i,b(i)
  enddo
close(5)

return
end

    
subroutine disassembly

use vars

allocate(u_b(2*nn)); allocate(t_b(2*nn)); allocate(R_t(2*nn))

do i=1,2*nn
  if (bcc(i).eq.1) then 
    u_b(i)=bc(i)
    t_b(i)=b(i)
  else if (bcc(i).eq.2) then 
    t_b(i)=bc(i)
	u_b(i)=b(i)
  endif
enddo

return
end

