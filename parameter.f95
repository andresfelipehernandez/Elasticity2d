subroutine parameters
use vars

E=2e11                                  ! Modulus E, Pa
nu=0.32                                 ! Poisson ratio, dimensionless
lambda=nu*E/((1+nu)*(1-2*nu))           !Lame's first constant Pa
mu=E/(2*(1+nu))							!Lame's second constant rigid modulus Pa
Sy=2.568                                !Yield strength, Pa

end


subroutine boundarycalc
use vars
allocate(ux_b(nn),uy_b(nn),tx_b(nn),ty_b(nn))
! boundary displacements

R_t=matmul(H,u_b)-matmul(G,t_b)


 do i=1,nn                                                           
    ux_b(i)=u_b(2*(i-1)+1)
    uy_b(i)=u_b(2*(i-1)+2)
    tx_b(i)=t_b(2*(i-1)+1)
    ty_b(i)=t_b(2*(i-1)+2)
enddo

  do i=1,nn
  open(2,file='boundary displacement.txt',status='unknown')
  write(2,*)i,ux_b(i),uy_b(i)
  enddo
close(2)  

end


subroutine internalcalc
use vars
allocate(ux_i(ni),uy_i(ni))
!! internal displacements
 do i=1,ni                                                           
    ux_i(i)=u_i(2*(i-1)+1)
    uy_i(i)=u_i(2*(i-1)+2)
enddo

 do i=1,ni
  open(3,file='Internal displacement.txt',status='unknown')
  write(3,*)i,ux_i(i),uy_i(i)
  enddo
close(3)  

end


subroutine totu
use vars
allocate(ux_t(nn+ni),uy_t(nn+ni),despl(nn+ni,2),coordt(nn+ni,2))

do i=1,nn
  ux_t(i)=ux_b(i)
  uy_t(i)=uy_b(i)
  despl(i,1)=ux_b(i)
  despl(i,2)=uy_b(i)
  coordt(i,1)=coord(i,1)
  coordt(i,2)=coord(i,2)
enddo
do j=1,ni
  ux_t(j+nn)=ux_i(i)
  uy_t(j+nn)=uy_i(i)
  despl(j+nn,1)=ux_i(i)
  despl(j+nn,2)=uy_i(i)
  coordt(j+nn,1)=cordi(i,1)
  coordt(j+nn,2)=cordi(i,2)
enddo


  open(4,file='complete displacement.txt',status='unknown')
  do i=1,nn
  write(4,10)i,coordt(i,1),coordt(i,2),ux_t(i),uy_t(i)
 enddo
   do i=1,ni
  write(4,10)i+nn,coordt(i+nn,1),coordt(i+nn,2),ux_t(i+nn),uy_t(i+nn)
 enddo 
 10 format(i12,6e20.6)
close(4)  

end
  