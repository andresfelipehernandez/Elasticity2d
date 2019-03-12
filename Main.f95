module vars
integer a1,nn,ne,npg,ni,nt
real*8 E,lambda,nu,mu,Sy,d
real*8,allocatable :: coord(:,:),bc(:),intq(:,:),intt(:,:),ma(:,:),u_b(:),t_b(:),cordi(:,:),u_i(:),R_t(:),coordt(:,:)
real*8,allocatable :: Q(:,:),H(:,:),G(:,:),ux_b(:),uy_b(:),tx_b(:),ty_b(:),ux_i(:),uy_i(:),ux_t(:),uy_t(:),despl(:,:),esfue(:,:)
real*8,allocatable::b(:)
integer,allocatable :: conect(:,:),bcc(:),indx(:)
endmodule

!Programa

program Mainprogram

use vars

!!! read boundary mesh and internal mesh
call inputdata
call inter

!!! Parameter
call parameters
!!!!

npg=8 !! Gauss Points
call boundarysol
call assembly
allocate(indx(2*nn))
!!!  solve system Ax=b
call ludcmp(ma,2*nn,2*nn,indx,d)
call lubksb(ma,2*nn,2*nn,indx,b)

call disassembly 
!$$$$$$ !!! boundary calcules
call  boundarycalc
 
!!! internal displacements
call internalsolve
call internalcalc
call totu
nt=nn+ni
call mecanic

stop
end