subroutine mecanic
use vars
real*8::xm(nt),ym(nt),u(nt),v(nt),br(nt+1,nt+1),dfu(nt+1,nt+1),dfv(nt+1,nt+1)
real*8::alfau(nt+1),alfav(nt+1),du(nt+1),dv(nt+1),duy(nt+1),dvx(nt+1),rxx(nt+1),rxy(nt+1),ryy(nt+1)
real*8::m(nt+1,nt+1),exx(nt),eyy(nt),exy(nt),svm(nt)
real*8 det
integer st
allocate(esfue(nt,5))

! initialize vectors
alfau=0
alfav=0
du=0
dv=0
duy=0
dvx=0
rxx=0
rxy=0
ryy=0
exx=0
eyy=0
exy=0
rxx=0
ryy=0
rxy=0
svm=0



do i=1,nt
xm(i)=coordt(i,1)
ym(i)=coordt(i,2)
u(i)=despl(i,1)
v(i)=despl(i,2)
enddo

! Rbf calculation
call fn(nt,xm,ym,br)
call GINV(nt+1,br,m,0,1,det,st)
call mm(m,nt,u,nt,alfau)
call mm(m,nt,v,nt,alfav)
call dfn(nt,xm,ym,dfu,dfv)
! Derivatives Calculation
call mm(dfu,nt,alfau,nt,du)
call mm(dfv,nt,alfau,nt,duy)
call mm(dfu,nt,alfav,nt,dvx)
call mm(dfv,nt,alfav,nt,dv)

 open(4,file='Displacement.dat',status='unknown')
  do i=1,nt
  write(4,*)xm(i),ym(i),du(i),duy(i),dvx(i),dv(i)
 enddo
 close(4)

do i=1,nt
exx(i)=du(i)
eyy(i)=dv(i)
exy(i)=0.5*(duy(i)+dvx(i))
rxx(i)=2*mu*nu/(1-2*nu)*(du(i)+dv(i))+2*mu*(du(i))
ryy(i)=2*mu*nu/(1-2*nu)*(du(i)+dv(i))+2*mu*(dv(i))
rxy(i)=mu*(duy(i)+dvx(i))
svm(i)=dsqrt(rxx(i)**2+rxy(i)**2-ryy(i)*ryy(i)+3*rxy(i)**2)
enddo

 open(4,file='Strain_stress.dat',status='unknown')
 do i=1,nt
 write(4,*)i,exx(i),eyy(i),exy(i),rxx(i),ryy(i),rxy(i)
 enddo
 close(4)



end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        funcion 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fn(n,x,y,br)
integer::n,i,j
real*8::x(n),y(n),r(n,n),f(n,n),br(n+1,n+1)
br=1.0
br(n+1,n+1)=0
do i=1,n
  do j=1,n
    r(i,j)=dsqrt((x(j)-x(i))**2+(y(j)-y(i))**2)
    f(i,j)=1+r(i,j)
    br(i,j)=f(i,j)
  enddo
enddo
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! derivadas de la función
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfn(n,x,y,dfxx,dfyy)
integer::n,i,j
real*8::x(n),y(n),r(n,n),dfx(n,n),dfy(n,n),dfxx(n+1,n+1),dfyy(n+1,n+1)
dfxx=0.0d0
dfyy=0.0d0
do i=1,n
  do j=1,n
	r(i,j)=dsqrt((x(j)-x(i))**2+(y(j)-y(i))**2)
    if (r(i,j).eq.0) then
      dfx(i,j)=0
      dfy(i,j)=0
      else
    dfx(i,j)=(x(j)-x(i))/r(i,j)
    dfy(i,j)=(y(j)-y(i))/r(i,j)
  
    dfxx(i,j)=dfx(i,j)
    dfyy(i,j)=dfy(i,j)
    
    endif
  enddo
enddo
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiplicar matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mm(a,n,b,nb,d)
integer::n,nb
real*8::a(n,nb),b(nb),d(n)
d=matmul(a,b)
return 
end


