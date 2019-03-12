subroutine boundarysol

!H matrix

use vars

real*8 x(1,2),y(3,2),phi(3),dphi(3),psi(npg),ypsi(2),dypsi(2),xb(2),r,jac,no(2),pi,wg(npg)
real*8 drdx,drdy,drdn,uxx,uxy,uyx,uyy,txx,txy,tyx,tyy,us(2,2),ts(2,2)
allocate(H(2*nn,2*nn),G(2*nn,2*nn))

pi=3.1415926535897932384626433832795

psi(1)=-0.960289856497537
psi(2)=-0.796666477413627
psi(3)=-0.525532409916329
psi(4)=-0.183434642495650
psi(5)=0.183434642495650
psi(6)=0.525532409916329
psi(7)=0.796666477413627
psi(8)=0.960289856497536

wg(1)=0.101228536290376
wg(2)=0.222381034453375
wg(3)=0.313706645877887
wg(4)=0.362683783378362
wg(5)=0.362683783378362
wg(6)=0.313706645877887
wg(7)=0.222381034453374
wg(8)=0.101228536290376    


H=0.0
G=0.0

do j=1,nn
  x(1,1:2)=coord(j,1:2)                                           
  do k=1,ne                                                                 
    do ii=1,3                                         
      y(ii,1:2)=coord(conect(k,ii),1:2)
    enddo
    
    do i=1,npg
      phi(1)=-0.5*psi(i)*(1-psi(i))
      phi(2)=(1+psi(i))*(1-psi(i))
      phi(3)=0.5*psi(i)*(1+psi(i))
      dphi(1)=psi(i)-0.5
      dphi(2)=-2*psi(i)
      dphi(3)=psi(i)+0.5
      
      ypsi=0.0 ; dypsi=0.0
      
      !Creates parametrics function that describes element in Gauss points terms
      
      do ii=1,3                                             
        ypsi(1:2)=ypsi(1:2)+y(ii,1:2)*phi(ii)
        dypsi(1:2)=dypsi(1:2)+y(ii,1:2)*dphi(ii)
      enddo
      
      !Factors to find integral.

      xb(1:2)=x(1,1:2)-ypsi(1:2)                      
      r=dsqrt(xb(1)**2+xb(2)**2)
      jac=dsqrt(dypsi(1)**2+dypsi(2)**2)
      no(1)=dypsi(2)/jac
      no(2)=-dypsi(1)/jac
      
	  !mechanical components
      
      do ii=1,3
       	drdx=(y(ii,1)-x(1,1))/r
        drdy=(y(ii,2)-x(1,2))/r
        drdn=drdx*no(1)+drdy*no(2)                                  
        uxx=1/(8*pi*mu*(1-nu))*((3-4*nu)*log(1/r)+drdx**2)
        uxy=1/(8*pi*mu*(1-nu))*drdx*drdy
        uyx=1/(8*pi*mu*(1-nu))*drdx*drdy
        uyy=1/(8*pi*mu*(1-nu))*((3-4*nu)*log(1/r)+drdy**2)
        txx=-1/(4*pi*(1-nu)*r)*drdn*((1-2*nu)+2*drdx**2)      
        txy=-1/(4*pi*(1-nu)*r)*(2*drdx*drdy*drdn+(1-2*nu)*(drdy*no(1)-drdx*no(2)))      
        tyx=-1/(4*pi*(1-nu)*r)*(2*drdx*drdy*drdn-(1-2*nu)*(drdy*no(1)-drdx*no(2)))      
        tyy=-1/(4*pi*(1-nu)*r)*drdn*((1-2*nu)+2*drdy**2)

        !! ustar and tstar---------------------------------------------------------------------
        us(1,1)=uxx
        us(1,2)=uxy
        us(2,1)=uyx
        us(2,2)=uyy
        
        ts(1,1)=txx
        ts(1,2)=txy
        ts(2,1)=tyx
        ts(2,2)=tyy

       H(2*j-1:2*j,2*conect(k,ii)-1:2*conect(k,ii))=H(2*j-1:2*j,2*conect(k,ii)-1:2*conect(k,ii))+ts*jac*wg(i)*phi(ii) 
       G(2*j-1:2*j,2*conect(k,ii)-1:2*conect(k,ii))=G(2*j-1:2*j,2*conect(k,ii)-1:2*conect(k,ii))+us*jac*wg(i)*phi(ii) 

      enddo
    enddo
  enddo
  
enddo

do i=1,2*nn
  open(11,file='matrixG.txt',status='unknown')
  write(11,*)'gama number',i,G(i,:)
  enddo
close(11)

!! Helps calculate the diagonal coefficients indirectly
    
do i=1,nn                        
  H(2*i-1:2*i,2*i-1:2*i)=0.0
  do j=1,nn
    if (i/=j) then
       H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)-H(2*i-1:2*i,2*j-1:2*j)
    end if
  end do
end do  

do i=1,2*nn
  open(9,file='matrixH.txt',status='unknown')
  write(9,*)'gama number',i,H(i,:)
  enddo
close(9)

return
end