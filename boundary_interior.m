
close all                                                                  
clear all                                                                   
clc                                                                         

%--------------------------------------------------------------------------
%Carga de datos
nodos=load('Nodes.txt');                                                    
nx=nodos(:,2);                                                              
ny=nodos(:,3);                                                              

 nodosi=load('Internal nodes.txt');                                         
 nxi=nodosi(:,2);                                                           
 nyi=nodosi(:,3);                                                            

%Gráfica de la malla de la frontera
figure(1)                                                                   
set(1,'NumberTitle','off','Name','Malla de la frontera',...
    'Position',[50 075 560 420])                                           
plot(nx,ny,'ok')                                                            
xlabel('x (cm)')                                                             
ylabel('y (cm)') 
title('Boundary Nodes')
axis equal                                                                 
grid on                                                                    

figure(2)                                                                   
set(2,'NumberTitle','off','Name','Internal mesh',...
    'Position',[100 100 560 420])                                           
plot(nxi,nyi,'.k')                                                          
xlabel('x (cm)')                                                            
ylabel('y (cm)')
title('Internal Nodes')
axis equal                                                                  
grid on                                                                    
