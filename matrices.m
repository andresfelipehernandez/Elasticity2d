close all                                                                   %Cierra ventanas previas
clear all                                                                   %Borra variables previas
clc                                                                         %Borra la ventana de comandos

%--------------------------------------------------------------------------
%Datos nodos de la frotera
nodos=224;                                                                  %Número de nodos en la frontera
x=1:1:nodos;                                                                %Vector de prueba
y=x;                                                                        %Vector de prueba
[X,Y]=meshgrid(x,y);                                                        %Genera una malla reticulada en XY. lim_min:espaciamiento:lim_máx

%--------------------------------------------------------------------------
%Carga de datos
fid=fopen('matrixH.txt.','r');                                              %Abre el archivo de texto, sólo lectura
for k=1:nodos                                                               %Inicia el ciclo en función del número de nodos
string=textscan(fid,'%30c',1);                                              %Lee la etiqueta de encabezado
fila=textscan(fid,'%f', nodos, 'delimiter','\t');                           %Lee los 108 elementos de la fila y los mete en una celda
fila = fila{:};                                                             %Desencapsula el vector contenido en la celda
H(k,:)=fila;                                                                %Arma la matriz Q fila por fila
end                                                                         %Termina el ciclo
fclose(fid);                                                                %Cierra el archivo de texto

fid=fopen('matrixG.txt.','r');                                              %Abre el archivo de texto, sólo lectura
for k=1:nodos                                                               %Inicia el ciclo en función del número de nodos
string=textscan(fid,'%30c',1);                                              %Lee la etiqueta de encabezado
fila=textscan(fid,'%f', nodos, 'delimiter','\t');                           %Lee los 108 elementos de la fila y los mete en una celda
fila = fila{:};                                                             %Desencapsula el vector contenido en la celda
G(k,:)=fila;                                                                %Arma la matriz T fila por fila
end                                                                         %Termina el ciclo
fclose(fid);                                                                %Cierra el archivo de texto

%--------------------------------------------------------------------------
%Gráfica matriz Q
figure(1)
set(1,'NumberTitle','off','Name','Matriz Q',...
    'Position',[200 125 560 420])%Cambia propiedades de la ventana de la gráfica
surf(X,Y,H,'edgecolor','none')                                                            %Crea gráfica 3D con las coordenas X, Y, Q
xlabel('x Nodes')
ylabel('y Nodes')
zlabel('Matrix H')
colorbar

%Gráfica matriz T
figure(2)
set(2,'NumberTitle','off','Name','Matriz T',...
    'Position',[250 150 560 420])                                           %Cambia propiedades de la ventana de la gráfica
surf(X,Y,G,'edgecolor','none')                                                           %Crea gráfica 3D con las coordenas X, Y, T
xlabel('x Nodes')
ylabel('y Nodes')
zlabel('Matrix G')
colorbar
