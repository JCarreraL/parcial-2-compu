%% Configuración de la conexión
clear
%Revisa el puerto al que se conectó el Arduino 
disp("Puertos disponibles:");
disp(serialportlist);

%% Configuración de la conexión
port = "COM6"; %Reemplaza COM4 con el puerto detectado en tu computador
baudRate = 9600; %El arduino está configurado en 9600 baudios
ArduinoUno = serialport(port, baudRate);
pause(3); %tiempo de espera para evitar errores

%% Su código aquí

sols_j = load("sols_j.mat").sols_j;
tiempos = zeros(length(sols_j),1);
tic;

for i = 1:length(sols_j)
    mensaje = sprintf("%.2f,%.2f,%.2f\n", sols_j(i,1), sols_j(i,2), sols_j(i,3));
    writeline(ArduinoUno, mensaje);
    pause(1)
    tiempos(i) = toc;
end
pause(2)

% USE LA SIGUIENTE ESTRUCTURA PARA MANDAR LOS ÁNGULOS ITERATIVAMENTE
% Mensaje en formato "theta1,theta2,theta3"

%% Al terminar regresa a la posición inicial y cierra la conexión

%posición inicial horizontal
theta1 = 20;
theta2 = 0;
theta3 = 0;

mensaje = sprintf("%.2f,%.2f,%.2f\n", theta1, theta2, theta3);
writeline(ArduinoUno, mensaje);

%%
%termina la conexion con estos comandos
flush(ArduinoUno);
delete(ArduinoUno);
clear ArduinoUno;