clear all
% close all
clc

%% Parámetros Generales

Nos = 16;           %Tasa de sobremuestreo
tau = 5e-9;         %Ancho del pulso [s]
fs = (1/tau)*Nos;   %Frecuencia de muestreo [1/s]
Po = 500;           %Potencia del pulso [W]
Ao = sqrt(Po);      %Amplitud del pulso
Lpulse = tau*fs;    %Longitud del pulso
f0 = 77e9;          %Frec. de la portadora
c = 3e8;            %Velocidad de la luz [m/s]
Gt = 20;            %Ganancia direccional de antena TX [dBi]
Gr = 20;            %Ganancia direccional de antena RX [dBi]
lambda = c/f0;      %Long. de onda [m]
cross_s = 1;        %Radar cross section [m^2]
R = 2000;           %Rango o distancia [m]
rango_max = 2500;   %Rango máximo del radar [m]

%% Pulso

x_t = [ones(Lpulse,1); zeros(Nos,1)]*Ao; %Senal transmitida
% figure(1); 
% plot(x_t); grid on; ylim([0 Ao+1]); title('Señal transmitida x_t');
% figure(2);
% plot(abs(real(x_t)).^2); grid on; title('Potencia de la señal transmitida');


%% Parámetros del Canal
GtLin = 10^(Gt/10); %Linealizo Gt
GrLin = 10^(Gr/10); %Linealizo Gr

%Atenuacion
alphaLin = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*R^4)); %Atenuación lineal
alplhadBi = 20*log10(alphaLin); %Atenuacion en dBi

%Delay
delay = 2*R/c;                  %Retardo temporal
delayDiscreto = round(delay*fs);%Se lo lleva a un número entero de muestras
delayReal = delayDiscreto/fs; 

%Phase change
phase = exp(1j*2*pi*f0*delayReal);   %Cambio de fase
  
%% Modelado del Canal

h_t = phase.*alphaLin.*[zeros(delayDiscreto,1); x_t];
% figure
% plot(real(h_t)); grid on; title('Fase de la señal recibida')
% hold on;
% plot(imag(h_t)); legend('Real', 'Imaginaria', 'Location', 'northwest')
% figure;
% plot(abs(h_t).^2); title('Potencia del canal'); grid on;


%% Potencias

% Comprobación de la atenuación en dB
Ptx = 10*log10(Po/1e-3);
Prx = 10*log10(max(abs(h_t).^2)/1e-3);
atenuacion_de_canal_medida = Prx-Ptx;

%% Rango máximo

max_range = 2500;
muestras_max_range = round(2*max_range/c*fs);
muestras = length(h_t);
zeros_left = muestras_max_range-muestras;

H_t = [h_t; zeros(zeros_left,1)];
hold all
plot(abs(H_t).^2); title('Señal recibida respecto del rango máximo');grid on;
%legend('R = 500', 'R = 1000', 'R = 1500', 'R = 2000');

%% Front-end

fe_gain = 1;
noise_power = 1e-26*fs;     %noise one-side
noise = sqrt(noise_power/2)*randn(size(H_t)) + ...
     1j.*sqrt(noise_power/2)*randn(size(H_t));
fe_output = fe_gain*H_t+noise;

figure
plot(abs(fe_output).^2); grid on; title('Front-end');

%% Receptor óptimo (Match Filter)

r_t = fe_output;            % Entrada del MF
h_mf = conj(x_t(end:-1:1)); % Match Filter
y_mf = conv(h_mf, r_t);     % Salida del detector

%Gráfico en función del tiempo
figure
tline = 1/fs*(0:length(y_mf)-1);
plot(tline,abs(y_mf)); grid on; title('Salida del MF');
hold all; xlabel('tiempo [s]'); 

rline = tline*c/2;
figure
plot(rline,abs(y_mf)); grid on; title('Salida del MF');
hold all; xlabel('rango [m]'); 





