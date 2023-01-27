clear all
close all
clc
%% Parametros Generales

Nos = 16;           %Tasa de sobremuestreo
tau = 7.5e-9;       %Ancho del pulso [s]

fs = 1/tau*Nos;     %Frecuencia de muestreo [1/s]
Po = 5e3;           %Potencia del pulso [W]
Ao = sqrt(Po);      %Amplitud del pulso
Lpulse = tau*fs;    %Longitud del pulso
f0 = 77e9;          %Frec. de la portadora
c = 3e8;            %Velocidad de la luz [m/s]
Gt = 35;            %Ganancia direccional de antena TX [dBi]
Gr = 35;            %Ganancia direccional de antena RX [dBi]
lambda = c/f0;      %Long. de onda [m]
cross_s = 1;        %Radar cross section [m^2]
R = 1250;           %Rango o distancia [m]
rango_max = 2500;   %Rango maximo del radar [m]

fe_gain = 5;
input_voltage_noise_density = 0e-9; %[V/sqrt(z)]
No = input_voltage_noise_density^2; % PSD del ruido one-side
noise_power = No*fs;

%% Pulso
x_t = [ones(round(Lpulse),1); zeros(Nos,1)]*Ao; %Senal transmitida

%% Parametros del Canal
GtLin = 10^(Gt/10); %Linealizo Gt
GrLin = 10^(Gr/10); %Linealizo Gr

Range_Res = c*tau/2;    %Resolucion de rango
Delta_R = Range_Res*1.1;    %Distancia del segundo target

%Atenuacion
alphaLin0 = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*R^4)); %Atenuacion lineal
alplhadBi0 = 20*log10(alphaLin0); %Atenuacion en dBi

alphaLin1 = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*(R+Delta_R)^4)); %Atenuacion lineal
alplhadBi1 = 20*log10(alphaLin1); %Atenuacion en dBi

%Delay
delay0 = 2*R/c;                      %Retardo temporal
delayDiscreto0 = round(delay0*fs);    %Se lo lleva a un numero entero de muestras
delayReal0 = delayDiscreto0/fs;

delay1 = 2*(R+Delta_R)/c;             %Retardo temporal
delayDiscreto1 = round(delay1*fs);    %Se lo lleva a un numero entero de muestras
delayReal1 = delayDiscreto1/fs;

%Phase change
phase0 = exp(1j*2*pi*f0*delayReal0);   %Cambio de fase
phase1 = exp(1j*2*pi*f0*delayReal1);   %Cambio de fase

%% Modelado del Canal

h_t0 = phase0.*alphaLin0.*[zeros(delayDiscreto0,1); x_t];   % Salida del canal
h_t1 = phase0.*alphaLin1.*[zeros(delayDiscreto1,1); x_t];   % Misma fase

%% Potencias

% Comprobacion de la atenuacion en dB
Ptx = 10*log10(Po/1e-3);
Prx0 = 10*log10(max(abs(h_t0).^2)/1e-3);
atenuacion_de_canal_medida0 = Prx0-Ptx;

%% Rango maximo

muestras_max_range = round(2*rango_max/c*fs);
muestras0 = length(h_t0);
zeros_left0 = muestras_max_range-muestras0;

muestras1 = length(h_t1);
zeros_left1 = muestras_max_range-muestras1;

H_t0 = [h_t0; zeros(zeros_left0,1)];
H_t1 = [h_t1; zeros(zeros_left1,1)];

H_t = H_t0+H_t1;        %Salida del canal con ambos targets

%% Front-end
y_mf_accum = 0;
Nexp = 500;
Xi_range = zeros(Nexp,1);

noise = sqrt(noise_power/2)*randn(size(H_t)) + 1j.*sqrt(noise_power/2)*randn(size(H_t));
fe_output = fe_gain*H_t+noise;

    r_t = fe_output;            % Entrada del MF
    h_mf = conj(x_t(end:-1:1)); % Match Filter
    y_mf = conv(h_mf, r_t);     % Salida del detector

%Lineas de tiempo y rango
tline = 1/fs*(0:length(y_mf)-1);
rline = tline*c/2;

%Graficos
hold on
plot(rline, abs(y_mf).^2); grid on;
title('Potencia salida MF'); xlabel('Rango[m]');