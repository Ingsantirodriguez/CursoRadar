% Generar ruido blanco gaussiano para simular el ruido de un amplificador 
% analogico.
% Ancho de banda de la senal de datos = 1GHz.
% Factor de sobremuestreo = 8. => fs=8GHz.
% Potencia del ruido termico medido en el ancho de banda de 8GHz es de
% 100mW.

clear all; close all; clc

fs8 = 8e6; %frec de muestreo 8 [GHz]
t8 = 0:1/fs8:1;

fs16 = 16e6; %frec de muestreo 16 [GHz]
t16 = 0:1/fs16:1;

%% Genero el ruido con randn (8GHz)

rbg = randn(size(t8))*sqrt(100e-3);
plot(rbg); grid on;

%% Media y varianza
media_rbg = mean(rbg);
varianza_rbg = var(rbg); % Potencia = 100 mW

%% PSD del ruido blanco gaussiano (rbg)
NFFT = 2048;
%En el parametro fs colocamos nuestra frecuencia de muestreo.
[pxx, freq] = pwelch(rbg-media_rbg, hanning(NFFT/2), 0, NFFT, fs8, 'twoside'); 
figure;
plot(freq, pxx); grid on; title('PSD')
% hold on;

mediapxx8 = mean(pxx);
potenciaPXX8 = mediapxx8*2*pi;   %potencia de la psd
pot8 = sum(pxx)*(2*pi/NFFT); %integral de la psd

%% %Genero el ruido con randn (16 GHz)

rbg = randn(size(t16))*sqrt(100e-3);
figure
plot(rbg); grid on;

%% PSD del ruido blanco gaussiano (16 GHz)
NFFT = 2048;
%En el parametro fs colocamos nuestra frecuencia de muestreo.
[pxx, freq] = pwelch(rbg-media_rbg, hanning(NFFT/2), 0, NFFT, fs16, 'twoside');
%figure;
plot(freq, pxx); grid on; title('PSD');
legend('fs = 8 GHz', 'fs = 16 GHz')
mediapxx16 = mean(pxx);
potenciaPXX16 = mediapxx16*2*pi; %potencia de la psd
pot16 = sum(pxx)*(2*pi/NFFT); %integral de la psd

%% Ruido Complejo Blanco Gaussiano

% x = randn(size(t8))*sqrt(50e-3);   %Procesos estocasticos independientes
% xc = randn(size(t8))*sqrt(50e-3);
% rbgc = x + 1j*xc;
% plot(rbgc); grid on;

%% Media y varianza (rbgc)
media_rbgC = mean(rbgc);
varianza_rbgC = var(rbgc); % Potencia = 100 mW por ser procesos independientes

%% PSD del ruido blanco gaussiano complejo (rbgc)
NFFT = 2048;
%En el parametro fs colocamos nuestra frecuencia de muestreo.
[pxx, freq] = pwelch(rbgc-media_rbgC, hanning(NFFT/2), 0, NFFT, fs8, 'twoside'); 
figure;
plot(freq, pxx); grid on; title('PSD')
hold on

mediapxx8C = mean(pxx);
potenciaPXX8C = mediapxx8*2*pi;   %potencia de la psd
pot8C = sum(pxx)*(2*pi/NFFT); %integral de la psd












