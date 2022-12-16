%Ejercicio 6

clc; close all; clear all;

t = 1:50e6;
rbg = sqrt(1)*randn(size(t));
plot(rbg);grid on;

%% Media y Varianza
media = mean(rbg);
varianza = var(rbg);

%% PSD

NFFT = 1024;
[pxx, freq] = pwelch(rbg-media, [], 0, NFFT, 50e6, 'power');
figure
plot(freq, pxx);grid on;