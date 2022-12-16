clc; close all; clear all;

% Generacion de PAM2 
imin = 1;
imax = 2;
t = 1:100e3;
X_n = randi([imin, imax],[1, length(t)]);

for i=1:length(t)
    if X_n(i) == 2
        X_n(i) = -1;
    end
end

figure;
hist(X_n); grid on; xlim([-3;3]); title('X(n)');

%% Ruido blanco gaussiano

potencia = 0.1;
N_n = sqrt(potencia)*randn(size(t));
figure;
plot(N_n); grid on; title('Ruido Blanco Gaussiano N(n)');
figure;
hist(N_n, 100); grid on; title('N(n)');

%% Señal AWGN
Y_n = X_n + N_n;

figure;
hist(Y_n, 100); grid on; title('Y(n)');

%% Calculo de SNR

Ps = 1
Pn = var(N_n)

SNR = Ps/Pn

%% Toma de decision en el receptor

Y_decidido = Y_n;
Y_decidido(Y_n > 0) = 1;
Y_decidido(Y_n < 0) = -1;
figure;
hist(Y_decidido); grid on; xlim([-3;3]); title('Simbolos detectados');

F = 0;
M = length(t);
for i = 1:length(t)
    if X_n(i) ~= Y_decidido(i)
       F = F + 1; 
    end
end 

BER = F/M
