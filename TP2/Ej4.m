%% Teoria
%LA PSD ES LA TRANSFORMADA DE FOURIER DE LA AUTOCORRELACION
%Sx(f)=int(Rx(tau)*e(-j 2 pi f tau))dtau ---> TIEMPO CONTINUO

%s(ejw)= sum(Rx(m)*e(-jmw))  ----> TIEMPO DISCRETO

%SI la autocorrelacion la valuo en cero Rx(0)=var(x)= potencia

%El area abajo de la PSD es la potencia 
%var(x)=int(Sx(f)df) ---> TIEMPO CONTINUO
%var(x) = 1/(2 pi)* int(pi y -pi)(Sx(ejw))dw

clear all; close all; clc

tline = 1:100e3;
rbg = wgn(length(tline),1,3)+1;     %Ruido blanco gaussiano
% hist(rbg,100)
% plot(rbg);

media = mean(rbg);
varianza = var(rbg);

%% Ruido Blanco Gaussiano Filtrado y(x)

h = 1/25*ones(1,25);
rbgf = filter(h,1,rbg);

mediaf = mean(rbgf);
varianzaf = var(rbgf);

%% PSD del Ruido Blanco Gaussiano x(n)

NFFT = 1024;
[pxx,freqv] = pwelch(rbg(:)-media,[],0,NFFT,'twoside');
figure;
plot(freqv,pxx);
grid on;
xlabel('Discrete Frequency [rad]');
ylabel('PSD [V^2/Hz]');


%% Potencia RBG 

potencia = mean(pxx)*2*pi;

% Integral de la PSD del RBG
potenciaPSD = sum(pxx)*(2*pi/NFFT);

%% PSD del RBG Filtrado y(n)

[pxxf,freqvf] = pwelch(rbgf(:)-media,[],0,NFFT,'twoside');
figure;
plot(freqvf,pxxf);
grid on;
xlabel('Discrete Frequency [rad]');
ylabel('PSD [V^2/Hz]');

%% Potencia RBGF

varF = mean(pxxf)*2*pi;

% Integral de la PSD del RBGF
potenciaPSDF = sum(pxxf)*(2*pi/NFFT);

%% SY(f) = Sx(f)|H(f)|'2

%Defino omega
NPOINTS = 25;
% omega = -pi:2*pi/NPOINTS:pi-2*pi/NPOINTS;
omega = [-1:2/NPOINTS:1-2/NPOINTS].*pi;

NFFT = 1024;

[Sxx, freqvs] = pwelch(rbg-media, [],0, NFFT, 'twoside', 'centered');
plot(freqvs, Sxx); grid on;

figure;
[Syy, freqvsy] = pwelch(rbgf-media, hanning(NFFT/2), 0,NFFT, 'twoside', 'centered');
plot(freqvsy, Syy); grid on;

Hf = fftshift(fft(h));
figure;
plot(omega, abs(Hf).^2); grid on;


Syy2 = Sxx.*Hf.*conj(Hf);
figure;
plot(Syy2);grid on;
