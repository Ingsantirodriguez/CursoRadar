clear all; close all; clc

tline = 1:100e3;
rbg = wgn(length(tline),1,3)+1;
h = 1/25*ones(1,25);
rbgf = filter(h,1,rbg);
figure
plot(rbg);
hold all
plot(rbgf);
legend('RBG','RBG Filtrado'); grid on;

%Media y varianza del ruido blanco
media = mean(rbg)
varianza = var(rbg)

%Media y esperanza de la señal filtrada (y)
mediaf = mean(rbgf)
esperanzaf = sum(rbgf)/length(tline)

%% Autocorrelacion
t = -99999:1:99999;
Rx = xcorr(rbg, 'biased');
figure
plot(t,Rx); grid on;% title('Autocorrelacion')
hold on

%% Autocovarianza
Cx = xcov(rbg, 'biased');
% figure
plot(t,Cx); grid on;% title('Autocovarianza')
legend('Autocorrelacion', 'Autocovarianza')

%% Autocorr Filtrado
t = -99999:1:99999;
Ry = xcorr(rbgf, 'biased');
figure
plot(t,Ry); grid on;% title('Autocorrelacion')
hold on

%% Funcion del filtro
t = -100023:1:100023;
hCI = conj(h(end:-1:1));
g = conv(hCI,h); 
Rxconv = conv(Rx,g);

% figure
plot(t,Rxconv); grid on;
legend('Ry(m)', 'Rx(m)*h(m)*h*(-m)')
