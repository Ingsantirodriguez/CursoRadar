clear all
% close all
clc
%% Parámetros Generales

Nos = 16;           %Tasa de sobremuestreo
tau = 5e-9;       %Ancho del pulso [s]
fs = (1/tau)*Nos;   %Frecuencia de muestreo [1/s]
Po = 2.5e3;          %Potencia del pulso [W]
Ao = sqrt(Po);      %Amplitud del pulso
Lpulse = tau*fs;    %Longitud del pulso
f0 = 77e9;          %Frec. de la portadora
c = 3e8;            %Velocidad de la luz [m/s]
Gt = 35;            %Ganancia direccional de antena TX [dBi]
Gr = 35;            %Ganancia direccional de antena RX [dBi]
lambda = c/f0;      %Long. de onda [m]
cross_s = 1;        %Radar cross section [m^2]
R = 1100;           %Rango o distancia [m]
rango_max = 2500;   %Rango máximo del radar [m]

fe_gain = 5;
input_voltage_noise_density = 0.8e-9; %[V/√Hz]
No = input_voltage_noise_density^2; %noise one-side
noise_power = No*fs;

% Pulso
x_t = [ones(round(Lpulse),1); zeros(Nos,1)]*Ao; %Senal transmitida

% Parámetros del Canal
GtLin = 10^(Gt/10); %Linealizo Gt
GrLin = 10^(Gr/10); %Linealizo Gr


% Atenuacion
alphaLin = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*R^4)); %Atenuación lineal
alplhadBi = 20*log10(alphaLin); %Atenuacion en dBi

% Delay
delay = 2*R/c;                      %Retardo temporal
delayDiscreto = round(delay*fs);    %Se lo lleva a un número entero de muestras
delayReal = delayDiscreto/fs;

% Posición
pos_real = delayDiscreto+Lpulse;

% Phase change
phase = exp(1j*2*pi*f0*delayReal);   %Cambio de fase

% Modelado del Canal

h_t = phase.*alphaLin.*[zeros(delayDiscreto,1); x_t];

% Potencias

% Comprobación de la atenuación en dB
Ptx = 10*log10(Po/1e-3);
Prx = 10*log10(max(abs(h_t).^2)/1e-3);
atenuacion_de_canal_medida = Prx-Ptx;

% Rango máximo

muestras_max_range = round(2*rango_max/c*fs);
muestras = length(h_t);
zeros_left = muestras_max_range-muestras;

H_t = [h_t; zeros(zeros_left,1)];

% Front-end
y_mf_accum = 0;
Nexp = 500;

COI = 1+ceil(pos_real/Nos); % Índice de la celda de interés
vector_COI = zeros(Nexp,1); % Vector de valores de las COI
vector_inter = zeros(Nexp,1); % Vector de valores de las demás celdas

%% Simulacion
ii = 0
for Nexp=1:Nexp
    noise = sqrt(noise_power/2)*randn(size(H_t)) + 1j.*sqrt(noise_power/2)*randn(size(H_t));
    fe_output = fe_gain*H_t+noise;

    r_t = fe_output;            % Entrada del MF
    h_mf = conj(x_t(end:-1:1)); % Match Filter
    y_mf = conv(h_mf, r_t);     % Salida del detector
    
    offset = mod(pos_real-1,Nos); % Factor que varía para encontrar la mitad de la celda
    y_mf_celdas = y_mf(1+offset:Nos:end); % Vector que contiene los centros de las celdas
%     plot(abs(y_mf));
%     hold on
%     plot([1+offset:Nos:length(y_mf)], abs(y_mf_celdas), 'o');
    
    vector_COI(Nexp) = y_mf_celdas(COI);    % Le asigno el valor de la COI para este exp.
    
    aux = y_mf_celdas;  % Creo una copia
    aux(COI) = [];      % Quito la COI
    if vector_inter==0  % Sólo primera iteración
        vector_inter = zeros(Nexp, length(aux)); % Le asigno la cantidad de celdas no COI por exp.
        vector_inter(Nexp,:) = aux; % Al exp. actual lo guardo en la matriz
    else
        vector_inter(Nexp,:) = aux;
    end
    ii=ii+1
end

noise_samples = reshape(vector_inter, [], 1); %Convierto matriz en vector

disp('Simulaciones Terminadas')

% % Histogramas

[counts, centers] = hist(abs(noise_samples), 100);
figure
plot(centers,counts/length(noise_samples));
hold all
[counts, centers] = hist(abs(vector_COI), 100);
plot(centers,counts/length(vector_COI));grid on;

%% Thresholds

max_threshold = max(abs(vector_COI))*1.1;
Nthrs = 500;
thresholds = (max_threshold/Nthrs:max_threshold/Nthrs:max_threshold);
fn_vector = zeros(Nthrs,1);
fp_vector = zeros(Nthrs,1);
tp_vector = zeros(Nthrs,1);
tn_vector = zeros(Nthrs,1);

% Calculo de thresholds

X = (abs(vector_COI)>thresholds).'; % Ths = filas, exp = columnas
Y = (abs(noise_samples)>thresholds).';

disp('Thresholds calculados')

%% Calculo de ROC

PD = zeros(length(thresholds));
PFA = zeros(length(thresholds));

ths_stats = zeros(4,length(thresholds)); % TPs, FNs, FPs, TNs = filas en ese orden ; Ths = columnas

for i=1:1:length(thresholds)
    ths_stats(1,i) = sum(X(i,:)); % TPs para th i
    ths_stats(2,i) = length(X(i,:))-ths_stats(1,i); % FNs para th i
    ths_stats(3,i) = sum(Y(i,:)); % FPs para th i
    ths_stats(4,i) = length(Y(i,:))-ths_stats(3,i); % TNs para th i
    
    PD(i) = (ths_stats(1,i))/length(vector_COI);
    PFA(i) = (ths_stats(3,i))/length(noise_samples);
    i
end

% Gráficos ROC

figure
semilogx(PFA,PD);grid on;xlabel("PFA");ylabel("PD");title("ROC");legend('250 [m]','500 [m]','1000 [m]','1500 [m]');

snr_teo = (fe_gain*max(abs(H_t)))^2*tau/No;
snr_teo_dB = 10*log10(snr_teo);

[Pd_teo,Pfa_teo] = rocsnr(snr_teo_dB,SignalType='NonFluctuatingNonCoherent');
hold on
semilogx(Pfa_teo,Pd_teo);grid on;xlabel("PFA");ylabel("PD");legend('ROC Computada', 'ROC Teórica');