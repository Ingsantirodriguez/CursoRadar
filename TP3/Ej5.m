clear all
close all
clc
%% Parámetros Generales

Nos = 16;           %Tasa de sobremuestreo
tau = 5e-9;       %Ancho del pulso [s]
fs = (1/tau)*Nos;   %Frecuencia de muestreo [1/s]
Po = 5e3;          %Potencia del pulso [W]
Ao = sqrt(Po);      %Amplitud del pulso
Lpulse = tau*fs;    %Longitud del pulso
f0 = 77e9;          %Frec. de la portadora
c = 3e8;            %Velocidad de la luz [m/s]
Gt = 35;            %Ganancia direccional de antena TX [dBi]
Gr = 35;            %Ganancia direccional de antena RX [dBi]
lambda = c/f0;      %Long. de onda [m]
cross_s = 1;        %Radar cross section [m^2]
R = 1500;           %Rango o distancia [m]
rango_max = 2500;   %Rango máximo del radar [m]

fe_gain = 5;
input_voltage_noise_density = 0.4e-9; %[V/√Hz]
No = input_voltage_noise_density^2; %noise one-side
noise_power = No*fs;

%% Pulso

x_t = [ones(round(Lpulse),1); zeros(Nos,1)]*Ao; %Senal transmitida

%% Parámetros del Canal
GtLin = 10^(Gt/10); %Linealizo Gt
GrLin = 10^(Gr/10); %Linealizo Gr

ii = 1;

for R=1800:100:rango_max

    %Atenuacion
    alphaLin = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*R^4)); %Atenuación lineal
    alplhadBi = 20*log10(alphaLin); %Atenuacion en dBi

    %Delay
    delay = 2*R/c;                      %Retardo temporal
    delayDiscreto = round(delay*fs);    %Se lo lleva a un número entero de muestras
    delayReal = delayDiscreto/fs;
    
    %Posición
    pos_real = delayDiscreto+Lpulse;

    %Phase change
    phase = exp(1j*2*pi*f0*delayReal);   %Cambio de fase

    %% Modelado del Canal

    h_t = phase.*alphaLin.*[zeros(delayDiscreto,1); x_t];

    %% Potencias

    % Comprobación de la atenuación en dB
    Ptx = 10*log10(Po/1e-3);
    Prx = 10*log10(max(abs(h_t).^2)/1e-3);
    atenuacion_de_canal_medida = Prx-Ptx;

    %% Rango máximo

    muestras_max_range = round(2*rango_max/c*fs);
    muestras = length(h_t);
    zeros_left = muestras_max_range-muestras;

    H_t = [h_t; zeros(zeros_left,1)];

    %% Front-end
    y_mf_accum = 0;
    Nexp = 1000;

    COI = 1+ceil(pos_real/Nos); % Índice de la celda de interés
    vector_COI = zeros(Nexp,1); % Vector de valores de las COI
    vector_inter = zeros(Nexp,1); % Vector de valores de las demás celdas

    for Nexp=1:Nexp
        noise = sqrt(noise_power/2)*randn(size(H_t)) + 1j.*sqrt(noise_power/2)*randn(size(H_t));
        fe_output = fe_gain*H_t+noise;

        r_t = fe_output;            % Entrada del MF
        h_mf = conj(x_t(end:-1:1)); % Match Filter
        y_mf = conv(h_mf, r_t);     % Salida del detector
        
        offset = mod(pos_real-1,Nos); % Factor que varía para encontrar la mitad de la celda
%         plot(abs(y_mf));
        y_mf_celdas = y_mf(1+offset:Nos:end); % Vector que contiene los centros de las celdas
%         hold on
%         plot(abs(y_mf_celdas));
%         plot([1+offset:Nos:length(y_mf)], abs(y_mf_celdas), 'o');
        
        vector_COI(Nexp) = y_mf_celdas(COI);    % Le asigno el valor de la COI para este exp.
        
        aux = y_mf_celdas;  % Creo una copia
        aux(COI) = [];      % Quito la COI
        if vector_inter==0  % Sólo primera iteración
            vector_inter = zeros(Nexp, length(aux)); % Le asigno la cantidad de celdas no COI por exp.
            vector_inter(Nexp,:) = aux; % Al exp. actual lo guardo en la matriz
        else
            vector_inter(Nexp,:) = aux;
        end

    end

    %Lineas de tiempo y rango
    tline = 1/fs*(0:length(y_mf)-1);
    rline = tline*c/2;

    ii=ii+1;
    ii
end
%% Gráficos

rangeq = [100:100:rango_max];