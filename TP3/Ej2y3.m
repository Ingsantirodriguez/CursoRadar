clear all
% close all
clc
%% Parámetros Generales

Nos = 16;           %Tasa de sobremuestreo
tau = 2.5e-9;       %Ancho del pulso [s]
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
R = 1760;           %Rango o distancia [m]
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

for R=100:100:rango_max

    %Atenuacion
    alphaLin = sqrt((GtLin*GrLin*lambda^2*cross_s)/((4*pi)^3*R^4)); %Atenuación lineal
    alplhadBi = 20*log10(alphaLin); %Atenuacion en dBi

    %Delay
    delay = 2*R/c;                      %Retardo temporal
    delayDiscreto = round(delay*fs);    %Se lo lleva a un número entero de muestras
    delayReal = delayDiscreto/fs;

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
    Xi_range = zeros(Nexp,1);

    snr_teo(ii) = (fe_gain*max(abs(H_t)))^2*tau/No;
    snr_teo_dB(ii) = 10*log10(snr_teo(ii));

    for Nexp=1:Nexp
        noise = sqrt(noise_power/2)*randn(size(H_t)) + 1j.*sqrt(noise_power/2)*randn(size(H_t));
        fe_output = fe_gain*H_t+noise;

        r_t = fe_output;            % Entrada del MF
        h_mf = conj(x_t(end:-1:1)); % Match Filter
        y_mf = conv(h_mf, r_t);     % Salida del detector

        if y_mf_accum == 0
            y_mf_accum = abs(y_mf).^2/Nexp;
        else
            y_mf_accum = abs(y_mf).^2/Nexp + y_mf_accum;   %Promedio
        end

        [M,I] = max(abs(y_mf).^2);
        Xi_time(Nexp) = I;
        max_index = I*rango_max/length(y_mf);
        Xi_range(Nexp) = max_index;
        Xi_olf_range = remove_outliers(Xi_range, R, 50);
        
        hold on
        plot(y_mf_accum)

    end

    %Lineas de tiempo y rango
    tline = 1/fs*(0:length(y_mf)-1);
    rline = tline*c/2;

    %SNR
    [PRX_peak,i] = max(y_mf_accum);
    if R<1000                       %Se calcula el piso de ruido dependiendo del rango.
        noise_floor = mean(y_mf_accum(i+1000:end));
    else
        noise_floor = mean(y_mf_accum(1:i-1000));
    end

    snr_comp(ii) = (PRX_peak-noise_floor)/noise_floor;  % Computada
    snr_comp_dB(ii) = 10*log10(snr_comp(ii));

    precision(ii) = std(Xi_range);  % Con outliers
    precision2(ii) = std(Xi_olf_range); % Sin outliers

    if precision(ii)<1e-2       % Con outliers
        precision(ii) = [];
    end
    
    if precision2(ii)<1e-2      % Sin outliers
        precision2(ii) = [];
    end

    ii=ii+1;
    ii
end
%% Gráficos

% Precisión con y sin outliers
rangeq = [100:100:rango_max];
hold on;
semilogy(rangeq, precision);grid on; % Se grafica la precisión a distintos rangos con el eje y en escala log
xlabel('Rango[m]'); ylabel('Presicion');
% hold on
% semilogy(rangeq, precision2);grid on;title("Sin outliers")
legend('\tau = 7.5[nseg]', '\tau = 5[nseg]', '\tau = 2.5[nseg]');

% SNR Computada y Teórica
figure
semilogy(rangeq, snr_comp_dB);grid on;title("SNR")
hold on
semilogy(rangeq, snr_teo_dB);grid on;title("SNR Teórica")
legend({'SNR Computada','SNR Teórica'},'Location', 'northeast')