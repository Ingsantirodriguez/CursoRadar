clear all
% close all
clc

%% Parámetros del LiDAR
PTX = 1e-3*10^(15/10);  % Potencia transmitida [W] 15 dBm
chirp_bw = 2e9;         % Ancho de banda del chirp
max_range = 350;        % Rango máximo
c = 3e8;
Twait = 2*max_range/c;
q_elect = 1.6e-19;              % Carga del electrón, para calcular ruido
Tmeas = 2.5e-6;                 % Tiempo de medición
Tmod = Tmeas+Twait;             % Tiempo de modulación
chirp_slope = chirp_bw/Tmod;    % Pendiente del chirp

% Parámetros del canal
range = 200;            %m
ARX=pi*(2.5e-2/2)^2;    % Apertura de 1in de diametro
rho = 0.1;              % Reflectividad
lambda0=1550e-9;        % m
omega0 = 2*pi*3e8/lambda0; 
RPD = 0.7;              % A/W Responsitividad del diodo

% Muestreo
NOS = 4;
fs=NOS*chirp_bw;          % Frec. de muestreo de Matlab
Ncells = ceil(fs*Tmeas);
fs=Ncells/Tmeas;        % Frecuencia de muestreo
Ts = 1/fs;              % Período de muestreo

% Transmisor
Lsim = ceil(Tmod*fs);                       % Largo de la sim.
tline = Ts.*(0:Lsim-1)';                    % Linea de tiempo de modulación
insta_freq = chirp_slope.*tline;            % Frecuencia instantánea
insta_phase = 2*pi*cumsum(insta_freq).*Ts;  % Fase instantánea
chirp_tx = exp(1j*insta_phase);             % Chirp unitario

s_t = sqrt(PTX).*chirp_tx;                  % Señal transmitida

% Canal
tau = 2*range/c;              % Delay del canal
delay_samples = round(tau*fs);
real_tau = delay_samples*Ts;
real_range = real_tau*c/2;    % Rango real

power_gain = rho*ARX/(4*pi*range.^2);   % Potencia
atten = sqrt(power_gain);               % Atenuación
delta_phase = 2*pi*c/(lambda0*real_tau);         % Cambio de fase

ch_out = atten.*[zeros(delay_samples,1); s_t(1:end-delay_samples)].*exp(-1j*delta_phase);
% Salida del canal (atten, fase y delay)

%% Receptor
wait_samples = ceil(Twait*fs);
mixer = conj(ch_out.*conj(chirp_tx)); % Salida del detector (conv.) con delay
mixer_v = mixer(1+wait_samples:end);  % Salida del mixer ventaneada
ch_out_v = ch_out(1+wait_samples:end);% Salida del canal ventaneada

% Matched Filter FFT

Ncells = ceil(fs*Tmeas);    % numero de bines de la FFT minimo
FFT_NOS = 32;               % Sobremuestreo de la FFT
NFFT = FFT_NOS*Ncells;
fvec = (0:NFFT-1)*(fs/NFFT);            % Vector de frecuencia
max_range_FFT = fs*c/(chirp_slope*2);   % Máximo rango que puede ser representado por la FFT
rvec = (0:NFFT-1)*(max_range_FFT/NFFT); % Vector de rango

% División en celdas

fbeat = chirp_slope*2*real_range/c; % Frecuencia (teórica) del tono resultante
COI = ceil(fbeat*Tmeas);            % Índice de la celda de interés en fvec_dec
fbeat_index = round(fbeat/fs*NFFT)+1;  % Posición teórica (muestras) del fbeat en la transformada
fft_dec_offset = mod(fbeat_index-1,FFT_NOS);    % Posición dentro de la celda donde está el máximo
fvec_dec = fvec(1+fft_dec_offset:FFT_NOS:end);   % Vector de frecuencia decimado, con máximo en medio de la celda

% Experimentos
noise_power = q_elect/RPD*fs;     % Potencia del ruido
Nexp = 1000;
vector_COI = zeros(Nexp,1); % Vector de valores de las COI
vector_inter = zeros(Nexp,1); % Vector de valores de las demás celdas (ruido)
t=0;
for i=1:Nexp
    
    noise = sqrt(noise_power/2).*(randn(size(ch_out))+1j.*randn(size(ch_out))); % Señal de ruido
    noise_v = noise(1+wait_samples:end);    % Ruido Ventaneado
    fe_output = mixer_v + noise_v;          % Salida del detector con ruido

    y_mf = abs(fft(fe_output, NFFT)).^2;
    y_mf_dec = y_mf(1+fft_dec_offset:FFT_NOS:end);   % Vector y_mf decimado
    
    vector_COI(i) = y_mf_dec(COI);

    aux = y_mf_dec;
    aux(COI) = [];

    if vector_inter==0  % Sólo primera iteración
        vector_inter = zeros(Nexp, length(aux)); % Le asigno la cantidad de celdas no COI por exp.
        vector_inter(i,:) = aux; % Al exp. actual lo guardo en la matriz
    else
        vector_inter(i,:) = aux;
    end

    if y_mf_dec(COI) < max(y_mf_dec)
        t = t+1;
    end

    i

end

noise_samples = reshape(vector_inter, [], 1); %Convierto matriz en vector
noise_samples = noise_samples(1:length(noise_samples)/5); % Recorta
% vector de ruido para agilizar los cálculos
mean(noise_samples)


%% Thresholds

Nthrs = 1000;   % Cantidad de thresholds
max_threshold = max(vector_COI)*1.05;  % Threshold máximo
min_threshold = max_threshold/Nthrs;    % Threshold mínimo
thresholds = (min_threshold:min_threshold:max_threshold);
fn_vector = zeros(Nthrs,1);
fp_vector = zeros(Nthrs,1);
tp_vector = zeros(Nthrs,1);
tn_vector = zeros(Nthrs,1);

% Calculo de thresholds

X = (abs(vector_COI)>thresholds).'; % Ths = filas, exp = columnas
Y = (abs(noise_samples)>thresholds).';

disp('Thresholds calculados')

%% Calculo de ROC

PD = zeros(length(thresholds),1);   % PD para cada thr
PFA = zeros(length(thresholds),1);  % PFA para cada thr

ths_stats = zeros(4,length(thresholds)); % TPs, FNs, FPs, TNs = filas en ese orden ; Ths = columnas

for i=1:1:length(thresholds)
    ths_stats(1,i) = sum(X(i,:));                   % TPs para th i
    ths_stats(2,i) = length(X(i,:))-ths_stats(1,i); % FNs para th i
    ths_stats(3,i) = sum(Y(i,:));                   % FPs para th i
    ths_stats(4,i) = length(Y(i,:))-ths_stats(3,i); % TNs para th i
    
    PD(i) = (ths_stats(1,i))/length(vector_COI);
    PFA(i) = (ths_stats(3,i))/length(noise_samples);
    i
end

%% Gráficos
range = 200;
power_gain = rho*ARX/(4*pi*range.^2);
prx_theo = PTX*power_gain;              % Potencia teórica (power gain es atenuación)
theo_snr = prx_theo*Tmeas/(q_elect/RPD);
theo_snr_dB = 10*log10(theo_snr);

[Pd_teo,Pfa_teo] = rocsnr(theo_snr_dB,SignalType='NonFluctuatingNonCoherent');
% figure
hold on
semilogx(Pfa_teo,Pd_teo);grid on;xlabel("PFA");ylabel("PD");title("ROC");xlim([2.5e-7 1])   % ROC computada

% hold on
% semilogx(PFA,PD);grid on;xlabel("PFA");ylabel("PD");


