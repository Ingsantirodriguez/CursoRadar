clear all
close all
clc

%%
PTX = 1e-3*10^(15/10);  % Potencia transmitida [W]
chirp_bw = 2e9;         % Ancho de banda del chirp
max_range = 300;        % Rango máximo
c = 3e8;
Twait = 2*max_range/c;
q_elect = 1.6e-19;      % Carga del electrón, para calcular ruido
Tmeas = 2.5e-6;         % Tiempo de medición
Tmod = Tmeas+Twait;     % Tiempo de modulación
chirp_slope = chirp_bw/Tmod;    % Pendiente del chirp

% Parámetros del canal
range = 200;                %m
ARX=pi*(2.5e-2/2)^2;        % Apertura de 1in de diametro
rho = 0.1;                  % Reflectividad
lambda0 = 1550e-9;          % m
omega0 = 2*pi*c/lambda0;    % Frecuencia inicial
RPD = 0.7;                  % A/W Responsitividad del diodo

% Muestreo
NOS = 4;
fs=NOS*chirp_bw;          % Frec. de muestreo de Matlab
Ncells = ceil(fs*Tmeas);
fs=Ncells/Tmeas;        % Frecuencia de muestreo
Ts = 1/fs;              % Período de muestreo

%% Transmisor
Lsim = ceil(Tmod*fs);                       % Largo de la sim.
tline = Ts.*(0:Lsim-1)';                    % Linea de tiempo de modulación
insta_freq = chirp_slope.*tline;            % Frecuencia instantánea
insta_phase = 2*pi*cumsum(insta_freq).*Ts;  % Fase instantánea
chirp_tx = exp(1j*insta_phase);             % Chirp unitario

s_t = sqrt(PTX).*chirp_tx;                  % Señal transmitida

%% Canal
tau = 2*range/c;              % Delay del canal
delay_samples = round(tau*fs);
real_tau = delay_samples*Ts;
real_range = real_tau*c/2;    % Rango real

power_gain = rho*ARX/(4*pi*range.^2);       % Potencia
atten = sqrt(power_gain);                   % Atenuación
delta_phase = 2*pi*c/(lambda0*real_tau);    % Cambio de fase

ch_out = atten.*[zeros(delay_samples,1); s_t(1:end-delay_samples)].*exp(-1j*delta_phase);
    % Salida del canal (atten, fase y delay)
% figure
% plot(tline, real(s_t));
% hold all
plot(tline, real(ch_out));

%% Receptor
wait_samples = ceil(Twait*fs);
mixer_out = ch_out.*conj(chirp_tx);  % Salida del detector (mixer) con delay

noise_power = q_elect/RPD * fs;     % Potencia del ruido



