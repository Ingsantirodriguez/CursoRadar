clear all
% close all
clc

%%
PTX = 1e-3*10^(15/10);  % Potencia transmitida [W] 15 dBm
chirp_bw = 2e9;         % Ancho de banda del chirp
max_range = 350;        % Rango máximo
c = 3e8;
Twait = 2*max_range/c;
q_elect = 1.6e-19;      % Carga del electrón, para calcular ruido
Tmeas = 3.5e-6;         % Tiempo de medición
Tmod = Tmeas+Twait;     % Tiempo de modulación
chirp_slope = chirp_bw/Tmod;    % Pendiente del chirp

% Parámetros del canal
range = 350;            %m
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

%% Transmisor
Lsim = ceil(Tmod*fs);                       % Largo de la sim.
tline = Ts.*(0:Lsim-1)';                    % Linea de tiempo de modulación
insta_freq = chirp_slope.*tline;            % Frecuencia instantánea
insta_phase = 2*pi*cumsum(insta_freq).*Ts;  % Fase instantánea
chirp_tx = exp(1j*insta_phase);             % Chirp unitario

s_t = sqrt(PTX).*chirp_tx;                  % Señal transmitida
j = 1;
start = 50;
step = 25;
theo_snr = zeros(length(start:step:max_range),1);
theo_snr_dB = zeros(length(start:step:max_range),1);
snr_comp = zeros(length(start:step:max_range),1);
snr_comp_dB = zeros(length(start:step:max_range),1);
for range=start:step:max_range
    j
    %% Canal
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
    mixer_v = mixer(1+wait_samples:end);
    ch_out_v = ch_out(1+wait_samples:end);

    max2 = max(abs(mixer_v).^2);
    [max1,I1] = max(abs(ch_out_v).^2);

    % Matched Filter FFT

    Ncells = ceil(fs*Tmeas);    % numero de bines de la FFT minimo
    FFT_NOS = 32;               % Sobremuestreo de la FFT
    NFFT = FFT_NOS*Ncells;
    fvec = (0:NFFT-1)*(fs/NFFT);            % Vector de frecuencia
    max_range_FFT = fs*c/(chirp_slope*2);   % Máximo rango que puede ser representado por la FFT
    rvec = (0:NFFT-1)*(max_range_FFT/NFFT); % Vector de rango

    % Experimentos
    noise_power = q_elect/RPD*fs;     % Potencia del ruido
    Nexp = 500;
    f_vals = zeros(Nexp,1);
    r_vals_ol = zeros(Nexp,1);
    r_vals_sol = zeros(Nexp,1);

    for i=1:Nexp

        noise = sqrt(noise_power/2).*(randn(size(ch_out))+1j.*randn(size(ch_out))); % Señal de ruido
        noise_v = noise(1+wait_samples:end);
        fe_output = mixer_v + noise_v; %Salida del detector con ruido

        y_mf = abs(fft(fe_output, NFFT)).^2;

        if i==1
            y_mf_accum = y_mf/Nexp;
        else
            y_mf_accum = y_mf_accum+y_mf/Nexp;
        end

%         i

    end
%     hold on
%     plot(y_mf_accum);grid on

    % SNR Teórica
    prx_theo = PTX*power_gain;                  % Valor de la potencia recibida (teórico)
    theo_snr(j) = prx_theo*Tmeas/(q_elect/RPD); % SNR Teórica
    theo_snr_dB(j) = 10*log10(theo_snr(j));     % SNR Teórica [dB]

    % SNR Computada
    [max_val,max_pos] = max(y_mf_accum);  % Valor e índice del pico de la señal
    mean_noise = mean(y_mf_accum(max_pos+1000:end));    % Media del ruido
    snr_comp(j) = (max_val-mean_noise)/mean_noise; % Cálculo de SNR
    snr_comp_dB(j) = 10*log10(snr_comp(j));    % SNR logarítmico

    j = j+1;
end

% Gráficos

rline_std = (start:step:max_range);
% figure
hold on
semilogy(rline_std,theo_snr_dB);grid on;
xlabel('Rango [m]');
% semilogy(rline_std,snr_comp_dB);
% fbeat = chirp_slope*2*real_range/c;