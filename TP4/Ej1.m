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
figure
plot(tline, real(s_t));
hold all
plot(tline, real(ch_out));

%% Receptor
wait_samples = ceil(Twait*fs);
detector_out_noiseless = ch_out.*conj(chirp_tx);  % Salida del detector (mixer) con delay
dsp_input_noiseless = detector_out_noiseless(1+wait_samples:end);   % Salida del detector sin delay
noise_power = q_elect/RPD * fs;     % Potencia del ruido

prx_measured = mean(abs(dsp_input_noiseless).^2);   % Potencia medida de señal recibida
prx_theo = PTX*power_gain;                          % Potencia teórica
theo_snr = prx_theo*Tmeas/(q_elect/RPD);
theo_snr_dB = 10*log10(theo_snr);

Ncells = ceil(fs*Tmeas);    % numero de bines de la FFT minimo
FFT_NOS = 2;                % Sobremuestreo de la FFT
NFFT = FFT_NOS*Ncells;
fvec = (0:NFFT-1)*(fs/NFFT);

fbeat = chirp_slope*2*real_range/c;     % Frecuencia del tono resultante
theo_max_pos = round(fbeat/fs*NFFT)+1;
cell_of_interest = round(fbeat*Tmeas);  % Celda de interés
total_cells = fs*Tmeas;
fft_dec_phase = mod(theo_max_pos-1,FFT_NOS);
fvec_dec = fvec(1+fft_dec_phase:FFT_NOS:end);

tline2 = tline(1+wait_samples:end);

Niters = 1;
range_est_v=zeros(Niters,1);
for niter=1:Niters
    niter
    
    noise = sqrt(noise_power/2).*(randn(size(ch_out))+1j.*randn(size(ch_out))); % Señal de ruido
    detector_out = detector_out_noiseless+noise;    % Salida del detector con ruido
    dsp_input = detector_out(1+wait_samples:end);   % Salida del detector con ruido sin delay

    % 
%     figure
%     plot(tline, real(detector_out_noiseless))
%     hold all
%     plot(tline, real(noise))
%     plot(tline, real(detector_out))
%     plot(tline2, real(dsp_input))

    y_mf = abs(fft(dsp_input, NFFT)).^2;
    
    if niter==1
        y_mf_accum = y_mf/Niters;
    else
        y_mf_accum = y_mf_accum+y_mf/Niters;
    end
    
    [~,max_pos] = max(y_mf);
    fbeat_meas = fvec(max_pos);
    range_est = (fbeat_meas/chirp_slope)*3e8/2;
    range_est_v(niter)= range_est;
    
    
    y_mf_dec = y_mf(1+fft_dec_phase:FFT_NOS:end);
    
    figure
    plot(fvec,y_mf,'-x');
    hold all
    plot(fvec_dec,y_mf_dec,'o');
    
end

figure
plot(y_mf_accum)

% Calculamos la SNR y dio bien

%%
range_error = range_est_v - real_range;
figure
hist(range_error, 100)