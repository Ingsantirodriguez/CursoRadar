clear all
close all
clc

%%
PTX = 1e-3*10^(15/10);  % Potencia transmitida [W], Aparecen muchos outliers si se baja mucho
chirp_bw = 2e9;         % Ancho de banda del chirp
max_range = 300;        % Rango máximo
c = 3e8;                % Velocidad de la luz
Twait = 2*max_range/c;  % Tiempo de espera a un target r = r_max
q_elect = 1.6e-19;      % Carga del electrón, para calcular ruido
Tmeas = 2.5e-6;         % Tiempo de medición
Tmod = Tmeas+Twait;     % Tiempo de modulación
chirp_slope = chirp_bw/Tmod;    % Pendiente del chirp

% Parámetros del canal
range = 200;                % [m]
ARX=pi*(2.5e-2/2)^2;        % Apertura de 1in de diametro
rho = 0.1;                  % Reflectividad
lambda0 = 1550e-9;          % Ancho de banda inicial [m]
omega0 = 2*pi*c/lambda0;    % Frecuencia inicial
RPD = 0.7;                  % A/W Responsitividad del diodo (para calcular ruido)

% Muestreo
NOS = 4;                    % Sobremuestreo para gráficos con más detalle
fs=NOS*chirp_bw;            % Frec. de muestreo de Matlab
Ncells = ceil(fs*Tmeas);    % Número de celdas de la FFT (para threshold y CFAR)
fs=Ncells/Tmeas;        % Frecuencia de muestreo como numero entero  
Ts = 1/fs;              % Período de muestreo

%% Transmisor
Lsim = ceil(Tmod*fs);                       % Largo de la sim. en muestras
tline = Ts.*(0:Lsim-1)';                    % Linea de tiempo de modulación (Lsim-1 es para que tengan la misma longitud los 2)
insta_freq = chirp_slope.*tline;            % Frecuencia instantánea
insta_phase = 2*pi*cumsum(insta_freq).*Ts;  % Fase instantánea !Ver porqué no usa exponencial
chirp_tx = exp(1j*insta_phase);             % Chirp unitario

s_t = sqrt(PTX).*chirp_tx;                  % Señal transmitida

%% Canal
tau = 2*range/c;              % Delay del canal
delay_samples = round(tau*fs);
real_tau = delay_samples*Ts;
real_range = real_tau*c/2;    % Rango real

power_gain = rho*ARX/(4*pi*range.^2);       % Potencia (ganancia menor a 1 = atenuación)
atten = sqrt(power_gain);                   % Atenuación
delta_phase = 2*pi*c/(lambda0*real_tau);    % Cambio de fase

ch_out = atten.*[zeros(delay_samples,1); s_t(1:end-delay_samples)].*exp(-1j*delta_phase);
    % Salida del canal (atten, delay y fase)
figure
plot(tline, real(s_t));
hold all
plot(tline, real(ch_out));

%% Receptor
wait_samples = ceil(Twait*fs);  % Twait en muestras
detector_out_noiseless = conj(ch_out.*conj(chirp_tx));  % Salida del detector (mixer), se conjuga todo porque sino da frecuencias negativas
dsp_input_noiseless = detector_out_noiseless(1+wait_samples:end);   % Salida del detector sin Twait

noise_power = q_elect/RPD * fs;     % Potencia del ruido

prx_measured = mean(abs(dsp_input_noiseless).^2);   % Potencia medida de señal recibida
prx_theo = PTX*power_gain;                          % Potencia teórica (power gain es atenuación)
theo_snr = prx_theo*Tmeas/(q_elect/RPD)
theo_snr_dB = 10*log10(theo_snr)    % Aparecen outliers de 14dB para abajo

Ncells = ceil(fs*Tmeas);    % numero de bines/celdas de la FFT minimo | Cuántas celdas me entrarían dentro de la FFT en el ancho de banda que estoy planteando
FFT_NOS = 16;                % Sobremuestreo de la FFT, para poder medir precisión, usar valor mayor a 8 al medir precisión
NFFT = FFT_NOS*Ncells;      % Tamaño de la transformada
fvec = (0:NFFT-1)*(fs/NFFT);    % Vector de frecuencia, la cantidad de frecuencias es el tamaño de la transformada, el paso de frecuencia es fs/NFFT
                                % fs es el valor máx de la FFT, entonces el
                                % paso es fs/NFFT
fbeat = chirp_slope*2*real_range/c;     % Frecuencia teórica del tono resultante
cell_of_interest = round(fbeat*Tmeas);  % Celda de interés = fbeat/(1/Tmeas)
% total_cells = fs*Tmeas;   % El ancho de cada celda es 1/Tmeas

% Como el numero total de celdas es fs*Tmeas, y la frecuencia mayor es fs
% entonces tiene sentido que la COI sea fbeat*Tmeas.

% Decimado/División en celdas
theo_max_pos = round(fbeat/fs*NFFT)+1;  % Posición teórica del máximo (fbeat) en la transformada
fft_dec_phase = mod(theo_max_pos-1,FFT_NOS);    % Posición de la celda donde está el máximo
fvec_dec = fvec(1+fft_dec_phase:FFT_NOS:end);   % Vector de frecuencia decimado, con máximo en medio de la celda

% fft_dec_phase Es un número tal que si empiezo a decimar desde ese punto
% con pasos FFT_NOS entonces siempre voy a tener el máximo en el centro de
% la celda.

tline2 = tline(1+wait_samples:end);

Niters = 1;
range_est_v=zeros(Niters,1);
for niter=1:Niters
    niter;
    
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

    y_mf = abs(fft(dsp_input, NFFT)).^2;    % Salida del match filter (fft) potencia
    
    if niter==1
        y_mf_accum = y_mf/Niters;
    else
        y_mf_accum = y_mf_accum+y_mf/Niters;
    end
    
    [~,max_pos] = max(y_mf);    % Valor máximo de la fft
    fbeat_meas = fvec(max_pos); % Ubicación del valor máximo en el vector de frecuencia (fbeat medido)
    range_est = (fbeat_meas/chirp_slope)*c/2;   % Rango estimado con la fbeat medida
    range_est_v(niter)= range_est;  % Vector de rangos estimados (cálculo de precisión)
    
    
    y_mf_dec = y_mf(1+fft_dec_phase:FFT_NOS:end);   % Vector y_mf decimado
    
    figure
    plot(fvec,y_mf,'-x');
    hold all
    plot(fvec_dec,y_mf_dec,'o');
    
end

figure
plot(fvec/1e6, y_mf_accum)

[max_val,max_pos] = max(y_mf);  % Valor e índice del pico de la señal
mean_noise = mean(y_mf_accum(max_pos+1000:end));    % Media del ruido
snr_comp = (max_val)/mean_noise % Cálculo de SNR
snr_comp_dB = 10*log10(snr_comp)    % SNR logarítmico

% Calculamos la SNR y dio bien

%%
range_error = range_est_v - real_range;
figure
hist(range_error, 100)