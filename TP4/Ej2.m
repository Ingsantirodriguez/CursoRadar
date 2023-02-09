clear all
% close all
clc

%%
PTX = 1e-3*10^(25/10);  % Potencia transmitida [W] 15 dBm
chirp_bw = 2e9;         % Ancho de banda del chirp
max_range = 350;        % Rango máximo
c = 3e8;
Twait = 2*max_range/c;
q_elect = 1.6e-19;      % Carga del electrón, para calcular ruido
Tmeas = 2.5e-6;         % Tiempo de medición
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
for range=start:step:max_range
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
        [~,I] = max(y_mf);
        f_vals(i) = I*(fs/NFFT);  % Fbeat simulado
        r_vals_ol(i) = f_vals(i)*c/(chirp_slope*2);  % Rango medido

        r_vals_sol = remove_outliers(r_vals_ol, range, 1);

        i;

    end

    % frec_std = std(f_vals);
    % range_std = std(r_vals);
    range_err = r_vals_ol-real_range;
    range_err_sol = r_vals_sol-real_range;
    range_std(j) = std(range_err);
    range_std_sol(j) = std(range_err_sol);

    j = j+1
end

% figure
% histogram(range_err,100);grid on;title('Con Outliers');
% 
% figure
% histogram(range_err_sol,100);grid on;title('Sin Outliers');
rline_std = (start:step:max_range);
hold on
semilogy(rline_std,range_std_sol);grid on;title('Precisión sin outliers')
ylabel('Desvío estándar');xlabel('Rango [m]');
fbeat = chirp_slope*2*real_range/c;
% cell_of_interest = round(fbeat*Tmeas);  % Celda de interés
% total_cells = fs*Tmeas;
% fft_dec_phase = mod(theo_max_pos-1,FFT_NOS);
% fvec_dec = fvec(1+fft_dec_phase:FFT_NOS:end);