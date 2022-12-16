clear all;
close all;

%% White noise
x_white_ = load('/home/abrila/Documents/Beca/TP2/process_data_white.m');
%x_colored_ = load('process_data_colored.mat');

noise_to_process = x_white_.x;
%noise_to_process = x_colored_.y;

figure 
hist(noise_to_process(:), 100);
%%
media = mean(noise_to_process(:));
varianza = var(noise_to_process(:));

%% Esperanza por definicion
size_x = size(noise_to_process);
Nexp = size_x(2);
Nsamp = size_x(1);
esperanza_definicion = sum(noise_to_process,2)/Nexp;

figure
plot(esperanza_definicion)
% ylim([0,2])
hold all
plot(1:length(esperanza_definicion), media*ones(size(esperanza_definicion)),'-r')

%% Varianza por definicion
aux_mtx = abs(noise_to_process-media).^2;
varianza_definicion = sum(aux_mtx,2)/Nexp;
figure
plot(varianza_definicion)
%ylim([0,2])
hold all
plot(1:length(varianza_definicion), varianza*ones(size(varianza_definicion)),'-r')

%% Auto-correlacion
%% Calculo autocorrelacion
% Rx(tau) = E{ x(t) x*(t-tau) }
% defino g(tau) = x(t) x*(t-tau)
g_tau = zeros(2*length(noise_to_process(:,1))-1,Nexp); % Guardo el de g(tau)

for nexp=1:Nexp
    % 
    %g_tau(:,nexp)= 1/Lsamp*conv(x_white(:,nexp),conj(x_white(end:-1:1,nexp)));
    g_tau(:,nexp)= xcorr(noise_to_process(:,nexp)-media,'unbiased'); % equivalente a lo anterior
    % Biased es para que normalice con la longitud del vector de entrada
end
corr_result = 1/Nexp * sum(g_tau,2);

figure
plot(g_tau(:,1:10))

figure
plot(corr_result,'-o')

max(corr_result)
