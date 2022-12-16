clear all; close all; clc

signal = load('signal_exercise_2.mat');

%Señal en el eje x
signal_to_process = signal.x;
 %plot(signal_to_process(1:3,:))

%Histogram
figure
hist(signal_to_process(:),100); grid on;

%% Calculo de media y varianza
media = mean(signal_to_process(:))
varianza = var(signal_to_process(:))

%% Media como funcion del tiempo
size_x = size(signal_to_process);
Nexp = size_x(1)  %10k
Nsamp = size_x(2)   %1k
esperanza_definicion = sum(signal_to_process,1)/Nexp;

%Plotea la espereranza por definicion y le superpone la media
figure
plot(esperanza_definicion); grid on
% ylim([-0.5,0.5])
hold all
plot(1:length(esperanza_definicion), media*ones(size(esperanza_definicion)),'-r')
legend('Esperanza', 'Media'); xlabel('Tiempo')
%% Varianza en funcion del tiempo
aux_mtx = abs(signal_to_process-media).^2;
varianza_definicion = sum(aux_mtx,1)/Nexp;
figure
plot(varianza_definicion);grid on;
% ylim([0,0.6])
hold all
plot(1:length(varianza_definicion), varianza*ones(size(varianza_definicion)),'-r')
legend('Varianza/tiempo', 'Varianza'); xlabel('Tiempo')
%% Calculo autocorrelacion
% Rx(tau) = E{ x(t) x*(t-tau) }
% defino g(tau) = x(t) x*(t-tau)
g_tau = zeros(2*length(signal_to_process(:,1))-1,Nexp); % Guardo el de g(tau)

for nexp=1:1000
    % 
    %g_tau(:,nexp)= 1/Lsamp*conv(x_white(:,nexp),conj(x_white(end:-1:1,nexp)));
    g_tau(:,nexp)= xcorr(signal_to_process(:,nexp)-media,'unbiased'); % equivalente a lo anterior
    % Biased es para que normalice con la longitud del vector de entrada
end
corr_result = 1/Nexp * sum(g_tau,2);

figure
plot(g_tau(:,1:10))

figure
plot(corr_result);grid on;title('Autocorrelacion');

max(corr_result)

% %% 
% t1 = 1; % indice de tiempo a barrer
% term1 = signal_to_process(:,t1); % Devuelve la muestra en t1 de todos los experimentos
% conjugado = conj(signal_to_process);
% z = term1.*conjugado;
% size(term1) %10k
% size(z); %2 dimensiones
% correlation_for_t1 = sum(z,1)/nexp;
% plot(correlation_for_t1);
% % hold all
% figure
% autoC = autocorr(term1);
% plot(autoC)

%%
figure
plot(signal_to_process(1,:))
hold on
plot(signal_to_process(2,:))
hold on
plot(signal_to_process(3,:))
hold on
plot(signal_to_process(4,:))
hold on
plot(signal_to_process(5,:))



