clear all; close all; clc

%columnas = samples
%filas = experimentos

signal = load('signal_exercise_1.mat');

%Señal a procesar
signal_to_process = signal.x;

%% Calculo de media y varianza
media = mean(signal_to_process(:))
varianza = var(signal_to_process(:))

%% Media como funcion del tiempo
size_x = size(signal_to_process);
Nexp = size_x(1)
Nsamp = size_x(2)

esperanza_definicion = sum(signal_to_process,1)/Nexp;

%Plotea la espereranza por definicion y le superpone la media
figure
plot(esperanza_definicion);
ylim([-0.3,0.3]);
hold all;
plot(1:length(esperanza_definicion), media*ones(size(esperanza_definicion)),'-r');
grid on;
legend('Esperanza','Media'); 
xlabel('Tiempo')

%% Varianza en funcion del tiempo
aux_mtx = (signal_to_process-media).^2;
varianza_definicion = sum(aux_mtx,1)/Nexp;    %preg
figure
plot(varianza_definicion);
ylim([0.4,0.6]);
hold all;
plot(1:length(varianza_definicion), varianza*ones(size(varianza_definicion)),'-r');
legend('Varianza/tiempo','Varianza'); xlabel('Tiempo')
%% Calculo autocorrelacion   
figure
for t1=1:250:1000
%     t1 = 1; % indice de tiempo a barrer
    term1 = signal_to_process(:,t1); % Devuelve la muestra t1 de todos los experimentos
    z = term1.*conj(signal_to_process);
    correlation_for_t1 = sum(z,1)/Nexp;
    plot(correlation_for_t1);
    hold on
end
legend('t = 1', 't = 250', 't = 500', 't = 750')