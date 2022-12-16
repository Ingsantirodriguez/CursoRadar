clear all; close all; clc

%columnas = samples
%filas = experimentos


signal = load('signal_exercise_1.mat');

%Señal en el eje x
signal_to_process = signal.x;
% plot(signal_to_process)

%Histogram
% figure
% hist(signal_to_process(:),100); grid on;

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
plot(esperanza_definicion); grid on; 
xlabel('tiempo')
ylim([-0.3,0.3])
hold all
plot(1:length(esperanza_definicion), media*ones(size(esperanza_definicion)),'-r')
legend('Esperanza','Media')

%% Varianza en funcion del tiempo
aux_mtx = (signal_to_process-media).^2;
varianza_definicion = sum(aux_mtx,1)/Nexp;    %preg
figure
plot(varianza_definicion); grid on;
ylim([0.4,0.6]);
hold all;
plot(1:length(varianza_definicion), varianza*ones(size(varianza_definicion)),'-r')
legend('Varianza/tiempo', 'Varianza','Location','northeast'); xlabel('tiempo');

%% Calculo autocorrelacion
% % % Rx(tau) = E{ x(t) x*(t-tau) }
% % % defino g(tau) = x(t) x*(t-tau)
% g_tau = zeros(2*length(signal_to_process(:,1))-1,Nexp ); % Guardo el de g(tau)
% 
% for nexp=1:1000
%     % 
%     %g_tau(:,nexp)= 1/Lsamp*conv(x_white(:,nexp),conj(x_white(end:-1:1,nexp)));
%     g_tau(:,nexp) = xcorr(signal_to_process(:,nexp) - media,'unbiased'); % equivalente a lo anterior
%     % Biased es para que normalice con la longitud del vector de entrada
% end
% corr_result = sum(g_tau,2)/Nsamp;
% 
% % figure
% % plot(g_tau(:,1:1000)); grid on; 
% % % ylim([-1,1])
% 
% figure
% plot(corr_result); grid on;title('Autocorrelacion')
% ylim([-0.3,0.6;])
% max(corr_result)



t1=1; % indice de tiempo a barrer
term1=signal_to_process(:,t1);
z = term1.*conj(signal_to_process);
correlation_for_t1 = sum(z,1)/Nexp;
figure
plot(correlation_for_t1); grid on; %ylim([-0.5;0.5])
% title('Autocorrelacion para t1=1');
hold all

t1=100; % indice de tiempo a barrer
term1=signal_to_process(:,t1);
z = term1.*conj(signal_to_process);
correlation_for_t1 = sum(z,1)/Nexp;
plot(correlation_for_t1); grid on; %ylim([-0.5;0.5])
% title('Autocorrelacion para t1=1');
hold all

t1=250; % indice de tiempo a barrer
term1=signal_to_process(:,t1);
z = term1.*conj(signal_to_process);
correlation_for_t1 = sum(z,1)/Nexp;
plot(correlation_for_t1); grid on; %ylim([-0.5;0.5])
% title('Autocorrelacion para t1=1');
hold all

t1=500; % indice de tiempo a barrer
term1=signal_to_process(:,t1);
z = term1.*conj(signal_to_process);
correlation_for_t1 = sum(z,1)/Nexp;
plot(correlation_for_t1); grid on; %ylim([-0.5;0.5])
% title('Autocorrelacion para t1=1');
hold all

t1=750; % indice de tiempo a barrer
term1=signal_to_process(:,t1);
z = term1.*conj(signal_to_process);
correlation_for_t1 = sum(z,1)/Nexp;
plot(correlation_for_t1); grid on; %ylim([-0.5;0.5])
% title('Autocorrelacion para t1=1');
legend('t=1', 't=100', 't=250', 't=500', 't=750')

%OPCION 2
% % Rx(tau) = E{ x(t) x*(t-tau) }
%esperanza_definicion = sum(signal_to_process,1)/Nexp;

% z=0;
% for t1=1:1000
%     xt = signal_to_process(:,t1);
%     
%     %xc = signal_to_process(:,t1-z);  
%     g_aux = xt.*conj(signal_to_process);
%     z=z+1;
% end
% 
% 
% esperanzaE = sum(g_aux,1)/Nexp;
% plot(esperanzaE, 'r');

%OPCION 3
% figure
% for t1=1:1000
% %     t1=1; % indice de tiempo a barrer
%     term1=signal_to_process(:,t1);
%     z = term1.*conj(signal_to_process);
%     correlation_for_t1 = sum(z,1)/10000;
% %     plot(correlation_for_t1);
% %     hold all
% end

%% 

% plot(signal_to_process(1,:)); grid on;
% hold on;
% plot(signal_to_process(50,:));
% hold on;
% plot(signal_to_process(10000,:)); legend('e=1', 'e=50', 'e=10k')
