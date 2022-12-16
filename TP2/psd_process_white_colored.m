clear all;
close all;

%% White noise
x_white_ = load('process_data_white.mat');
x_colored_ = load('process_data_colored.mat');

white_process = x_white_.x;
colored_process = x_colored_.y;

%% 
NFFT=2048;
media = mean(white_process(:));
[pxx,freqv] = pwelch(white_process(:)-media,hanning(NFFT/2),0,NFFT);
figure
plot(freqv,pxx)
grid on
xlabel('Discrete Frequency [rad]')
ylabel('PSD [V^2/Hz]')


aaa = colored_process;
aaa2 = aaa(:);

NFFT=2048;
media2 = mean(aaa2);
[pxx,freqv] = pwelch(aaa2-media2,hanning(NFFT/2),0,NFFT,'oneside');
hold all
plot(freqv,pxx)
grid on
xlabel('Discrete Frequency [rad]')
ylabel('PSD [V^2/Hz]')

%% 
varianza = var(colored_process(:));
varianza2 = sum(pxx)*(2*pi/NFFT);