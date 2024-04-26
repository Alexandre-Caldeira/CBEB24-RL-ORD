clearvars; 
close all;clc;

snr = 5;
max_length = 400;
signal_freq_bins= [82    84    86    88    90    92    94    96];
noise_freq_bins=round(signal_freq_bins.*pi)+1; 
res = checkForHarmonics(signal_freq_bins, noise_freq_bins);
if res~= false
    error(['Signal and Noise frequencies contain harmonics in noise_freq=',...
        num2str(res(1))])
end

% fixed sampling rate...
FS = 1000; 
NFFT = FS; % rever isso aqui ! motivo/causas/impactos comp. 
Njanelas = max_length;
NpontosTotal = NFFT*Njanelas;    % Numero total de pontos de cada sinal
tempo = (0:NpontosTotal-1)/FS;   % Vetor de tempo utilizado para gerar o sinal    

% Gerando e observando sinais ruidosos no tempo e frequencia
sigma_n = 2/NFFT;
snr = 10^(snr/10);

% Constante a ser multiplicada ao ruido e ao sinal para configurar a relação sinal ruido desejada
SNRs = sqrt(4*sigma_n*snr/NFFT);
SNRn = sqrt(sigma_n);

% cria matriz para receber a simulacao
signals = nan(numel(signal_freq_bins), ...
                NFFT*max_length);   

% cria os sinais: (size = 8x1000xmax_length)
for signal_idx = 1:numel(signal_freq_bins)
    SFREQ = signal_freq_bins(signal_idx)-1;
    sinal_puro = SNRs*sin(2*pi*SFREQ*tempo+rand()*2*pi);
    sinal_ruidoso = sinal_puro;
    signals(signal_idx ,:) = sinal_puro;
end

ruido = randn(1,NpontosTotal);  % Gera um ruido gaussiano, teoricamente de variancia unitaria e media nula
ruido = ruido-mean(ruido);      % Força a media nula
ruido = ruido/std(ruido)*SNRn;  % Força a variancia desejada para o sinal
signals = sum(signals,1)+ruido;
signals = signals./std(signals);
signals = reshape(signals,NFFT, Njanelas);

% calcula transformada e estados contínuos:
SIGNALS = fft(signals);
SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid

MSC = nan(numel(f),max_length);
for idx_episodio = 1:max_length 

    % para cada frequencia em f
    for idx_f = 1:1:numel(f)
        
        X_atual = SIGNALS(idx_f,1:idx_episodio);

        num_msc = abs(sum(X_atual))^2;
        den_msc = idx_episodio*sum(abs(X_atual).^2);
        MSC(idx_f,idx_episodio) = num_msc/den_msc;
   
    end
end

% cria matriz dos estados continuos:
all_freqs = [signal_freq_bins,noise_freq_bins];
abs_fft = rescale( abs(SIGNALS(all_freqs,:)),   0,1 );  
ang_fft = rescale( angle(SIGNALS(all_freqs,:)),-1,1 );
re_fft  = rescale( real(SIGNALS(all_freqs,:)), -1,1 );
im_fft  = rescale( imag(SIGNALS(all_freqs,:)), -1,1 );  
msc = MSC(all_freqs,:);
epp = repmat(100.*[1:max_length]./max_length,16,1);

c_states = nan(numel(signal_freq_bins)+numel(noise_freq_bins), ... % freqs
                6,... % states
                max_length); % windows

c_states(:,1,:) = abs_fft;
c_states(:,2,:) = ang_fft;
c_states(:,3,:) = re_fft;
c_states(:,4,:) = im_fft;
c_states(:,5,:) = msc;
c_states(:,6,:) = epp;

% discretiza a matriz de estados
resolucao = 10; % 10 niveis discretos para cada estado
states =  nan(size(c_states));
min_vals = [0,-1,-1,-1,0,0]; % min for abs, ang, re, im, msc, ep%
max_vals = [1,1,1,1,1,100];  % max for abs, ang, re, im, msc, ep%

for idx_freq = 1:numel(all_freqs)
    for idx_state = 1:6
        values = c_states(idx_freq,idx_state,:);
        min_val = min_vals(idx_state);
        max_val = max_vals(idx_state);

        states(idx_freq,idx_state,:) = discretize_val( ...
                                        values, min_val, max_val, resolucao);
    end
end
%%
% tic()
% signal_freq_bins =  [82    84    86    88    90    92    94    96]
% noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5
% snr = 15;
% max_length = 400;
% resolution = 10
% size(rlord_gen_states(signal_freq_bins, ...
%                         noise_freq_bins, ...
%                         snr, ...
%                         max_length, ...
%                         resolution))
% 
% toc()

%% validate 
f = FS/2*linspace(0,1,NFFT/2+1)';
fsins = abs(SIGNALS);
% fsins = fsins(1:floor(end/2)+1,:);
figure(1)
stem(f,fsins(:,3))
hold on 
stem(f(signal_freq_bins),fsins(signal_freq_bins,3),'r')
hold off
disp(c_states(1,6,1:20))
figure(2)
stem(f(all_freqs),c_states(:,1,1))