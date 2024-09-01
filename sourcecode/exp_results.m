%% CLEAN ENV
clearvars;  
% close all;
clc;
dbpath = 'C:\PPGEE\SBEB_CBA_24\ASSR - Coleta OFFLINE';
addpath(genpath(dbpath))
%% DEFINE AGENT
for cont_vol = 1:11
resolution = 10;
sssize = [resolution,resolution,resolution,resolution,2];
% qtable = zeros(sssize);
% qtable = load('longertests_ws_snr_-6_numep_500_dt_23-Apr-2024-14-26-09','qtable');
% filename = 'efeito-aumento-alpha_ws_snr_-6_numep_500_dt_23-Apr-2024-14-26-09';
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-2_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-11_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-15_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-12_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-3_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-8_qtable.mat')
% load('C:\PPGEE\EEE935 Aprendizado por Reforço\RLORD\treino_snr_-10_qtable.mat')
load('log10states_ws_snr_-8_numep_500_dt_23-Apr-2024-14-26-09.mat')

% qtable = load([filename,'.mat'],'qtable');
% qtable = qtable.qtable;

% alpha = 5e-1;
% gamma = 0.3;
epsilon = @(episode) 0.2*exp(-0.5*log10(episode));
% plot(0.2*exp(-0.5*log10(1:5*1e3)))

test_q = qtable;
test_alpha =  1;
test_gamma = 0.95;

%% GET DATA
% Vintensidade = {'70';'60';'50';'40';'30'};

% wtf happens at
% cont_vol = 9;
% cont_int = 3;

% vol 3, int 4 ficou bom na log10states_ws_snr_-8_numep_500_dt_23-Apr-2024-14-26-09

% cont_vol = 3
% cont_int = 1:5 %-> 4 e 2 -> 2

% cont_vol = 1
cont_int = 2;

signal_freq_bins =  [82 82   84    86    88    90    92    94    96];
noise_freq_bins = round(signal_freq_bins.*exp(1)*sqrt(2))+12;
noise_freq_bins = [noise_freq_bins,round(signal_freq_bins.*exp(1)*sqrt(2))-30];
% noise_freq_bins = [noise_freq_bins,round(signal_freq_bins.*exp(1)*sqrt(2))-100]
% noise_freq_bins = [noise_freq_bins,round(signal_freq_bins.*exp(1)*sqrt(2))+100]
% noise_freq_bins = [noise_freq_bins,round(signal_freq_bins.*exp(1)*sqrt(2))+35]
% noise_freq_bins = [noise_freq_bins,round(signal_freq_bins.*exp(1)*sqrt(2))-11,round(signal_freq_bins.*exp(1)*sqrt(2))+10,signal_freq_bins-15]

noise_freq_bins = noise_freq_bins(4);
signal_freq_bins = signal_freq_bins(4);

% noise_freq_bins = [1:500];
% noise_freq_bins([1:125,signal_freq_bins,signal_freq_bins*2, ...
    % signal_freq_bins*3,signal_freq_bins*4,signal_freq_bins*5, ...
    % 60*2,60*3,60*4,60*5, 480:500])=[];
% noise_freq_bins() = [];
% noise_freq_bins(signal_freq_bins*3) = [];
% noise_freq_bins(signal_freq_bins*4) = [];
% noise_freq_bins(signal_freq_bins*4) = [];


all_freq_bins = [signal_freq_bins,noise_freq_bins];

% size = 16 x 6 x max_length
all_states= db_to_dstates(signal_freq_bins, ...
                                noise_freq_bins, ...
                                resolution, ...
                                cont_vol, ...
                                cont_int);
% all_states(isnan(all_states))=1;

max_episode_length = size(all_states,3);
max_length = max_episode_length;
max_length_teste = max_length;

%% TESTE in EXP

tdr = 0;
tfp = 0;

resultados_teste =  nan(max_episode_length-1,3);
rewards_vec_teste  = nan(1,max_episode_length);
mean_tdr_vec_teste = nan(1,max_episode_length);
mean_tfp_vec_teste = nan(1,max_episode_length);


% initialize episode parameters
terminal_state = false;
total_reward_for_this_episode = zeros(1,max_length_teste-1);
tdr_hist= zeros(1,max_length_teste-1);
tfp_hist= zeros(1,max_length_teste-1);
is_freq_undecided = ones(1,numel(all_freq_bins));



for current_window = 1:max_length_teste-1  
    % for count = 1:10

    % take action for each freq
    dr = 0;
    fp = 0;
    decisions= zeros(numel(all_freq_bins),5);
    
    for idx_freq = 1:numel(all_freq_bins)
        if is_freq_undecided(idx_freq)

            current_freq = all_freq_bins(idx_freq);

            % size = 1 x 6
            current_states = all_states(idx_freq,:,current_window);
            next_states = all_states(idx_freq,:,current_window+1);

            % size = 1x2
            sa1 = sub2ind(sssize, ...
                current_states(1),current_states(2),...
                current_states(3),current_states(4),...
                1);
            sa2 = sub2ind(sssize, ...
                current_states(1),current_states(2),...
                current_states(3),current_states(4),...
                2);

            current_q_sa = [test_q(sa1),test_q(sa2)];
             
            nexts_a1 = sub2ind(sssize, ...
                next_states(1),next_states(2),...
                next_states(3),next_states(4),...
                1);

            nexts_a2 = sub2ind(sssize, ...
                next_states(1),next_states(2),...
                next_states(3),next_states(4),...
                2);

            next_q_sa = reshape(test_q([nexts_a1,nexts_a2]),1,2);
            [old_q, current_action] = max(current_q_sa);  
            max_q = max(next_q_sa(current_action));

            decisions(idx_freq,1:3) = [old_q,max_q, current_action];
            decisions(idx_freq,4:5) = current_q_sa;

            % if current_freq not in noise_freq_bins, should detect
            should_detect = isempty(find(noise_freq_bins==current_freq));
            if current_action==2 || current_window==max_length_teste

                % calculate FP and DR 
                if should_detect
                    dr = dr+1;

                else
                    fp = fp+1;

                end

                
            end

        end

    end

    tdr = 100*dr/numel(signal_freq_bins);
    tfp = 100*fp/numel(noise_freq_bins);
    el = 100*current_window/max_episode_length;
    reward = +((tfp)^2)/(-100) ...
             +((tdr)^2)/(100)-el/1000;   
    
    
    test_reward = +((tfp)^2)/(-100) ...
             +((tdr)^2)/(100)-el/1000;
    
    
    for idx_freq = 1:numel(all_freq_bins)
    
        % if current_freq not in noise_freq_bins, should detect
        should_detect = isempty(find(noise_freq_bins==current_freq));
    
        % if ~should_detect %is_freq_undecided(idx_freq)
    
        % decisions(idx_freq,1:3) = [old_q,max_q, current_action];
        % decisions(idx_freq,4:5) = current_q_sa;
        
        old_q = decisions(idx_freq,1);
        max_q = decisions(idx_freq,2);
        current_action = decisions(idx_freq,3);
        current_states = all_states(idx_freq,:,current_window);
    
        new_q =  old_q + test_alpha * (test_reward + test_gamma * max_q - old_q);
    
        selected_sa = sub2ind(sssize, ...
            current_states(1),current_states(2),...
            current_states(3),current_states(4),...
            current_action);
    
        test_q(selected_sa) = new_q;
        % end
    
    end
    % end
    
    
    % accumulate reward
    total_reward_for_this_episode(current_window) = reward;
    tdr_hist(current_window) = tdr;
    tfp_hist(current_window) = tfp;


    % accumulate reward
    total_reward_for_this_episode(current_window) = reward;
    tdr_hist(current_window) = tdr;
    tfp_hist(current_window) = tfp;

end

% mean_tdr_vec_teste(idx_window_teste) = mean(tdr_hist);
% mean_tfp_vec_teste(idx_window_teste) = mean(tfp_hist);
% rewards_vec_teste(idx_window_teste) = mean(total_reward_for_this_episode);
% resultados_teste(idx_window,3)=mean(mean_tdr_vec_teste);
% resultados_teste(idx_window,2)=mean(mean_tfp_vec_teste);
% resultados_teste(idx_window,1)=mean(rewards_vec_teste);

mean_tdr_vec_teste= tdr_hist;
mean_tfp_vec_teste = tfp_hist;
rewards_vec_teste = total_reward_for_this_episode;

resultados_teste(:,3)=mean_tdr_vec_teste;
resultados_teste(:,2)=mean_tfp_vec_teste;
resultados_teste(:,1)=rewards_vec_teste;

figure(cont_int+cont_vol)
subplot(2,2,[1 3])
plot(1:max_length-1,resultados_teste(:,1),'.', ...
    'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
title(['voluntario = ',num2str(cont_vol),'- average reward during test episodes ', num2str(mean(resultados_teste(:,1),'all'))])
hold on
plot(1:max_length-1,movmean(resultados_teste(:,1),[10 0],'omitnan'), 'linewidth',3)
grid on
xlabel('training epoch')
hold off


% figure(3)
subplot(2,2,2)
plot(1:max_length-1,resultados_teste(:,3),'.', ...
    'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
hold on
title(['voluntario = ',num2str(cont_vol),'- mean tdr during test episodes ', num2str(mean(resultados_teste(:,3),'all'))])
plot(1:max_length-1,movmean(resultados_teste(:,3),[10 0],'omitnan'), 'linewidth',3)
grid on
xlabel('training epoch')
hold off

% figure(2)
subplot(2,2,4)
plot(1:max_length-1,resultados_teste(:,2),'.', ...
    'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
hold on
title(['voluntario = ',num2str(cont_vol),'- mean tfp during test episodes ', num2str(mean(resultados_teste(:,2),'all'))])
plot(1:max_length-1,movmean(resultados_teste(:,2),[10,0],'omitnan'), 'linewidth',3)
grid on
xlabel('training epoch')
hold off


drawnow 

disp('____________________________')
final_mean_ftd(cont_vol) = mean(mean(resultados_teste(:,3),'all'))
final_mean_fp(cont_vol) = mean(mean(resultados_teste(:,2),'all'))

end
% save(['experiment_ws_snr_',num2str(snr),'_numep_',num2str(n_total_episodes),'_dt_','23-Apr-2024-18-42-09','.mat'])





