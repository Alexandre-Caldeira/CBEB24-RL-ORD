% clearvars;  
% close all;
% clc;

%% CLEAN ENV
clearvars -except snr_vec idx_snr  

snr = -8;

%% DEFINE PARAMS
%%% env
% n_total_episodes = 1e3;
n_total_episodes = 500;
max_episode_length = 40;

%%% learning
alpha = 5e-3;
gamma = 0.5;
epsilon = @(episode) 0.2*exp(-0.5*log10(episode));
% plot(0.2*exp(-0.5*log10(1:5*1e3)))

rewards_vec  = nan(1,n_total_episodes);
mean_tdr_vec = nan(1,n_total_episodes);
mean_tfp_vec = nan(1,n_total_episodes);
%% DEFINE AGENT

resolution = 10;
sssize = [resolution,resolution,resolution,resolution,2];
qtable = zeros(sssize);

%% TRAIN AGENT
signal_freq_bins =  82; %[82   90    84    86    88    90    92    94    96];
noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5;
all_freq_bins = [signal_freq_bins,noise_freq_bins];

max_length = max_episode_length;
tdr = 0;
tfp = 0;

% size = 16 x 6 x max_length
snr_atual= snr;
all_states = rlord_gen_log_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length,...
                      resolution);
all_states(isnan(all_states))=1;

resultados_teste =  nan(n_total_episodes,3);

% for each episode
for idx_window = 1:n_total_episodes

 
    % initialize episode parameters
    terminal_state = false;
    total_reward_for_this_episode = zeros(1,max_episode_length-1);
    tdr_hist= zeros(1,max_episode_length-1);
    tfp_hist= zeros(1,max_episode_length-1);
    is_freq_undecided = ones(1,numel(all_freq_bins));

    for current_window = 1:max_episode_length-1  

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
    
                current_q_sa = [qtable(sa1),qtable(sa2)];
                 
                nexts_a1 = sub2ind(sssize, ...
                    next_states(1),next_states(2),...
                    next_states(3),next_states(4),...
                    1);
    
                nexts_a2 = sub2ind(sssize, ...
                    next_states(1),next_states(2),...
                    next_states(3),next_states(4),...
                    2);
    
                next_q_sa = reshape(qtable([nexts_a1,nexts_a2]),1,2);
    
                % get eps greedy action
                if epsilon(idx_window)<=rand()
                    current_action = randi([1,2],1);
                    old_q = current_q_sa(current_action);
                else
                    [old_q, current_action] = max(current_q_sa);
                end
    
    
                
                max_q = max(next_q_sa(current_action));
    
                decisions(idx_freq,1:3) = [old_q,max_q, current_action];
                decisions(idx_freq,4:5) = current_q_sa;
    
    
                % if current_freq not in noise_freq_bins, should detect
                should_detect = isempty(find(noise_freq_bins==current_freq));
                if current_action==2 || current_window==max_episode_length
    
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
        el = current_window/max_episode_length;

        reward = +(((tfp)^2)/(-100) ...
                 +((tdr)^2)/(100))*el;

        for idx_freq = 1:numel(all_freq_bins)
            if is_freq_undecided(idx_freq)

            % decisions(idx_freq,1:3) = [old_q,max_q, current_action];
            % decisions(idx_freq,4:5) = current_q_sa;
            
            old_q = decisions(idx_freq,1);
            max_q = decisions(idx_freq,2);
            current_action = decisions(idx_freq,3);
            current_states = all_states(idx_freq,:,current_window);

            new_q =  old_q + alpha * (reward + gamma * max_q - old_q);

            selected_sa = sub2ind(sssize, ...
                current_states(1),current_states(2),...
                current_states(3),current_states(4),...
                current_action);
  
            qtable(selected_sa) = new_q;
            end

        end
        

        % accumulate reward
        total_reward_for_this_episode(current_window) = reward;
        tdr_hist(current_window) = tdr;
        tfp_hist(current_window) = tfp;

    end

    % salva resultados 
    mean_tdr_vec(idx_window) = mean(tdr_hist);
    mean_tfp_vec(idx_window) = mean(tfp_hist);
    rewards_vec(idx_window) = mean(total_reward_for_this_episode);

    % testa por 100 episodios
    if mod(idx_window,5)==0
        %% TESTE
        n_total_episodes_teste =100;
        rewards_vec_teste  = nan(1,n_total_episodes_teste);
        mean_tdr_vec_teste = nan(1,n_total_episodes_teste);
        mean_tfp_vec_teste = nan(1,n_total_episodes_teste);
        
        % size = 16 x 6 x max_length
        snr_atual= snr;
        signal_freq_bins =  [82 82   90    84    86    88    90    92    94    96];
        noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5;
        all_freq_bins = [signal_freq_bins,noise_freq_bins,];
        
        max_length_teste = 40;
        all_states = rlord_gen_log_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length_teste,...
                      resolution);
        all_states(isnan(all_states))=1;
       
        
        % for each episode
        for idx_window_teste = 1:n_total_episodes_teste
        
         
            % initialize episode parameters
            terminal_state = false;
            total_reward_for_this_episode = zeros(1,max_length_teste-1);
            tdr_hist= zeros(1,max_length_teste-1);
            tfp_hist= zeros(1,max_length_teste-1);
            is_freq_undecided = ones(1,numel(all_freq_bins));
    
            test_q = qtable;
            test_alpha = 1e-2;
            test_gamma = 0.5;
            
        
            for current_window = 1:max_length_teste-1  
        
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
            
                        current_q_sa = [qtable(sa1),qtable(sa2)];
                         
                        nexts_a1 = sub2ind(sssize, ...
                            next_states(1),next_states(2),...
                            next_states(3),next_states(4),...
                            1);
            
                        nexts_a2 = sub2ind(sssize, ...
                            next_states(1),next_states(2),...
                            next_states(3),next_states(4),...
                            2);
            
                        next_q_sa = reshape(qtable([nexts_a1,nexts_a2]),1,2);
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
                reward = +(((tfp)^2)/(-100) ...
                         +((tdr)^2)/(100))*el;   

    
            test_reward = +(((tfp)^2)/(-100) ...
                     +((tdr)^2)/(100))*el;

            
            for idx_freq = 1:numel(all_freq_bins)

                % if current_freq not in noise_freq_bins, should detect
                should_detect = isempty(find(noise_freq_bins==current_freq));

                if ~should_detect %is_freq_undecided(idx_freq)
    
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
                end
    
            end
        

        % accumulate reward
        total_reward_for_this_episode(current_window) = reward;
        tdr_hist(current_window) = tdr;
        tfp_hist(current_window) = tfp;
        
        
                % accumulate reward
                total_reward_for_this_episode(current_window) = reward;
                tdr_hist(current_window) = tdr;
                tfp_hist(current_window) = tfp;
        
            end
        
            mean_tdr_vec_teste(idx_window_teste) = mean(tdr_hist);
            mean_tfp_vec_teste(idx_window_teste) = mean(tfp_hist);
            rewards_vec_teste(idx_window_teste) = mean(total_reward_for_this_episode);
        
            
            signal_freq_bins =  [82 82   90    84    86    88    90    92    94    96];
            noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5;
            all_freq_bins = [signal_freq_bins,noise_freq_bins];

            all_states = rlord_gen_log_states(signal_freq_bins, ...
                          noise_freq_bins, ...
                          snr_atual, ...
                          max_length_teste,...
                          resolution);
            all_states(isnan(all_states))=1;
       
        
        end
    
        resultados_teste(idx_window,3)=mean(mean_tdr_vec_teste);
        resultados_teste(idx_window,2)=mean(mean_tfp_vec_teste);
        resultados_teste(idx_window,1)=mean(rewards_vec_teste);

        if mod(idx_window,25)==0
        figure(1)
        subplot(2,2,[1 3])
        plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
            'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
        title(['snr = ',num2str(snr),'- average reward during test episodes'])
        hold on
        plot(1:n_total_episodes,movmean(resultados_teste(:,1),[10 0],'omitnan'), 'linewidth',3)
        grid on
        xlabel('training epoch')
        hold off
    
        % figure(2)
        subplot(2,2,2)
        plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
            'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
        hold on
        title(['snr = ',num2str(snr),'- mean tfp during test episodes'])
        plot(1:n_total_episodes,movmean(resultados_teste(:,2),[10,0],'omitnan'), 'linewidth',3)
        grid on
        xlabel('training epoch')
        hold off
    
        % figure(3)
        subplot(2,2,4)
        plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
            'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
        hold on
        title(['snr = ',num2str(snr),'- mean tdr during test episodes'])
        plot(1:n_total_episodes,movmean(resultados_teste(:,3),[10 0],'omitnan'), 'linewidth',3)
        grid on
        xlabel('training epoch')
        hold off
        
        drawnow 
        end
    end


    % cria nova serie temporal
    signal_freq_bins =  82; %[82   90    84    86    88    90    92    94    96];
    noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5;
    all_freq_bins = [signal_freq_bins,noise_freq_bins];

    all_states = rlord_gen_log_states(signal_freq_bins, ...
                  noise_freq_bins, ...
                  snr, ...
                  max_length,...
                  resolution);
    all_states(isnan(all_states))=1;

end

% save(['log10states_ws_snr_',num2str(snr),'_numep_',num2str(n_total_episodes),'_dt_','23-Apr-2024-14-26-09','.mat'])




