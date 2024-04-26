%% CLEAN ENV
clearvars;  
% close all;
clc;

%% DEFINE PARAMS
%%% env
n_total_episodes = 2*1e3;
max_episode_length = 40;

%%% learning
alpha = 0.8;
gamma = 0.8;
epsilon = @(episode) 0.2*exp(-0.5*log10(episode));
% plot(0.2*exp(-0.5*log10(1:5e3)))

rewards_vec  = nan(1,n_total_episodes);
mean_tdr_vec = nan(1,n_total_episodes);
mean_tfp_vec = nan(1,n_total_episodes);
%% DEFINE AGENT

resolution = 10;
sssize = [resolution,resolution,resolution,resolution,2];
% sssize = [resolution,resolution,resolution,resolution,resolution,resolution,2];
qtable = rand(sssize);

%% TRAIN AGENT
signal_freq_bins =  [82    90    84    86    88    90    92    94    96];
noise_freq_bins = round(signal_freq_bins.*exp(1)/2)+5;
all_freq_bins = [signal_freq_bins,noise_freq_bins];

snr = @(episode) 5; %randi([-40,5]);

max_length = max_episode_length;
tdr = 0;
tfp = 0;

% size = 16 x 6 x max_length
snr_atual= snr(0);
all_states = rlord_gen_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length,...
                      resolution);
all_states(isnan(all_states))=1;

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
        % reward = +((5*tfp)^2)/(-100) ...
        %          +((tdr)^2)/(100)-el;   
        
        reward = +((tfp)^2)/(-100) ...
                 +((tdr)^2)/(100)-el;

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

    mean_tdr_vec(idx_window) = mean(tdr_hist);
    mean_tfp_vec(idx_window) = mean(tfp_hist);
    rewards_vec(idx_window) = mean(total_reward_for_this_episode);

    if mod(idx_window,100)==0
        figure(1)
        plot(1:n_total_episodes,rewards_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('average reward during episode')
        hold on
        plot(1:n_total_episodes,movmean(rewards_vec,[100 0]), 'linewidth',3)
        grid on
        hold off

        figure(2) % [0.2 0.5 0.9 0.2]
        plot(1:n_total_episodes,mean_tfp_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('mean tfp during episode')
        hold on
        plot(1:n_total_episodes,movmean(mean_tfp_vec,[100,0]), 'linewidth',3)
        grid on
        hold off

        figure(3)
        plot(1:n_total_episodes,mean_tdr_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('mean tdr during episode')
        hold on
        plot(1:n_total_episodes,movmean(mean_tdr_vec,[100 0]), 'linewidth',3)
        grid on
        hold off

        drawnow

        decisions

    end

    if mod(idx_window,200)==0
        disp(idx_window)
        snr_atual= snr(0);
    end
    if mod(idx_window,1)==0
        % t= tic()
        all_states = rlord_gen_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length,...
                      resolution);
        all_states(isnan(all_states))=1;
        % t= toc(t)
    end

end


%% TESTE
rewards_vec  = nan(1,n_total_episodes);
mean_tdr_vec = nan(1,n_total_episodes);
mean_tfp_vec = nan(1,n_total_episodes);

% size = 16 x 6 x max_length
snr_atual= snr(0);
all_states = rlord_gen_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length,...
                      resolution);
all_states(isnan(all_states))=1;


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
                [old_q, current_action] = max(current_q_sa);  
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
        el = 100*current_window/max_episode_length;
        reward = +((tfp)^2)/(-100) ...
                 +((tdr)^2)/(100)-el/1000;   


        % accumulate reward
        total_reward_for_this_episode(current_window) = reward;
        tdr_hist(current_window) = tdr;
        tfp_hist(current_window) = tfp;

    end

    mean_tdr_vec(idx_window) = mean(tdr_hist);
    mean_tfp_vec(idx_window) = mean(tfp_hist);
    rewards_vec(idx_window) = mean(total_reward_for_this_episode);

    if mod(idx_window,100)==0
        figure(1)
        plot(1:n_total_episodes,rewards_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('average reward during episode')
        hold on
        plot(1:n_total_episodes,movmean(rewards_vec,[100 0]), 'linewidth',3)
        grid on
        hold off

        figure(2) % [0.2 0.5 0.9 0.2]
        plot(1:n_total_episodes,mean_tfp_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('mean tfp during episode')
        hold on
        plot(1:n_total_episodes,movmean(mean_tfp_vec,[100,0]), 'linewidth',3)
        grid on
        hold off
        
        figure(3)
        plot(1:n_total_episodes,mean_tdr_vec, 'Color',[0.2 0.5 0.9 0.4])
        title('mean tdr during episode')
        hold on
        plot(1:n_total_episodes,movmean(mean_tdr_vec,[100 0]), 'linewidth',3)
        grid on
        hold off
 
        drawnow

        decisions

    end
    if mod(idx_window,5)==0
        snr_atual= snr(0);
    end
    if mod(idx_window,1)==0
        all_states = rlord_gen_states(signal_freq_bins, ...
                      noise_freq_bins, ...
                      snr_atual, ...
                      max_length,...
                      resolution);
        all_states(isnan(all_states))=1;
    end

end
