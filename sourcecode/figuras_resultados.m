%% SETUP
load('data_from_exp_results.mat')
lw = 2;
ms = 12;
resultados_teste = media_voluntarios; %squeeze(resultados_teste2(randi([1 11],1),:,:));
n_total_episodes = max_length-1;

%% PLOTS

figure(1)
subplot(2,2,[1 3])
% plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
% title(['snr = ',num2str(snr),'- average reward during test episodes'])

hold on
plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
'Markersize',ms,'Color',[0.5843    0.8157    0.9882])
plot(1:n_total_episodes,movmean(resultados_teste(:,1),[10 0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
grid on
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',18)
ylabel('Reward', 'Interpreter','latex', 'FontSize',18)
legend('Average results across 11 voluntaries at each 1s window','Moving average (last 10 values)', ...
    'Location', 'southeast', 'FontSize',12)
hold off

% figure(2)
% plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
subplot(2,2,2)
% title(['snr = ',num2str(snr),'- mean tfp during test episodes'])

hold on
plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
'Markersize',ms,'Color',[0.5843    0.8157    0.9882])

x = 1:n_total_episodes;
idx = find(movmean(resultados_teste(:,2),[10,0],'omitnan')<=5);
plot(x(idx),resultados_teste(idx,2),'.', ...
'Markersize',ms,'Color',[1    0.1    0.1])
plot(1:n_total_episodes,movmean(resultados_teste(:,2),[10,0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
grid on
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',18)
ylabel('False Positive Rate [\%]', 'Interpreter','latex', 'FontSize',18)
hold off

% figure(3)
% plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
subplot(2,2,4)
% title(['snr = ',num2str(snr),'- mean tdr during test episodes'])

hold on
plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
'Markersize',ms,'Color',[0.5843    0.8157    0.9882] )
plot(1:n_total_episodes,movmean(resultados_teste(:,3),[10 0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
grid on
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',18)
% ylabel('Detection Rate ', 'Interpreter','latex', 'FontSize',18)
ylabel('Detection Rate [\%]', 'Interpreter','latex', 'FontSize',18)
hold off



%%
% 
% time_below_5pct = nan(1,11);
% for ii= 1:size(resultados_teste2,1)
%     time_below_5pct(ii) = sum(resultados_teste2(ii,:,2)<=1)./numel(resultados_teste2(ii,:,2));
% end

% % Mean time below 5% false positive threshold in seconds:
mean(time_below_5pct*max_length)

% % of exam duration with valid responses:
100*mean(time_below_5pct*max_length)/max_length

% First valid detection:
idx_valid = find(media_voluntarios(:,2)<=1);
idx_valid(1)


%%
% 
% time_below_5pct = nan(1,11);
% for ii= 1:size(resultados_teste2,1)
%     time_below_5pct(ii) = sum(resultados_teste2(ii,:,2)<=1)./numel(resultados_teste2(ii,:,2));
% end

% % Mean time below 5% false positive threshold in seconds:
mean(time_below_5pct*max_length)

% % of exam duration with valid responses:
100*mean(time_below_5pct*max_length)/max_length

% First valid detection:
idx_valid = find(media_voluntarios(:,2)<=1);
idx_valid(1)

