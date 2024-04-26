lw = 2;
ms = 15;

figure(1)
% plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
title(['snr = ',num2str(snr),'- average reward during test episodes'])
plot(1:n_total_episodes,movmean(resultados_teste(:,1),[10 0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
hold on
plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
'Markersize',ms,'Color',[0.2 0.5 0.9 0.4])
grid on
xlabel('Training Episode', 'Interpreter','latex', 'FontSize',18)
ylabel('Reward', 'Interpreter','latex', 'FontSize',18)
legend('Moving average (last 10 values)','Average during 100 test episodes', ...
    'Location', 'southeast')
hold off
figure(2)
% plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])

title(['snr = ',num2str(snr),'- mean tfp during test episodes'])
plot(1:n_total_episodes,movmean(resultados_teste(:,2),[10,0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
hold on
plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
'Markersize',ms,'Color',[0.2 0.5 0.9 0.4])
grid on
xlabel('Training Episode', 'Interpreter','latex', 'FontSize',18)
ylabel('False Positive Rate [\%]', 'Interpreter','latex', 'FontSize',18)
hold off
figure(3)
% plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
% 'Markersize',10,'Color',[0.2 0.5 0.9 0.4])

title(['snr = ',num2str(snr),'- mean tdr during test episodes'])
plot(1:n_total_episodes,movmean(resultados_teste(:,3),[10 0],'omitnan'), 'linewidth',lw,...
    'Color',[0.8500    0.3250    0.0980 1])
hold on
plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
'Markersize',ms,'Color',[0.2 0.5 0.9 0.4])
grid on
xlabel('Training Episode', 'Interpreter','latex', 'FontSize',18)
% ylabel('Detection Rate ', 'Interpreter','latex', 'FontSize',18)
ylabel('Detection Rate [\%]', 'Interpreter','latex', 'FontSize',18)
hold off