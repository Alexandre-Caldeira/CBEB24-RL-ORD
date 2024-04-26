load('worskpace_res_pausado.mat')
figure(1)
subplot(2,2,[1 3])
plot(1:n_total_episodes,resultados_teste(:,1),'.', ...
'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
% title(['snr = ',num2str(snr),'- average reward during test episodes'])
hold on
plot(1:n_total_episodes,movmean(resultados_teste(:,1),[10 0],'omitnan'), 'linewidth',3,'Color',[0.8500    0.3250    0.0980 0.8])
grid on
xlabel('Training Episode', 'FontSize',14, 'Interpreter','latex')
ylabel('Average Reward', 'FontSize',14, 'Interpreter','latex')
legend('Average during 100 test episodes','Moving mean (last 10 values)', ...
    'Location', 'Southeast')
hold off
% figure(2)
subplot(2,2,2)
plot(1:n_total_episodes,resultados_teste(:,2),'.', ...
'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
hold on
% title(['snr = ',num2str(snr),'- mean tfp during test episodes'])
plot(1:n_total_episodes,movmean(resultados_teste(:,2),[10,0],'omitnan'), 'linewidth',3,'Color',[0.8500    0.3250    0.0980 0.8])
grid on
xlabel('Training Episode', 'FontSize',14, 'Interpreter','latex')
ylabel('Average  FPR [\%]', 'FontSize',14, 'Interpreter','latex')
hold off
% figure(3)
subplot(2,2,4)
plot(1:n_total_episodes,resultados_teste(:,3),'.', ...
'Markersize',10,'Color',[0.2 0.5 0.9 0.4])
hold on
% title(['snr = ',num2str(snr),'- mean tdr during test episodes'])
plot(1:n_total_episodes,movmean(resultados_teste(:,3),[10 0],'omitnan'), 'linewidth',3,'Color',[0.8500    0.3250    0.0980 0.8])
grid on
xlabel('Training Episode', 'FontSize',14, 'Interpreter','latex')
ylabel('Average DR [\%]', 'FontSize',14, 'Interpreter','latex')
hold off
drawnow