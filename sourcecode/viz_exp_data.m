%% GET DATA
clearvars;
% close all;
clc;



% Vintensidade = {'70';'60';'50';'40';'30'};
cont_int = 3;
cont_vol = 3;

signal_freq_bins =  [82   82   84    86    88    90    92    94    96];
noise_freq_bins = round(signal_freq_bins.*exp(1)*sqrt(2))+12;
all_freq_bins = [signal_freq_bins,noise_freq_bins];
resolution = 10;
% size = 16 x 6 x max_length
all_states= db_to_dstates(signal_freq_bins, ...
                                noise_freq_bins, ...
                                resolution, ...
                                cont_vol, ...
                                cont_int);


%%

last_t = size(all_states,3);
fist_t = last_t-60;
% fist_t = round(size(all_states,3)/3);

fontsize = 15;
figure(2)
subplot(131)
s = stem3(all_states(:,1,fist_t:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1);
hold on
stem3(all_states(1:9,1,fist_t:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5);


view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [CSM]', 'FontSize',fontsize, 'Interpreter','latex')
hold off



subplot(132)
stem3(all_states(:,2,fist_t:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1)
% [lines(1) .5]
hold on
stem3(all_states(1:9,2,fist_t:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5)

view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [GFT]', 'FontSize',fontsize, 'Interpreter','latex')
hold off


% figure(3)
% close all
subplot(133)
stem3(all_states(:,3,fist_t:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1)

hold on
stem3(all_states(1:9,3,fist_t:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5)

view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [MSC]', 'FontSize',fontsize, 'Interpreter','latex')
hold off


