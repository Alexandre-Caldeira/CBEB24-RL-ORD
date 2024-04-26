%% GET DATA
% Vintensidade = {'70';'60';'50';'40';'30'};


clearvars;
% close all;
clc;
% wtf happens at
% cont_vol = 9;
% cont_int = 3;

resolution = 10;
% last_t = 55;
cont_vol = 3;
cont_int = 1;

signal_freq_bins =  [82   82   84    86    88    90    92    94    96];
noise_freq_bins = round(signal_freq_bins.*exp(1)*sqrt(2))+12;
all_freq_bins = [signal_freq_bins,noise_freq_bins];
% 73 14
% 6818
% size = 16 x 6 x max_length
all_states= db_to_dstates(signal_freq_bins, ...
                                noise_freq_bins, ...
                                resolution, ...
                                cont_vol, ...
                                cont_int);

last_t = size(all_states,3);
% all_states(isnan(all_states))=1;

%%

fontsize = 15;
figure(2)
subplot(131)
s = stem3(all_states(:,1,1:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1);
% [lines(1) .5]
hold on
stem3(all_states(1:8,1,1:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5);

% close all
% stem3(all_states(:,1,1:last_t))
% hold on
% stem3(all_states(1:8,1,1:last_t),'r')
view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [CSM]', 'FontSize',fontsize, 'Interpreter','latex')



% figure(2)
% close all
subplot(132)
stem3(all_states(:,2,1:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1)
% [lines(1) .5]
hold on
stem3(all_states(1:8,2,1:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5)

view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [GFT]', 'FontSize',fontsize, 'Interpreter','latex')



% figure(3)
% close all
subplot(133)
stem3(all_states(:,3,1:last_t),':', ...
    'MarkerFaceColor',lines(1), ...
    'Color',[.7 .7 .7],'LineWidth',.1)
% [lines(1) .5]
hold on
stem3(all_states(1:8,3,1:last_t),':' ...
    ,'MarkerFaceColor','r', ...
    'Color',[.7 .7 .7],'LineWidth',.5)

% stem3(all_states(:,2,1:last_t),'MarkerFaceColor',lines(1))
% hold on
% stem3(all_states(1:8,2,1:last_t),'r')
% ylabel('Frequency bins [1-8 = stim (red), 9-16 = noise (blue)]', 'FontSize',fontsize)
view(73,14)
ylabel('Frequency bins', 'FontSize',fontsize, 'Interpreter','latex')
xlabel('Time [s]', 'FontSize',fontsize, 'Interpreter','latex')
zlabel('Discretized log10(Energy) [MSC]', 'FontSize',fontsize, 'Interpreter','latex')