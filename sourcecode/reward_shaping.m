clearvars;
close all;
clc;

n_tiles = 25;
fp = linspace(1e-4,100,n_tiles);
dr = linspace(1e-4,100,n_tiles);
el = [100];
fpdes  = 5;
% el = el(1);
c = lines(numel(el));
s = [];
r = [];

% quero r<0 se fp>fpdes , para todo dr
% -(abs(dr(ifp))/abs(fp(ifp)-fpdes ))
% -abs(fp(ifp)-fpdes)
% -(dr(idr)*fp(ifp)/(100*fpdes))

for iel = 1:numel(el)
    for idr = 1:numel(dr)
        for ifp = 1:numel(fp)
            % paraboloide
            r(idr,ifp) = +(  ((fp(ifp))^2)/(-100) ...
                 +((dr(idr))^2)/(100) -el(iel));

            % outras propostas
            r(idr,ifp) = +( -(((100-fpdes)^2)/(100)+(fpdes^2)/-100)...
                +((dr(idr))^2)/(-100) ...
                 +((fp(ifp))^2)/(100) ...
                 )*(el(iel)/100);
            % r(idr,ifp) =(+((dr(idr))^2)/(-100) ...
            %      +((fp(ifp)-fpdes)^2)/(100) ...
            %      )*(el(iel)/100);

            % r(idr,ifp) = +(  ((fp(ifp))^2)/(-100) ...
            %      +((dr(idr))^2)/(100) )*el(iel);
            
        end
    end
    s(iel) = surf(dr,fp,r,'FaceColor',c(iel,:), ...
        'EdgeColor','none',...
        'FaceAlpha',.7);
    hold on
    
end

xlabel('True Positive Rate [TPR \%]','Interpreter', 'Latex','Fontsize',14,'Rotation',25)
ylabel('False Positive Rate [FPR \%]','Interpreter', 'Latex','Fontsize',14,'Rotation',-18)
zlabel('Reward','Interpreter', 'Latex','Fontsize',14)

% legend([s(1),s(2),s(3),s(4),s(5)], ...
%     'Episode Progress = 0\%', ...
%     'Episode Progress = 20\%', ...
%     'Episode Progress = 50\%', ...
%     'Episode Progress = 80\%',...
%     'Episode Progress = 100\%', ...
%     'Interpreter', 'Latex',...
%     'Fontsize',10,...
%     'Location', 'northwest')

% ylim([0,20])


