% Spatial plots for pumping paper
% Helheim 11 sites pump/heal

% load Models/Helheim_11points_pump_heal.mat
% load Models/Helheim_refined_pump_heal.mat

load colormaps/batlow.mat
load colormaps/vik.mat

% % Pumping sites (11 points)
% xp=[308516 299807 293048 280612 294682 296884 290196 282422 286457 285247 302502]; 
% yp=[-2577550 -2576740 -2564710 -2578780 -2566250 -2572380 -2578010 -2579160 -2558610 -2579600 -2577520];

% Pumping site (confluence)
xp=299807;
yp=-2576740;

% % Plot initial velocity with points overlaid (rainbow sprinkles)
% % plotmodel(md,'data',md.results.TransientSolution(end).HydrologyBasalFlux)
% plotmodel(md,'data',md.results.TransientSolution(1).Vel,'unit','km','fontsize',14,...
%     'xlim',[275 311],'ylim',[-2585 -2555]);colormap(batlow)
% % hold on;plot(xp./1000,yp./1000,'w.',xp./1000,yp./1000,'ko')
% 
% for i=1:length(xp)
%     hold on;plot(xp(i)./1000,yp(i)./1000,'.',xp(i)./1000,yp(i)./1000,'wo','MarkerSize',20)
% end


% % Plot change in velocity
% plotmodel(md,'data',md.results.TransientSolution(180).Vel-md.results.TransientSolution(1).Vel,...
%     'unit','km','fontsize',14,'colormap',vik,'caxis',[-100 100],...
%     'xlim',[275 311],'ylim',[-2585 -2555]);
% for i=1:length(xp)
%     hold on;plot(xp(i)./1000,yp(i)./1000,'.',xp(i)./1000,yp(i)./1000,'wo','MarkerSize',20)
% end

% Plot change in velocity -- little caxis limits for confluence
plotmodel(md,'data',md.results.TransientSolution(180).Vel-md.results.TransientSolution(1).Vel,...
    'unit','km','fontsize',14,'colormap',vik,'caxis',[-1 1],...
    'xlim',[275 311],'ylim',[-2585 -2555]);
for i=1:length(xp)
%     hold on;plot(xp(i)./1000,yp(i)./1000,'.',xp(i)./1000,yp(i)./1000,'wo','MarkerSize',20)
hold on;plot(xp(i)./1000,yp(i)./1000,'w*','MarkerSize',10)
end

% % Plot change in effective pressure
% plotmodel(md,'data',(md.results.TransientSolution(180).EffectivePressure-md.results.TransientSolution(1).EffectivePressure)./1e6,...
%     'unit','km','fontsize',14,'colormap',vik,'caxis',[-2 2],...
%     'xlim',[275 311],'ylim',[-2585 -2555]);
% for i=1:length(xp)
% %     hold on;plot(xp(i)./1000,yp(i)./1000,'.',xp(i)./1000,yp(i)./1000,'wo','MarkerSize',20)
% hold on;plot(xp(i)./1000,yp(i)./1000,'w*','MarkerSize',10)
% end


% % Plot change in gap height
% plotmodel(md,'data',(md.results.TransientSolution(180).HydrologyGapHeight-md.results.TransientSolution(1).HydrologyGapHeight),...
%     'unit','km','fontsize',14,'colormap',vik,'caxis',[-.005 .005],...
%     'xlim',[275 311],'ylim',[-2585 -2555]);
% for i=1:length(xp)
%     hold on;plot(xp(i)./1000,yp(i)./1000,'.',xp(i)./1000,yp(i)./1000,'wo','MarkerSize',20)
% end
