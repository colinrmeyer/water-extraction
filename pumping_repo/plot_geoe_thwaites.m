% Plot difference in velocity and N after pumping

% plotmodel(md,'data',md.results.TransientSolution(end).Vel-md0.results.TransientSolution(end).Vel,...
%     'title','Change in Vel (m/yr)',...
%     'data',md.results.TransientSolution(end).EffectivePressure-md0.results.TransientSolution(end).EffectivePressure,...
%     'title','Change in N (Pa)',...
%     'caxis#1',[-17 17],'caxis#2',[-100000 100000],'colormap#all',vik,...
%     'data',(md.results.TransientSolution(end).Vel-md0.results.TransientSolution(end).Vel)./md0.results.TransientSolution(end).Vel.*100,...
%     'title','% Change in Vel',...
%     'data',(md.results.TransientSolution(end).EffectivePressure-md0.results.TransientSolution(end).EffectivePressure)./md0.results.TransientSolution(end).EffectivePressure.*100,...
%     'title','% Change in N',...
%     'caxis#3',[-1 1],'caxis#4',[-5 5],...
%     'unit#all','km')

% % Plot velocity and pumping sites
% plotmodel(md,'data',md0.results.TransientSolution(end).Vel,'title','Ice Velocity (m/yr)','caxis',[0 3000],'unit','km','colormap',batlow)
% xp=[-1456120 -1491640 -1444290 -1407350 -1488260 -1514910 -1448710 -1410100 -1325170 -1348310 -1385190]; % 11 points
% yp=[-486265 -448776 -454887 -429720 -474660 -480518 -415668 -482792 -444682 -382003 -354339]; % 11 points
% hold on
% plot(xp./1000,yp./1000,'w.')

% % To plot just velocity (Fig. 1)
% % load nuuk.mat
% % plotmodel(md,'data',md.initialization.vel,'unit','km','colormap',batlow,'caxis',[0 3000],'fontsize',20) Thwaites
% plotmodel(md,'data',md.initialization.vel,'unit','km','colormap',batlow,'caxis',[0 7000],'fontsize',20); % Helheim

% % Plot effective presure
% plotmodel(md,'data',md.results.TransientSolution(end).EffectivePressure./1e6,'unit','km','colormap',batlow,'caxis',[0 3],'fontsize',20)

% % Plot ice thickness
% plotmodel(md,'data',md.geometry.thickness,'unit','km','colormap',batlow,'caxis',[0 3000],'fontsize',20) % Thwaites
% plotmodel(md,'data',md.geometry.thickness,'unit','km','colormap',batlow,'caxis',[0 2000],'fontsize',20) % Helheim

% Plot basal shear stress tau_b (Budd)
tau_b=md.friction.coefficient.^2.*md.results.TransientSolution(end).EffectivePressure.*md.results.TransientSolution(end).Vel./md.constants.yts;
plotmodel(md,'data',tau_b./1e3,'unit','km','colormap',batlow,'fontsize',20); % Plot tau_b in kPa


% Plot water flux
% plotmodel(md,'data',md.results.TransientSolution(1).HydrologyBasalFlux,'unit','km','colormap',batlow,'caxis',[1e-6 1],'fontsize',20,'log',10)
% plotmodel(md,'data',md.results.TransientSolution(1).HydrologyBasalFlux,'unit','km','colormap',batlow,'caxis',[1e-6 0.01],'fontsize',20)

% % Plot basal melt rate
% L=3.34e5; % Latent heat of fusion J/kg
% % Calculate geothermal contribution to melt (kg/m2/s)
% G=max(md.basalforcings.geothermalflux);
% Gheat=G/L;
% for i=1:length(md.results.TransientSolution)
%     
%     % Get melt rate in kg/m2/s
%     meltrate(:,i)=md.results.TransientSolution(i).HydrologyMeltRate;
%     
%     % Get frictional heat from sliding (kg/m2/s)
%     frictionheat(:,i)=md.results.TransientSolution(i).HydrologyFrictionHeat ./L;
%     
%     % Get dissipation heat (kg/m2/s)
%     dissipation(:,i)=md.results.TransientSolution(i).HydrologyDissipation ./L;
%     
%     % Get fractions of total melt
%     fracG(:,i)=Gheat./meltrate(:,i);
%     fracfric(:,i)=frictionheat(:,i)./meltrate(:,i);
%     fracdiss(:,i)=dissipation(:,i)./meltrate(:,i);
%     
%     % Convert meltrate to m/yr
%     meltrate_myr(:,i)=meltrate(:,i)./md.materials.rho_ice*3600*24*365;
%     
% end
% 
% plotmodel(md,'data',meltrate_myr(:,end),'unit','km','colormap',batlow,'caxis',[0 1],'fontsize',20)
% % plotmodel(md,'data',md.results.TransientSolution(1).HydrologyBasalFlux,'unit','km','colormap',batlow,'caxis',[1e-6 0.01],'fontsize',20)