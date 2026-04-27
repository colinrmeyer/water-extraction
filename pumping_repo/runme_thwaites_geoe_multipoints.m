clear all;close all

load Models/thwaites_budd_1yr.mat

% Turn off inversion
md.inversion.iscontrol=0;

% Starting conditions
md.hydrology.head=md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end).HydrologyBasalFlux./1.787e-6;
md.hydrology.reynolds(md.hydrology.reynolds==0)=1;

md.initialization.vx=md.results.TransientSolution(end).Vx;
md.initialization.vy=md.results.TransientSolution(end).Vy;
md.initialization.vel=md.results.TransientSolution(end).Vel;

md.hydrology.englacial_input(:)=0.0;

% *** Negative meltwater inputs for pumping ***********************
%xp=-1456120; % Thwaites 1
%yp=-486265; % Thwaites 1
%[a,pos] = min(sqrt((md.mesh.x-xp).^2+(md.mesh.y-yp).^2));
% md.hydrology.moulin_input(pos)=-1.0; % m3/s
 

% Multiple points (Thwaites)
% xp=[-1456120 -1491640 -1444290 -1407350 -1488260 -1514910]; % 6 points
% yp=[-486265 -448776 -454887 -429720 -474660 -480518]; % 6 points
xp=[-1456120 -1491640 -1444290 -1407350 -1488260 -1514910 -1448710 -1410100 -1325170 -1348310 -1385190]; % 11 points
yp=[-486265 -448776 -454887 -429720 -474660 -480518 -415668 -482792 -444682 -382003 -354339]; % 11 points
for i=1:length(xp)
   [a,pos] = min(sqrt((md.mesh.x-xp(i)).^2+(md.mesh.y-yp(i)).^2));
   md.hydrology.moulin_input(pos)=-1.0; % m3/s
end
% ************************************************************************

% Time-stepping
md.timestepping.time_step=1*3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=365/365;
md.settings.output_frequency=24;

%md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation','HydrologyReynolds'};,

md.cluster=generic('np',30);
%md.cluster.interactive=0;
%md.settings.waitonlock=0;
md.verbose.solution=1;
md=solve(md,'Transient');

%f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
%f(f<0)=0;
%Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

% md=loadresultsfromcluster(md,'runtimename','thwaites-02-26-2026-11-51-26-4088267');
% save('Models/thwaites_C2000_multipoints','md')

