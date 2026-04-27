clear all;close all
%load thwaites_shakti_1yr.mat 
load Models/thwaites_budd_1d.mat
md.miscellaneous.name='thwaites';

md.inversion.iscontrol=0;

% Starting conditions
md.hydrology.head=md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end).HydrologyBasalFlux./1.787e-6;
md.hydrology.reynolds(md.hydrology.reynolds==0)=1;
md.friction.effective_pressure=md.results.TransientSolution(end).EffectivePressure;

 md.transient.isstressbalance=1; % Turn on velocity coupling
 md.initialization.vx=md.results.TransientSolution(end).Vx;
 md.initialization.vy=md.results.TransientSolution(end).Vy;
 md.initialization.vel=md.results.TransientSolution(end).Vel;

md.hydrology.englacial_input(:)=0.0;

% Time-stepping
md.timestepping.time_step=1*3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=365/365;
md.settings.output_frequency=24;

md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation'};

md.cluster=generic('np',30);
%md.cluster.interactive=0;
%md.settings.waitonlock=0;
md.verbose.solution=1;
md=solve(md,'Transient');

f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
f(f<0)=0;
Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

% md=loadresultsfromcluster(md,'runtimename','IS-07-23-2025-14-52-07-957867');
% save('Models/NW_IS_23July_1to2yr','md','-v7.3')

