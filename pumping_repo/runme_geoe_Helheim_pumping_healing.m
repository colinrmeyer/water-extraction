clear all;close all
load Models/Helheim_refined_winter_1to2yr.mat

% Starting conditions
md.hydrology.head=md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end).HydrologyBasalFlux./1.787e-6;
md.hydrology.reynolds(md.hydrology.reynolds==0)=1;
md.friction.effective_pressure=md.results.TransientSolution(end).EffectivePressure;

md.initialization.vx=md.results.TransientSolution(end).Vx;
md.initialization.vy=md.results.TransientSolution(end).Vy;
md.initialization.vel=md.results.TransientSolution(end).Vel;

md.hydrology.englacial_input(:)=0.0;

% % TURN OFF PUMPING ***
% md.hydrology.moulin_input(:)=0.0;

% Time-stepping
md.timestepping.time_step=1/4*3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=270/365;
md.settings.output_frequency=4*24;
timevec=0:md.timestepping.time_step:md.timestepping.final_time;

% % *** Negative meltwater inputs for geoengineering ***********************
% % xp=308516; yp=-2577550; % Helheim near terminus *A*
% xp=299807; yp=-2576740; % Helheim confluence *B*
% % xp=297134; yp=-2577740; % Helheim lower branch *C* (old)
% % xp=293048; yp=-2564710; % Helheim peak vel top branch *D* (new C)
% % xp=280612; yp=-2578780; % Helheim top of south branch *E* (new D)
%xp=274213; yp=-2563810; % Interior point for pumping test
%xp=5000; yp=5000; % bigsyn thinner pumping point
%[a,pos] = min(sqrt((md.mesh.x-xp).^2+(md.mesh.y-yp).^2));
% md.hydrology.moulin_input(pos)=-0.3; % m3/s
% % ************************************************************************
% All 11 points
xp=[299807 308516 293048 280612 294682 296884 290196 282422 286457 285247 302502];
yp=[-2576740 -2577550 -2564710 -2578780 -2566250 -2572380 -2578010 -2579160 -2558610 -2579600 -2577520];
% *** Time-varying pumping and healing at 11 points***
md.hydrology.moulin_input=zeros(md.mesh.numberofvertices+1,length(timevec));
md.hydrology.moulin_input(end,:)=timevec;
for i=1:length(xp)
    [a,pos] = min(sqrt((md.mesh.x-xp(i)).^2+(md.mesh.y-yp(i)).^2));
for tt=1:length(timevec)
    if timevec(tt)>=90/365 && timevec(tt)<180/365
        md.hydrology.moulin_input(pos,tt)=-1.0;
    end
end
end

% *** Time-varying pumping and healing *** at one point
%md.hydrology.moulin_input=zeros(md.mesh.numberofvertices+1,length(timevec));
%md.hydrology.moulin_input(end,:)=timevec;
%for tt=1:length(timevec)
%    if timevec(tt)>=90/365 && timevec(tt)<180/365
%        md.hydrology.moulin_input(pos,tt)=-1.0;
%    end
%end
% *** 


md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation','HydrologyPmpHeat'};

md.cluster=generic('np',30);
%md.cluster.interactive=0;
%md.settings.waitonlock=0;
md.verbose.solution=1;
%md=solve(md,'Transient');

f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
f(f<0)=0;
Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

 md=loadresultsfromcluster(md,'runtimename','Helheim-08-07-2025-11-37-57-1519959');
 save('Models/Helheim_11points_pump_heal_dt15min','md')

