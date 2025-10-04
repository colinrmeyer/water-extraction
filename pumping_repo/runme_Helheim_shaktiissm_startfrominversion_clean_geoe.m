% Coupled SHAKTI-ISSM simulation of Helheim, starting from inversion obtained using runme script (runme_Helheim_inversion)
clear all;close all

load Models/Helheim_refined_inversion.mat
% load Models/Helheim_refined_interior_inversion.mat
% load Models/Helheim_IS2_inversion.mat
% load Models/NW_IS2_inversion_17July.mat

% Turn off inversion
md.inversion.iscontrol=0;

% HYDROLOGY SPECIFIC PARAMETERIZATION:
% Change hydrology class to SHAKTI model
md.hydrology=hydrologyshakti();

% Define distributed englacial input to the subglacial system (m/yr)
md.hydrology.englacial_input = .0*ones(md.mesh.numberofvertices,1);

% Define initial water head such that water pressure is 80% of ice overburden pressure
md.hydrology.head = 0.8*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

% Initial subglacial gap height of 0.01m everywhere
md.hydrology.gap_height = 0.01*ones(md.mesh.numberofelements,1);

% Typical bed bump bump spacing
md.hydrology.bump_spacing = 1.0*ones(md.mesh.numberofelements,1);

% Typical bed bump height
md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);

% Initial Reynolds number (start at Re=1000 everywhere)
md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);

% Glacier front b.c.
md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);

% Set atmospheric pressure (head=0) everywhere <=10m ice thickness
pos=find(md.geometry.thickness<=10);
md.hydrology.spchead(pos)=0;

% Set pressure in fjord equal to hydrostatic pressure of fjord water
% NOTE: Can take this out in updated version of ISSM where SHAKTI ignores
% elements without ice
pos=find(md.mask.ice_levelset>0);
md.hydrology.spchead(pos)=0;

% % Deal with terminus b.c. in SHAKTI
% terminus=ContourToNodes(md.mesh.x,md.mesh.y,'Exp/terminus_big.exp',1);
% pos=find(terminus & md.mask.ice_levelset<=0);
% md.hydrology.spchead(pos)=0;

md.hydrology.moulin_input = zeros(md.mesh.numberofvertices,1); % No moulin inputs
md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1); % No-flux b.c. on boundary except outflow

% Coupling and friction
md.transient=deactivateall(md.transient);
md.transient.isstressbalance=1; % Solve for ice velocity
md.transient.ishydrology=1;

% Friction
Neff = md.materials.rho_ice*md.constants.g*md.geometry.thickness-md.materials.rho_water*md.constants.g*(md.hydrology.head - md.geometry.base);
md.friction.effective_pressure=Neff;
md.friction.coupling = 4; % Change this to couple (4) /uncouple with prescribed N (3) ***

% Specify that you want to run the model on your current computer
% Change the number of processors according to your machine
md.cluster=generic('np',8);

% Define the time stepping scheme
md.timestepping.time_step=3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=30/365; % Final time (in years)
md.settings.output_frequency=24; % To not save every time step in model output

% % Transient boundary condition at terminus (for termforce simulations)
% timevec=0:md.timestepping.time_step:md.timestepping.final_time;
% md.stressbalance.spcvx=NaN(md.mesh.numberofvertices+1,length(timevec));
% md.stressbalance.spcvx(end,:)=timevec;
% pos=find(terminus & md.mask.ice_levelset<=0);
% for i=1:length(timevec)
%     md.stressbalance.spcvx(pos,i)=-1000*sin(2*pi.*timevec(i))+md.initialization.vx(pos);
% end

md.verbose.solution=1;

md=solve(md,'Transient');
