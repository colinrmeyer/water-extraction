% Set up SHAKTI simulation on Thwaites
% Beginning from Thwaites model shared by Mathieu, 28 Jan 2026
% Change to Budd frictio with uniform coefficient to test

load Thwaites_fric2_Transient_Calibration.mat % Mathieu's model
% load thwaites_inversion_N_shakti.mat % Inversion using preliminary SHAKTI N
md.miscellaneous.name = 'thwaites';

% Begin with this
% md.friction.C = md.results.TransientSolution(1).FrictionC;
md.inversion.iscontrol = 0;
md.autodiff.isautodiff = 0;
md.toolkits = toolkits();
md.toolkits.DefaultAnalysis=bcgslbjacobioptions(); % Iterative solver instead of MUMPS
md.verbose=verbose('solution',true);
md.settings.waitonlock=1;
md.timestepping=timestepping(); % Turn off adaptive timestepping
md.basalforcings.geothermalflux=0.05*ones(md.mesh.numberofvertices,1); % Need to add geothermal

% Turn on SHAKTI and turn off other transient processes for now
md.transient=deactivateall(md.transient);
md.transient.isstressbalance=1; % Solve for ice velocity 1, turn off 0
md.transient.ishydrology=1;



% HYDROLOGY SPECIFIC PARAMETERIZATION:
% Change hydrology class to SHAKTI model
md.hydrology=hydrologyshakti();

% Define distributed englacial input to the subglacial system (m/yr)
md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);

% Define initial water head such that water pressure is 50% of ice overburden pressure
md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

% Initial subglacial gap height of 0.001m everywhere
md.hydrology.gap_height = 0.01*ones(md.mesh.numberofelements,1);

% Typical bed bump bump spacing
md.hydrology.bump_spacing = 1.0*ones(md.mesh.numberofelements,1);

% Typical bed bump height
md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);

% Initial Reynolds number (start at Re=1000 everywhere)
md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);

% Deal with boundary conditions
md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
md.hydrology.spchead(md.mask.ocean_levelset<=0)=0; % Set as 5 threshold to include the small spots under floating ice
% md.hydrology.spchead(md.geometry.thickness<=100)=0;

md.hydrology.moulin_input = zeros(md.mesh.numberofvertices,1); % No moulin inputs
md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1); % No-flux b.c. on boundary except outflow

% Change to Budd friction law
%Calculate basal stress from Schoof
N    = md.materials.rho_ice*md.constants.g*md.geometry.thickness+md.materials.rho_water*md.constants.g.*md.geometry.base;
vb   = sqrt(md.results.StressbalanceSolution.Vx.^2+md.results.StressbalanceSolution.Vy.^2)/md.constants.yts;
% vb   = md.initialization.vel/md.constants.yts;
m    = averaging(md,md.friction.m,0);
C    = averaging(md,md.friction.C,0);
Cmax = averaging(md,md.friction.Cmax,0);
taub = (C.^2 .* vb.^(m-1))./(1 + (C.^2./(Cmax.*N)).^(1./m).*vb).^(m) .* vb;

%Change friction from Schoof to Budd
md.friction=friction();
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);
md.friction.coefficient=sqrt(taub./(N.*vb));
md.friction.coefficient(isnan(md.friction.coefficient))=0;

md.friction.effective_pressure=N;
md.friction.effective_pressure_limit(:)=0;
md.friction.coupling = 4;
% md.friction.coefficient(md.mask.ocean_levelset<=0)=0;


md.cluster = generic('np',8);

% Define the time stepping scheme
md.timestepping.time_step=1*3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=1/365; % Final time (in years)
md.settings.output_frequency=24; % To not save every time step in model output

md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation','HydrologyReynolds'};

md.verbose.solution=1;

md=solve(md,'Transient');