clear all;close all

steps=[1:3];

if any(steps==1) 
	disp('	Step 1: Mesh');

% 	%Generate unstructured mesh on 2x10km box with typical element edge length of 100 m
% 	md=triangle(model,'Exp/geoe_box.exp',200);
% 
%     % Refine mesh according to radius around pumping site at x=2000, y=1000
%     r=sqrt((md.mesh.x-2000).^2+(md.mesh.y-1000).^2);
% %     h=NaN*ones(md.mesh.numberofvertices,1);
% %     h(r<=500)=10;
%     h=1/40.*r+20;
%     md=bamg(md,'err',50,'hmin',20,'hmax',200,'hVertices',h);

        % For round domain
    md=roundmesh(model,10000,100);
%md=roundmesh(model,1000,10);

        % Refine mesh according to radius around pumping site at x=0, y=0
    r=sqrt((md.mesh.x-0).^2+(md.mesh.y-0).^2);
    h=1/25.*r+10;
%    h=1/25.*r+1;
    md=bamg(md,'err',50,'hmin',10,'hmax',400,'hVertices',h);

	save MoulinMesh md
end 

if any(steps==2) 
	disp('	Step 2: Parameterization');
	md=loadmodel('MoulinMesh');

	md=setmask(md,'','');

	% Run parameterization script to set up geometry, velocity, material properties, etc.
	md=parameterize(md,'Par/geoe_box.par');

	% HYDROLOGY SPECIFIC PARAMETERIZATION:
	% Change hydrology class to Sommers' SHaKTI model
	md.hydrology=hydrologyshakti();

	% Define initial water head such that water pressure is 50% of ice overburden pressure
	md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

	% Initial subglacial gap height of 0.01m everywhere
	md.hydrology.gap_height = 0.1*ones(md.mesh.numberofelements,1);

	% Typical bed bump bump spacing (2m)
	md.hydrology.bump_spacing = 1*ones(md.mesh.numberofelements,1);

	% Typical bed bump height (0.1m)
	md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);

	% Define distributed englacial input to the subglacial system (m/yr)
	% Change the value 0.0 to add distributed input
	md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);

	% Initial Reynolds number (start at Re=1000 everywhere)
	md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);

	% Set Dirichlet boundary condition at outer edge of
	% domain (outflow, i.e. N=100kPa --> h=(1e5-rhoi*g*H)/(-rhow*g)+zb)
	md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
	pos=find(md.mesh.vertexonboundary);
%     md.hydrology.spchead(pos)=md.geometry.base(pos); % Atomspheric pressure
	md.hydrology.spchead(pos)=(2e5-md.materials.rho_ice*md.constants.g*md.geometry.thickness(pos))./(-md.materials.rho_freshwater*md.constants.g)+md.geometry.base(pos);

	save MoulinParam md;
end 

if any(steps==3) 
	disp('	Step 3: Solve!');
	md=loadmodel('MoulinParam');

	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;

	% Specify that you want to run the model on your current computer
	% Change the number of processors according to your machine (here np=4)
%	md.cluster=generic('np',8);

	% Define the time stepping scheme: run for 90 days with a time step of 1 hr
	md.timestepping.time_step=4*3600/md.constants.yts; % Time step (in years)
	md.timestepping.final_time=5*365/365;
    md.settings.output_frequency=6;

% 	%Add one pumping site with steady extraction at x=2000, y=1000
% 	[a,pos] = min(sqrt((md.mesh.x-2000).^2+(md.mesh.y-1000).^2));
%	time=0:md.timestepping.time_step:md.timestepping.final_time;
% 	md.hydrology.moulin_input = zeros(md.mesh.numberofvertices+1,numel(time));
% 	md.hydrology.moulin_input(end,:)=time;
% 	md.hydrology.moulin_input(pos,:)=-2;
md.hydrology.moulin_input=zeros(md.mesh.numberofvertices,1);
[a,pos] = min(sqrt((md.mesh.x-0).^2+(md.mesh.y-0).^2));
md.hydrology.moulin_input(pos)=-1;

	% Specify no-flux Type 2 boundary conditions on all edges (except
	% the Type 1 condition set at the outflow above)
%	md.hydrology.neumannflux=zeros(md.mesh.numberofelements+1,numel(time));
%	md.hydrology.neumannflux(end,:)=time;
md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);

    % Save melt rate and components
    md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation','HydrologyPmpHeat'};

	md.cluster=generic('np',20);
%md.cluster.interactive=0;
%md.settings.waitonlock=0;
md.verbose.solution=1;
verbose('all');
%md=solve(md,'Transient');

%f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
%f(f<0)=0;
%Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

md=loadresultsfromcluster(md,'runtimename','geoe_box-11-20-2024-12-11-43-78348');
%save('Models/Helheim_SI_brlr_1yr','md','f','Re')

	save MoulinTransient md
end 
