% runme script to perform a stress balance inversion

steps = [1:3];
clustername = oshostname();
loadonly=1;

%Cluster parameters{{{
if strcmpi(clustername,'pfe'),
	cluster=pfe('numnodes',1,'time',60,'processor','bro','cpuspernode',28,'queue','devel'); %max time is 120 (2hr) and max cpuspernode is 28 for 'bro'
	cluster=pfe('numnodes',1,'time',60,'processor','bro','cpuspernode',28,'queue','normal');
elseif strcmpi(clustername,'greenplanet'),
	cluster=greenplanet('numnodes',3,'cpuspernode',22,'port',0);
else
	cluster=generic('name',oshostname(),'np',15);
end
clear clustername 
%}}}
org=organizer('repository',['./Models'],'prefix',['Model_Helheim_'],'steps',steps); clear steps;

if perform(org,'Mesh'),% {{{

	md=triangle(model,['Exp/Helheim_CG.exp'],500);

	[velx vely]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
	vel  = sqrt(velx.^2+vely.^2);

% 	%refine mesh using surface velocities as metric
% 	md=bamg(md,'hmin',100,'hmax',1500,'field',vel,'err',5);
% 	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
% 	md.mesh.epsg=3413;

    % Refine mesh according to radius around pumping site at xp,yp
%     xp=299807; yp=-2576740; % Helheim confluence *B*
    xp=274213; yp=-2563810; % Interior point
    [a,pos] = min(sqrt((md.mesh.x-xp).^2+(md.mesh.y-yp).^2));
    r=sqrt((md.mesh.x-xp).^2+(md.mesh.y-yp).^2);
    h=1/25.*r+10;
    md=bamg(md,'err',50,'hmin',10,'hmax',400,'hVertices',h);


	savemodel(org,md);
end %}}}
if perform(org,'Param'),% {{{

	md=loadmodel(org,'Mesh');
	md=setflowequation(md,'SSA','all');
	
	md=setmask(md,'','');
	md=parameterize(md,'../shaktiissm/Par/Greenland.par');
	md.miscellaneous.name = 'Helheim';

	savemodel(org,md);
end%}}}
if perform(org,'Inversion_drag'),% {{{

	md=loadmodel(org,'Param');
	
	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);
	md.transient.amr_frequency = 0;
	
	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=5000;
	md.inversion.cost_functions_coefficients(:,2)=10;
	md.inversion.cost_functions_coefficients(:,3)=2*50^-3;
	pos=find(md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=100;
	md.inversion.maxiter =100;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=4000*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=cluster;

	md=solve(md,'sb');

	%Put results back into the model
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;

	savemodel(org,md);
end%}}}