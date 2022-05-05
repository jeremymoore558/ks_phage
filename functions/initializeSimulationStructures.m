function [GlobalParams, SimParams] = initializeSimulationStructures(GlobalParams, SimParams)
% a bit awkward, but must be called from the ks folder because of the way
% the paths are define using Matlab's pwd function

%% default parameter values
% if reSimulate = 1, the simulation will be run, and any previous
% results.mat file will be overwritten
if ~isfield(GlobalParams,'reSimulate')
    GlobalParams.reSimulate = 0;
end

% initial cell densituy (in units of optical density) if the default 
% initial condition is used
if ~isfield(GlobalParams,'OD0')
    GlobalParams.OD0 = 6;
end

% dimensionality of wave (bands vs rings). really only 1D is fully
% supported at the moment... step sizes need to be adjusted for 2D to be
% stable
if ~isfield(GlobalParams,'dim')
    GlobalParams.dim   = 1; %1D or 2D
end

% if set to 1, simulation will use a biophysical model for \chi, the
% chemotactic coefficient. otherwise, chi = alpha * mu, where mu is the
% diffusion coefficient, and alpha is an environment-dependent parameter
% but is the same for all phenotypes
if ~isfield(GlobalParams,'biophysicalChi')
    GlobalParams.biophysicalChi = 0;
end

% if biophysicalChi = 0, can choose to overwrite the default alpha = chi/mu
% with this parameter
if ~isfield(GlobalParams,'alpha')
    GlobalParams.alpha = [];
end

if ~isfield(GlobalParams,'growthRate')
    GlobalParams.growthRate = 0;
end

% what physical environment the wave travels through
if ~isfield(GlobalParams,'environment')
    GlobalParams.environment = 'liquid';
elseif ~any(strcmp(GlobalParams.environment,{'liquid','agar'}))
    disp('Error: GlobalParams.environment must be a string, ''liquid'' or ''agar''.')
    return
end

if isfield(GlobalParams,'pore_size')
    if ~isempty(GlobalParams.pore_size)
		if GlobalParams.pore_size<Inf
			GlobalParams.environment = 'agar';
		end
	end
else
	GlobalParams.pore_size = [];
end

% whether to simulate migration through a single environment or
% periodically changeing environments
if ~isfield(GlobalParams,'run_case')
    GlobalParams.run_case = 'single_environment';
elseif ~any(strcmp(GlobalParams.run_case,{'switch_environments','single_environment'}))
    disp('Error: GlobalParams.run_case must be a string, either ''switch_environments'' or ''single_environment''.')
    return
end

% run for a fixed amount of time (per environment) or fixed travel distance
if ~isfield(GlobalParams,'totTime')
    GlobalParams.totTime = 5*60*60; % s
end

% sets the spatial grid size by dividing the characteristic length scale of
% the problem by x_discretize_factor
if ~isfield(GlobalParams,'x_discretize_factor')
    GlobalParams.x_discretize_factor = 50;
elseif isempty(GlobalParams.x_discretize_factor)
    GlobalParams.x_discretize_factor = 50;
end
if GlobalParams.x_discretize_factor~=50
    xmeshName = ['_xmesh' num2str(GlobalParams.x_discretize_factor)];
else
    xmeshName = [];
end
    

% if run_case = 'switch_environment', how many times to alternate
% environments. the code enforces this to be an even number if nEnv>1
if ~isfield(GlobalParams,'nEnv')
    GlobalParams.nEnv = 1;
elseif GlobalParams.nEnv>1
    if ~mod(GlobalParams.nEnv,2)==0
        GlobalParams.nEnv = GlobalParams.nEnv+1;
    end
end

% parameters describing the tumble bias distribution
% NOTE: this is not actually the mean TB. it's the mode of the TB
% distribution.
if ~isfield(GlobalParams,'peakTB')
    GlobalParams.peakTB=0.140; 
end

% this is the standard deviation of the TB distribution.
if ~isfield(GlobalParams,'stdTB')
    GlobalParams.stdTB  = 0.065;
end

% uniform aspartate concentration in the environment
if ~isfield(GlobalParams, 'asp')
    GlobalParams.asp = 100; % uM , total aspartate concentration
end


% carrying capacity (units of OD); set large just to keep cells behind the wave from growing indefinitely
% not needed if the simulation has a nutrient
if GlobalParams.growthRate>0 && ~isfield(GlobalParams,'cap')
    GlobalParams.cap = 100; % max OD; 
end



%% phenotypes model
% phi is the correlation between mother and daughter phenotypes
if ~isfield(GlobalParams,'phi')
    if ~isfield(GlobalParams,'stdTB')
        % default case
        GlobalParams.phi = .7; % 0.65
    elseif GlobalParams.stdTB>0
        GlobalParams.phi = 0.7; % 0.65
    else
        % only case without diversity
        GlobalParams.phi = 0;
    end
elseif isempty(GlobalParams.phi)
    if ~isfield(GlobalParams,'stdTB')
        % default case
        GlobalParams.phi = 0.7; % 0.65
    elseif GlobalParams.stdTB>0
        GlobalParams.phi = 0.7; % 0.65
    else
        % only case without diversity
        GlobalParams.phi = 0;
    end
end

if GlobalParams.phi<0
    GlobalParams.phi = 0;
end


[GlobalParams,SimParams] = initializePhen(GlobalParams,SimParams);

peakTBname = ['peakTB=' num2str(round(GlobalParams.peakTB,3))];
stdTBname = ['stdTB=' num2str(round(GlobalParams.stdTB,3))];


%% attractant parameters
% aspartate parameters
% diffusivity
if ~isfield(GlobalParams,'DA')
    % should default value be 800 or 500?
    GlobalParams.DA  = 800; % um^2/s, aspartate diffusion coefficient, A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry, Volume 166, Issue 2, 1 November 1987, Pages 335-341
    DAname = ['DA=' num2str(GlobalParams.DA) 'um2s'];
elseif GlobalParams.DA==800
    DAname = ['DA=' num2str(GlobalParams.DA) 'um2s'];
else
    DAname = ['DA=' num2str(GlobalParams.DA) 'um2s'];
end
% consumption rate
if ~isfield(GlobalParams,'kA')
    GlobalParams.kA  = 9*0.84/60; % uM/OD/s, asp consumption rate when cell density is low (oxygen is maximal), 1e-12 umol/cell/min = 0.84/60 uM/OD/s
    kAname = [];
elseif GlobalParams.kA  == 9*0.84/60
    kAname = [];
else
    kAname = ['kA=' num2str(GlobalParams.DA) 'uMODs'];
end
% Michaelis-Menten constant of consumption
if ~isfield(GlobalParams,'KmA')
    GlobalParams.KmA = 0.5; % uM
    KmAname = ['KmA=' num2str(GlobalParams.KmA) 'uM'];
elseif isfield(GlobalParams,'KmA') && isempty(GlobalParams.KmA)
    GlobalParams.KmA = 0.5; % uM
    KmAname = ['KmA=' num2str(GlobalParams.KmA) 'uM'];
elseif GlobalParams.KmA == 0.5
    KmAname = ['KmA=' num2str(GlobalParams.KmA) 'uM'];
else
    KmAname = ['KmA=' num2str(GlobalParams.KmA) 'uM'];
end
% receptor dissocation constant in the inactive state
if ~isfield(GlobalParams,'KiA')
    GlobalParams.KiA  = 3.5; % uM % is this number just the half-occupancy ligand concentration from Mello and Tu 2003? not the same as Ki... Sourjik and Berg 2004 Nature estimate Ki=30uM for MeAsp and Ka=1mM for MeAsp. but those are different from Asp.
    KiAname = ['KiA=' num2str(GlobalParams.KiA) 'uM'];
elseif GlobalParams.KiA  == 3.5 % uM
    KiAname = ['KiA=' num2str(GlobalParams.KiA) 'uM'];
else
    KiAname = ['KiA=' num2str(GlobalParams.KiA) 'uM'];
end
    
% GlobalParams.KaA  = 55000; % uM, dissociation constant of asp with chemoreceptor. any good reference for this???
% receptor dissocation constant in the inactive state
if ~isfield(GlobalParams,'KaA')
    GlobalParams.KaA  = 55000; % uM % is this number just the half-occupancy ligand concentration from Mello and Tu 2003? not the same as Ki... Sourjik and Berg 2004 Nature estimate Ki=30uM for MeAsp and Ka=1mM for MeAsp. but those are different from Asp.
    KaAname = [];
elseif GlobalParams.KaA  == 55000 % uM
    KaAname = [];
else
    KaAname = ['KaA=' num2str(GlobalParams.KaA) 'uM'];
end


%% initial condition related
if ~isfield(GlobalParams,'IC')
    GlobalParams.IC = 'step'; % or 'sigmoid'
end

ICname = [GlobalParams.IC 'IC'];

%% set up simulation grid, time
SimParams = setupXT(GlobalParams,SimParams);


%% set up initial condition

if ~isfield(SimParams(end),'rho0') % not sure rho0 is currently used    
    SimParams = initializeDefaultIC(GlobalParams,SimParams);
end


%% names
% dimensionality
switch GlobalParams.dim
    case 1
        dStr = '1D';
    case 2
        dStr = '2D';
end

% top-most folder name
topDir = [dStr, '_', 'Asp', num2str(GlobalParams.asp), 'uM'];

if ~isempty(peakTBname)
    topDir = [topDir, '_', peakTBname];
end
if ~isempty(stdTBname)
    topDir = [topDir, '_', stdTBname];
end

if GlobalParams.biophysicalChi==1
    topDir = [topDir, '_biophysicalChi'];
elseif isfield(GlobalParams,'alpha')
    if ~isempty(GlobalParams.alpha)
        topDir = [topDir, '_alpha=' num2str(GlobalParams.alpha)];
    end
end

if ~isempty(xmeshName)
    topDir = [topDir, xmeshName];
end

% additional cases
subDir = '';
if GlobalParams.growthRate > 0
    
    subDir = ['r=', num2str(round(3600*GlobalParams.growthRate,3)), 'perHr'];
    
    nPhen = length(SimParams.F0);
    if nPhen > 0
        subDir = [subDir,'_phi=', num2str(round(GlobalParams.phi,3))];
    end
    
else
    switch GlobalParams.IC
        case 'sigmoid'
            subDir = [ICname, '_OD', num2str(round(GlobalParams.OD0,3))];
        case 'step'
            subDir = [ICname, '_OD', num2str(round(GlobalParams.OD0,3))];
        
    end
end

if ~isempty(KiAname)
    subDir = [subDir, '_', KiAname];
end
if ~isempty(KaAname)
    subDir = [subDir, '_', KaAname];
end
if ~isempty(DAname)
    subDir = [subDir, '_', DAname];
end
if ~isempty(kAname)
    subDir = [subDir, '_', kAname];
end
if ~isempty(KmAname)
    subDir = [subDir, '_', KmAname];
end



if ~isempty(subDir)
	if subDir(end) == '_'
	    subDir = subDir(1:end-1);
    end
end

% get the directory where this code is located
ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end

if ~strcmp(ks_dir(end),filesep)
    ks_dir = [ks_dir,filesep];
end

[~,ind2] = regexp(ks_dir,[filesep,'ks',filesep]);
ks_dir = ks_dir(1:ind2);

GlobalParams.dataDir = [ks_dir, 'data', filesep, topDir, filesep, subDir];

switch GlobalParams.run_case
    case 'switch_environments'
        if ~isfield(GlobalParams,'k')
            GlobalParams.k = 1;
        end
        
%         GlobalParams.dataDir = [GlobalParams.dataDir, '_nEnv=' num2str(GlobalParams.nEnv)];
        GlobalParams.dataDir = [GlobalParams.dataDir, '_switch_envs'];
        
        GlobalParams.dataDir = [GlobalParams.dataDir, '_rTEnv=' num2str(round(GlobalParams.totTime*GlobalParams.growthRate,2))];
        
        GlobalParams.simName = [GlobalParams.dataDir, filesep, 'results_env_' sprintf('%03d',GlobalParams.k) '.mat'];
    
    case 'single_environment'
        GlobalParams.dataDir = [GlobalParams.dataDir, '_', GlobalParams.environment];
        GlobalParams.simName = [GlobalParams.dataDir, filesep, 'results.mat'];
        % may need to extend ss simulations. should be easy to extend
        % these.
end

end