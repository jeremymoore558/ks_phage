function [GlobalParams,SimParams, flag] = prepExtendSimulationEnvironments(GlobalParams, SimParams, SimResults)
% replace Names
flag = 1;
firstInEnv = true;

% grab results from last block of simulation in previous environment
SimParams = SimParams(end);
SimResults = SimResults(end);

% change environment
switch GlobalParams.environment
    case 'liquid'
        GlobalParams.environment = 'agar';
    case 'agar'
        GlobalParams.environment = 'liquid';
end

% set up space, time
[SimParams,GlobalParams,flag] = setupXT(GlobalParams,SimParams,SimResults,firstInEnv);

SimParams = SimParams(end);
SimParams.time = 0;
SimParams.firstInEnv = true;

if flag == 0 
    disp('Error setting up next simulation.')
    return
end

% new name for the simulation same folder
[~,sind2] = regexp(GlobalParams.simName,'results_env_'); % results
GlobalParams.simName = [GlobalParams.simName(1:sind2), sprintf('%03d', GlobalParams.k), '.mat'];

end

