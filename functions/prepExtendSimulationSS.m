function [GlobalParams, SimParams, flag] = prepExtendSimulationSS(GlobalParams, SimParams, SimResults)

flag = 1;
    
% set up space, time
[SimParams,GlobalParams,flag] = setupXT(GlobalParams,SimParams,SimResults);

if flag == 0 
    disp('Error setting up next simulation.')
    return
end

end

