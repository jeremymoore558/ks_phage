function growthMatrix = computeGrowthMatrix(GlobalParams,SimParams)

nPhen = length(SimParams(end).F0);

if nPhen == 1
    growthMatrix = GlobalParams.growthRate;
end

if GlobalParams.growthRate>0 && nPhen > 1
    
    growthMatrix = zeros(nPhen);
    
    if ~isfield(GlobalParams,'phi')
        disp('Error: parameter tau in the growth model not specified.')
        return
    elseif isempty(GlobalParams.phi)
        disp('Error: parameter tau in the growth model not specified.')
        return
    end
    
    r = GlobalParams.growthRate;
    phi = GlobalParams.phi; % exp(-1/tau) is mother-daughter correlation
    meanF0 = GlobalParams.meanF0;
    stdF0 = GlobalParams.stdF0;
    
    for i = 1:nPhen
        growthMatrix(:,i) = normpdf(SimParams(end).F0, meanF0+(SimParams(end).F0(i)-meanF0)*phi, stdF0*sqrt((1-phi.^2)));
        growthMatrix(:,i) = growthMatrix(:,i)/sum(growthMatrix(:,i));
    end
    
    growthMatrix = r*(2*growthMatrix - eye(nPhen));
    
end

end

