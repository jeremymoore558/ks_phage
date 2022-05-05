function [GlobalParams,SimParams] = initializePhen(GlobalParams,SimParams)

% default values
% set mean and variance of TB, even if phenotype = 'F0'. then choose mean
% and variance of F0 if needed so that TB comes out with the right mean and
% variance.

peakTB = GlobalParams.peakTB; % really the peak TB
stdTB = GlobalParams.stdTB; % really the RMS variation around the peak....
% maximum value possible for stdTB is sqrt(1/12) -- uniform
% distribution of TB
    
if ~isfield(SimParams,'F0')
    
    if stdTB>sqrt(1/12)
        stdTB=sqrt(1/12);
        disp('stdTB set to maximum possible value.')
    end

    if stdTB>0
        % compute meanF0 and stdF0 from peakTB and stdTB
        options = optimoptions('fsolve','Display','off');
        F0params = fsolve(@(F0params) resFun(F0params,[peakTB,stdTB]), [log((1-peakTB)/peakTB);log(1)], options);
        
        meanF0 = F0params(1); % WT: 1.182; default: 1.9551
        stdF0 = exp(F0params(2)); % WT: 0.4583; default: 0.46
    else
        meanF0 = log((1-peakTB)/peakTB);
        stdF0 = 0;
    end
    
    GlobalParams.meanF0 = meanF0;
    GlobalParams.stdF0 = stdF0;
    
    if stdTB==0
        F0 = meanF0;
    else
        if GlobalParams.growthRate == 0
            F0 = (meanF0-6*stdF0):stdF0/10:(meanF0+6*stdF0); % similar range as TB case. could be 150?
        else
            F0 = (meanF0-6*stdF0):(stdF0*sqrt(1-GlobalParams.phi.^2))/10:(meanF0+6*stdF0); % required resolution depends on growth matrix
        end
        F0 = F0(:);
    end
    
    SimParams.F0 = F0(:);
    
else
    
    F0 = SimParams.F0(:);
    PF0 = SimParams.PF0(:);
    meanF0 = mean(F0);
    stdF0 = std(F0);
    
    GlobalParams.meanF0 = meanF0;
    GlobalParams.stdF0 = stdF0;

    TB = 1./(1+exp(F0));
    
    GlobalParams.peakTB = sum(PF0.*TB)/sum(PF0);
    GlobalParams.stdTB = sqrt(sum(PF0.*(TB-GlobalParams.peakTB).^2)/sum(PF0));
    
end

GlobalParams.F0 = F0(:);

if ~isfield(SimParams,'PF0')
    if stdTB>0
        GlobalParams.PF0 = normpdf(F0(:),meanF0,stdF0);
        GlobalParams.PF0 = GlobalParams.PF0(:)/sum(GlobalParams.PF0);
    else
        GlobalParams.PF0 = 1;
    end
else
    GlobalParams.PF0 = PF0(:);
end

end

function res = resFun(F0params,TBparams)

peakTB = TBparams(1);
stdTB = TBparams(2);
meanF0 = F0params(1);
stdF0 = exp(F0params(2));

F0 = linspace(meanF0-6*stdF0,meanF0+6*stdF0,10000);
F0 = F0(:);
PF0 = normpdf(F0(:),meanF0,stdF0);
PF0 = PF0(:);
% TB = 1./(1+exp(F0));
[TB,PTB] = convertPF0ToPTB(F0,PF0); % this is the problem, it seems
[~,mxi] = max(PTB);
peakTB_est = TB(mxi);
stdTB_est = sqrt(trapz(F0,(TB-peakTB_est).^2.*PF0)/trapz(F0,PF0));

res = [peakTB-peakTB_est;stdTB-stdTB_est];

end