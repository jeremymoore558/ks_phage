function [logZ,NF0_hat] = findLogZ(F0,NF0,dx,GlobalParams,makePlots,varargin)
%logZ = findLogZ(F0,NF0,GlobalParams,dx,varargin)

F0 = F0(:);

NF0(NF0<1e-10) = 0;
Ntot = sum(NF0);
PF0 = NF0/sum(NF0);
c = Ntot*GlobalParams.kA/GlobalParams.asp;

MFTmodel = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);
chi = MFTmodel.chi;
mu = MFTmodel.mu;

meanChi = sum(PF0.*chi);
meanMu = sum(PF0.*mu);
lambda = c/(meanChi-meanMu);

[~,mxi] = max(chi);

z = -8/lambda:(dx):3/lambda;

LB = -inf*ones(size(F0));
UB = inf*ones(size(F0));
options = optimoptions('fmincon','MaxIterations',1e5,'MaxFunctionEvaluations',1e5,'FiniteDifferenceType','central','OptimalityTolerance',1e-3,'StepTolerance',1e-3,'FunctionTolerance',1e-3);%,'Display','none');

temp = -eye(length(F0)) + diag(ones(length(F0)-1,1),1);
A = repmat(sign(diff(chi(:))),1,length(F0)).*temp(1:end-1,:); % A = temp(1:end-1,1);
B = zeros(length(F0)-1,1);

Aeq = zeros(1,length(F0)); Aeq(mxi) = 1; % set the phenotype with maximum chi to have logZ = 0
Beq = 0;

if isempty(varargin)
    logZ0 = zeros(size(F0));
    
    tic
    [logZ,resnorm,residuals] = fmincon(@(logZ) sum(computeRes(logZ,F0,NF0,GlobalParams,z).^2),logZ0,A,B,Aeq,Beq,LB,UB,[],options);
    toc
    
    NF0_hat = computeIntegral(logZ,F0,NF0,GlobalParams,z);
    
    if makePlots
    figure;hold on
    plot(F0,NF0)
    plot(F0,NF0_hat,'--')
    xlabel('F0');ylabel('N(F0)')
    end
    
else
    logZ0 = varargin{1};
    
    tic
    
    [logZ,resnorm,residuals] = fmincon(@(logZ) sum(computeRes(logZ,F0,NF0,GlobalParams,z).^2),logZ0,A,B,Aeq,Beq,LB,UB,[],options);
    
    toc
        
    NF0_hat = computeIntegral(logZ,F0,NF0,GlobalParams,z);
    
    if makePlots
    figure;hold on
    plot(F0,NF0)
    plot(F0,NF0_hat,'--')
    xlabel('F0');ylabel('N(F0)')
    end
end

    function residuals = computeRes(logZ,F0,NF0,GlobalParams,z)
        
        integrals = computeIntegral(logZ,F0,NF0,GlobalParams,z);
        residuals = (integrals(:)) - (NF0(:));
        residuals(NF0<1e-10) = 100;
        residuals(integrals<1e-10) = 100;
        
        residuals = [residuals; diff(logZ(:)).^2];
    end

    function integrals = computeIntegral(logZ,F0,NF0,GlobalParams,z)
        
        dz = z(2)-z(1);
        
        [rho0,s0] = KSpredict(F0,NF0,GlobalParams,dz,logZ,z);

        indz0 = 2;
        integrals = nansum(rho0(:,indz0-1:end),2)*dz;
        
    end
end