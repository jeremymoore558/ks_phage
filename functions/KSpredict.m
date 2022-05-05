function [rho0,s0,z,logZ] = KSpredict(F0,NF0,GlobalParams,dx,varargin)
% dx ensures grid size, and therefore integrals over position, are done the
% same way as what was done to produce NF0

F0 = F0(:);
NF0 = NF0(:);

% dF0 = F0(2)-F0(1);

Ntot = sum(NF0(:));
PF0 = NF0/sum(NF0);
asp0 = GlobalParams.asp;
c = Ntot*GlobalParams.kA/asp0;
ep = GlobalParams.KiA/asp0;

MFTmodel = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);

alpha = sum(PF0.*MFTmodel.chi./MFTmodel.mu);

meanChi = sum(PF0.*MFTmodel.chi);
meanMu = sum(PF0.*MFTmodel.mu);
lambda = c/(meanChi-meanMu);

chi = MFTmodel.chi;
mu = MFTmodel.mu;

% z = linspace(-6/lambda,3/lambda,10^4);
z = -8/lambda:dx:3/lambda;


if nargin>=5
    if ~isempty(varargin{1})
        logZ = varargin{1};
    end
else
    logZ = findLogZ(F0,NF0,GlobalParams,dx);
end

if nargin>=6
    if ~isempty(varargin{2})
        z = varargin{2};
        dx = z(2)-z(1);
    end
end

s0 = asp0*((1+ep)*1./(1+(alpha-1)*sum(repmat(PF0,1,length(z)).*exp(-c*repmat(1./mu,1,length(z)).*repmat(z,length(mu),1)-repmat(logZ,1,length(z))),1)).^(1/(alpha-1)) - ep);
rho0 = (1+ep)*Ntot.*c.*repmat(PF0./(chi-mu),1,length(z)).*(alpha-1).*exp(-c*repmat(1./mu,1,length(z)).*repmat(z,length(mu),1)-repmat(logZ,1,length(z)))./(1+(alpha-1)*sum(repmat(PF0,1,length(z)).*exp(-c*repmat(1./mu,1,length(z)).*repmat(z,length(mu),1)-repmat(logZ,1,length(z))),1)).^(alpha/(alpha-1));

s0(s0<=0) = 0;
rho0(:,s0==0)=0;


end

