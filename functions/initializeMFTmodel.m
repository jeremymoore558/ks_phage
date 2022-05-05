function MFTmodel = initializeMFTmodel(F0,environment,biophysicalChi,varargin)
% mean field theory model for chi. also model for mu.

F0 = F0(:);

% set default
alpha = [];
pore_size = [];
PF0 = [];

if nargin>=4
    alpha = varargin{1};
end
if nargin>=5
    pore_size = varargin{2};
end
if nargin>=6
    PF0 = varargin{3};
end

% motor parameters
MFTmodel.g1    = 40; % energy in units of kb*T, related to depth of potential wells between run/tumble states of the motor
MFTmodel.g2    = 40; % energy in units of kb*T, related to depth of potential wells between run/tumble states of the motor
MFTmodel.Kd    = 3.06; % uM, binding affinity of CheYp to FliM in the motor
MFTmodel.w     = 1.3; % 1/s, motor switching frequency, assuming 1 motor
MFTmodel.Yp_per_a = 6; % uM, conversion from fraction of receptors active to CheYp concentration
MFTmodel.N     = 6; % receptor cooperativity

% pathway and motor quantities
MFTmodel.F0       = F0(:); % log tumble odds -- NOT RECEPTOR FREE ENERGY
MFTmodel.TB     = 1./(1+exp(F0)); % this defines F0
MFTmodel.deltaG  = MFTmodel.F0/2; %MFTmodel.g1/4-MFTmodel.g2/2*(MFTmodel.Yp./(MFTmodel.Kd+MFTmodel.Yp)); % energy in units of kb*T, steady state free energy difference between motor run/tumble states. % TB = 1/(1+exp(2G))
MFTmodel.Yp      = (MFTmodel.g1*MFTmodel.Kd - 4*MFTmodel.Kd.*MFTmodel.deltaG)./(4.*MFTmodel.deltaG - MFTmodel.g1 + 2*MFTmodel.g2); % uM, steady state concentration of CheYp
MFTmodel.a0       = MFTmodel.Yp/MFTmodel.Yp_per_a; % steady state fraction of receptors active
MFTmodel.f0       = log((1-MFTmodel.a0)./MFTmodel.a0); % RECEPTOR free energy

MFTmodel.lambdaT0 = MFTmodel.w*exp(MFTmodel.deltaG); % 1/s, steady state rate of tumble termination
MFTmodel.lambdaR0 = MFTmodel.w*exp(-MFTmodel.deltaG); % 1/s, steady state rate of run termination

MFTmodel.dF0dTB = -1./(MFTmodel.TB.*(1-MFTmodel.TB));% if F0 were receptor free energy... % 4*MFTmodel.g2*MFTmodel.Yp_per_a./(MFTmodel.TB.*(1-MFTmodel.TB).*(MFTmodel.g1-2*log((1-MFTmodel.TB)./MFTmodel.TB)).*(-2*MFTmodel.g2*MFTmodel.Yp_per_a+MFTmodel.g1.*(MFTmodel.Kd+MFTmodel.Yp_per_a)-2.*(MFTmodel.Kd+MFTmodel.Yp_per_a).*log((1-MFTmodel.TB)./MFTmodel.TB)));
MFTmodel.dTBdF0 = -MFTmodel.TB.*(1-MFTmodel.TB);% if F0 were receptor free energy... % -(exp(MFTmodel.F0+MFTmodel.g2/2+MFTmodel.g2.*MFTmodel.Yp_per_a./(MFTmodel.Kd+exp(MFTmodel.F0)*MFTmodel.Kd+MFTmodel.Yp_per_a)).*MFTmodel.g2.*MFTmodel.Kd.*MFTmodel.Yp_per_a)./( (exp(MFTmodel.g1/2)+exp(MFTmodel.g2.*MFTmodel.Yp_per_a./(MFTmodel.Kd+exp(MFTmodel.F0).*MFTmodel.Kd+MFTmodel.Yp_per_a))).^2 .* (MFTmodel.Kd+exp(MFTmodel.F0).*MFTmodel.Kd+MFTmodel.Yp_per_a).^2 );

% physical parameters
MFTmodel.Drot  = 0.062; %0.062; % rad^2/s, rotational diffusion coefficient
MFTmodel.v0     = 25; %25*4/pi; % um/s, multiply by 4/pi to convert run speed of 2D projections to 3D run speed
MFTmodel.Theta = 0.33; % directional persistence, <cos(theta)>, where theta is the change in angle upon a tumble
MFTmodel.d     = 3; % # of dimensions

if isempty(pore_size)
    switch environment
        case 'agar'
            MFTmodel.pore_size = 30; %23; % um, typical size of a pore in the agar
        case 'liquid'
            MFTmodel.pore_size = inf;
        
    end
else
    MFTmodel.pore_size = pore_size;
end

x = MFTmodel.pore_size*MFTmodel.lambdaR0/MFTmodel.v0;

lambdaC = ((MFTmodel.d-1)*MFTmodel.Drot+(1-MFTmodel.Theta).*MFTmodel.lambdaR0); % direction decorrelation rate
% biophysical chi depends on how log(lambdaC) changes due to changes in
% RECEPTOR free energy in adapted state: d log(lambdaC) / d f0

% diffusion coefficient mu
switch environment
    case 'liquid'
        MFTmodel.mu = MFTmodel.v0^2*(1-MFTmodel.TB)./MFTmodel.d./lambdaC;
        
        if biophysicalChi
            MFTmodel.tau_m     = 20*exp(-6*MFTmodel.TB); % adaptation time. Xiongfei's originally tuned model
            
            dlog_lambdaC_dlambdaR0 = (1-MFTmodel.Theta)./lambdaC;
            dlambdaR0_dG = -MFTmodel.lambdaR0;
            dGdYp = -1/2*MFTmodel.g2.*MFTmodel.Kd./(MFTmodel.Kd+MFTmodel.Yp).^2;
            dYpda0 = MFTmodel.Yp_per_a;
            da0df0 = -MFTmodel.a0.*(1-MFTmodel.a0);
            
            dlog_lambdaC_df0 = dlog_lambdaC_dlambdaR0.*dlambdaR0_dG.*dGdYp.*dYpda0.*da0df0;
            
            MFTmodel.chi = MFTmodel.N .* MFTmodel.tau_m./(MFTmodel.tau_m+1./lambdaC) .* (-dlog_lambdaC_df0) .* MFTmodel.mu;
            
            % for transforming distributions of F0 to those of chi and mu
            MFTmodel.dmudF0 = [];
            MFTmodel.dchidF0 = [];
            
        else
            if isempty(alpha)
                alpha = 10; % roughly mean value of biophysical model chi/mu
                
                %MFTmodel.alpha = alpha;
                
                % prefactor gets the mean value of alpha right for base case
                MFTmodel.chi = 1.3537*(1-MFTmodel.Theta).*MFTmodel.lambdaR0./((1-MFTmodel.Theta).*MFTmodel.lambdaR0 + 2*MFTmodel.Drot) .* alpha .* MFTmodel.mu;
                
                if ~isempty(PF0)
                    MFTmodel.chi = MFTmodel.chi*alpha/sum(PF0(:).*MFTmodel.chi(:)./MFTmodel.mu(:));
                end
                
            else
                MFTmodel.chi = alpha * MFTmodel.mu;
            end
            
            % for transforming distributions of F0 to those of chi and mu
            MFTmodel.dmudF0 = [];
            MFTmodel.dchidF0 = [];
            
        end
        
        
    case 'agar'
        MFTmodel.mu_liquid = MFTmodel.v0^2*(1-MFTmodel.TB)./MFTmodel.d./lambdaC;
        
        % see Licata et al 2016 Biophysical Journal
        MFTmodel.mu = MFTmodel.mu_liquid.*(1-(1+x).*exp(-x)+MFTmodel.Theta.*exp(-x).*(exp(-x)-1+x));
        
        if isempty(alpha)
            alpha = 6.25;
            
            % prefactor gets the mean value of alpha right for base case
            MFTmodel.chi = 1.3538*(1-MFTmodel.Theta).*MFTmodel.lambdaR0./((1-MFTmodel.Theta).*MFTmodel.lambdaR0 + 2*MFTmodel.Drot) .* alpha .* MFTmodel.mu;
            
            if ~isempty(PF0)
                MFTmodel.chi = MFTmodel.chi*alpha/sum(PF0(:).*MFTmodel.chi(:)./MFTmodel.mu(:));
            end
            
        else
            
            MFTmodel.chi = alpha * MFTmodel.mu;
            
        end
        
        
end



end