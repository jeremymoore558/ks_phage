function [SimResults, SimParams] = simulateWave(GlobalParams, SimParams, varargin)

%% initialization
warning('off','signal:findpeaks:largeMinPeakHeight');
dense_factor = 10; % interpolate profile to get peak. increase sampling density by this factor

% check for output of previous simulation to start from
if nargin>=3 && ~isempty(varargin{1})
    SimResults = varargin{1};
end

% grid sizes
dx = SimParams(end).dx;
dt = SimParams(end).dt;

% channel set up
L = SimParams(end).L; % number of spatial index
x = SimParams(end).x; % um, position along channel, or radial position in 2D
nPhen = length(SimParams(end).F0);

% cell
MFTmodel = initializeMFTmodel(SimParams(end).F0, GlobalParams.environment, GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
Mu = MFTmodel.mu(:)*ones(1,L); % cell diffusion coefficients; matrix of size (nPhen x L)
Chi = MFTmodel.chi(:)*ones(1,L); % chemotactic diffusion coefficients; matrix of size (nPhen x L)

% initial conditions
rho = SimParams(end).rho0;

% boundary conditions
rho = cat(2,rho(:,2,:),rho,rho(:,end-1,:));

% aspartate parameters
DA   = GlobalParams.DA; % um^2/s, aspartate diffusion coefficient, A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry, Volume 166, Issue 2, 1 November 1987, Pages 335-341
kA   = GlobalParams.kA; % uM/OD/s, asp consumption rate when cell density is low
KmA   = GlobalParams.KmA; % uM, Michaelis-Menten constant in asp consumption

A = SimParams(end).asp0;

% boundary conditions
A = cat(2,A(2),A,A(end-1));

%% with growth
if GlobalParams.growthRate > 0
    growthMatrix = computeGrowthMatrix(GlobalParams,SimParams);
end


%% simulation steps
% initialize temporal profiles
rhot = zeros(nPhen, L, 1);
At = zeros(L, 1);


peakLocs = nan(1,1);

rho_tot = squeeze(sum(rho(:,2:end-1,:),1));
[pks,locs] = findpeaks(rho_tot,'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
if ~isempty(locs)
    startingPeakRho = pks(end);
else
    startingPeakRho = max(rho_tot(:,1));
end

backInd = find(A<1e-3,1,'last'); % important to make sure this is consistent across functions. <...'last or >...'first'
if isempty(backInd)
    backInd = 1;
end
NF0 = squeeze(sum(rho(:,backInd:end-1,1),2));
PF0_prev = NF0/sum(NF0,1);

tt        = 0; % output time index
if length(SimParams)==1
    nextOutTime = SimParams.time(end)+SimParams.outDt;
else
    nextOutTime = SimParams(end-1).time(end)+SimParams(end).outDt;
end

t1 = tic;
t = 0;
continueSimulating = 1;
step = length(SimParams);

while continueSimulating
    t = t+1;
    % Runge-Kutta RK4, https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    
    % Intermediate variable 1
    rho1 = rho;
    rho1(rho1<0) = 0;
    A1 = A;
    A1(A1<0) = 0;
    
    [drhodt_a, dAdt_a] = integrateStep(rho1, A1);
    
    % Intermediate variable 2
    rho1 = rho + drhodt_a*dt/2;
    rho1(rho1<0) = 0;
    A1 = A + dAdt_a*dt/2;
    A1(A1<0) = 0;
    
    [drhodt_b, dAdt_b] = integrateStep(rho1, A1);
    
    % Intermediate variable 3
    rho1 = rho + drhodt_b*dt/2;
    rho1(rho1<0) = 0;
    A1 = A + dAdt_b*dt/2;
    A1(A1<0) = 0;
    
    [drhodt_c, dAdt_c] = integrateStep(rho1, A1);
    
    % Intermediate variable 4
    rho1 = rho + drhodt_c*dt;
    rho1(rho1<0) = 0;
    A1 = A + dAdt_c*dt;
    A1(A1<0) = 0;
    
    [drhodt_d, dAdt_d] = integrateStep(rho1, A1);
    
    % Combine
    rho = rho + (drhodt_a + 2*drhodt_b + 2*drhodt_c + drhodt_d)*dt/6;
    A = A + (dAdt_a + 2*dAdt_b + 2*dAdt_c + dAdt_d)*dt/6;
    rho(:,1,:) = rho(:,2,:);
    rho(:,end,:) = rho(:,end-1,:);
    rho(rho<0) = 0;
    A(1) = A(2);
    A(end) = A(end-1);
    A(A<0) = 0;
    
    
    if any(isnan(rho(:))) || any(isnan(A(:)))
        disp(['Error: nans detected. Time step = ' num2str(t) '.'])
        
        break
    end
    
    % check for wave location
    backInd = find(A(2:end-1)<1e-3,1,'last');
    if isempty(backInd)
        backInd = 1;
    end
    rho_tot = squeeze(sum(rho(:,2:end-1,:),1));
    [pks,locs] = findpeaks(rho_tot(:,backInd:end,:),'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
    if ~isempty(pks)
        pks = pks(end);
        locs = locs(end);
        locs = locs+backInd-1;
    else
        pks = startingPeakRho;
    end
    
    % for checking the edges of the F0 distribution
    if t>1
        PF0_prev = PF0;
    end
    backInd = find(A<1e-3,1,'last'); % important to make sure this is consistent across functions. <...'last or >...'first'
    NF0 = squeeze(sum(rho(:,backInd:end-1,1),2));
    PF0 = NF0/sum(NF0,1);
    
    % output profiles at regular time intervals
    % absolute time is (t-1)*dt
    if step==1
        absTime = (t-1)*dt;
    else
        absTime = SimParams(end-1).time(end)+(t-1)*dt;
    end
    
    if absTime>=nextOutTime || (t==1 && step==1)
        disp(['Simulation time = ' num2str(round(absTime/60,2)) ' min.']) % of expected ' num2str(round(GlobalParams.totTime/60,2)) ' min.'])
        
        tt = tt+1;
        
        % record time point
        rhot(:,:,tt,:) = rho(:,2:end-1,:);
        
        At(:,tt)   = A(2:end-1);
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        rho_dense = interp1(x,rho_tot,x_dense,'pchip');
        [~,locs2] = findpeaks(rho_dense,'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
        
        if ~isempty(locs2)
            peakLocs(tt) = x_dense(locs2(end));
            disp(['Wave peak position = ' num2str(peakLocs(tt)/1000) ' mm.'])
        else
            peakLocs(tt) = nan;
        end
        
        SimParams(end).time(tt) = absTime;
        
        if ~(t==1 && step==1)
            nextOutTime = nextOutTime+SimParams(end).outDt;
        end
        
        toc(t1)
        
    end
    
    % need to know when to stop iterating: peak of wave has reached x >=
    % totLength, or simulation time>totTime
    if absTime>=GlobalParams.totTime
        continueSimulating = 0;
        
        disp(['Simulation time = ' num2str(round(absTime/60,2)) ' min.'])% of expected ' num2str(round(GlobalParams.totTime/60,2)) ' min.'])
        
        tt = tt+1;
        
        % record time point
        rhot(:,:,tt,:) = rho(:,2:end-1,:);
        
        At(:,tt)   = A(2:end-1);
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        rho_dense = interp1(x,rho_tot,x_dense,'pchip');
        [~,locs2] = findpeaks(rho_dense,'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
        
        if ~isempty(locs2)
            peakLocs(tt) = x_dense(locs2(end));
            disp(['Wave peak position = ' num2str(peakLocs(tt)/1000) ' mm.'])
        else
            peakLocs(tt) = nan;
        end
        
        SimParams(end).time(tt) = absTime;
        
        % collect results from this step
        SimResults(step).rho            = rhot; % nPhen * x * GlobalParams.time (* 2) array of cell density in OD
        SimResults(step).asp            = At; % x * GlobalParams.time array of asp concentration in uM
        
        SimResults(step).peakLocs = peakLocs;
        
        toc(t1)
        
        % save current progress
        disp('Saving current progress.')
        save(GlobalParams.simName, 'SimParams', 'SimResults', 'GlobalParams', '-v6')
        matObj = matfile(GlobalParams.simName);
        vars = who(matObj);
        
        % check whether file save. if not, use v7.3
        if all(cellfun(@isempty, regexp(vars, 'SimResults')))
            disp('Saving variables in .mat -v7.3 format.')
            save(GlobalParams.simName, 'SimParams', 'SimResults', 'GlobalParams', '-v7.3')
        end
        
        break
    end
    
    
    
    % check whether to start next block of simulation
    aspCond = A(end-2)<(1-1e-6)*GlobalParams.asp;
    peakCond = abs(log(pks/startingPeakRho))>0.15;
    lowF0Cond = (PF0(1)/max(PF0)>1e-4 && (PF0(1)-PF0_prev(1))>0 && GlobalParams.growthRate>0);
    highF0Cond = (PF0(end)/max(PF0)>1e-4 && (PF0(end)-PF0_prev(end))>0 && GlobalParams.growthRate>0);
    
    if (aspCond || peakCond || lowF0Cond || highF0Cond) && continueSimulating
        if aspCond
            flag = 1;
            disp(['Flag = ' num2str(flag) '; aspartate at the right boundary has dropped below threshold.']);
        elseif peakCond
            flag = 2;
            disp(['Flag = ' num2str(flag) '; peak cell density has changed significantly.']);
        elseif lowF0Cond
            flag = 3;
            disp(['Flag = ' num2str(flag) '; need to extend phenotype grid to include lower values of F0.']);
        elseif highF0Cond
            flag = 4;
            disp(['Flag = ' num2str(flag) '; need to extend phenotype grid to include higher values of F0.']);
        end
        
        % store last point of this block
        disp(['Simulation time = ' num2str(round(absTime/60,2)) ' min.'])% of expected ' num2str(round(GlobalParams.totTime/60,2)) ' min.'])
        
        tt = tt+1;
        
        % record time point
        rhot(:,:,tt,:) = rho(:,2:end-1,:);
        
        At(:,tt)   = A(2:end-1);
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        rho_dense = interp1(x,rho_tot,x_dense,'pchip');
        [~,locs2] = findpeaks(rho_dense,'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
        
        if ~isempty(locs2)
            peakLocs(tt) = x_dense(locs2(end));
            disp(['Wave peak position = ' num2str(peakLocs(tt)/1000) ' mm.'])
        else
            peakLocs(tt) = nan;
        end
        
        SimParams(end).time(tt) = absTime;
        
        % collect results from this step
        SimResults(step).rho            = rhot; % nPhen * x * GlobalParams.time (* 2) array of cell density in OD
        SimResults(step).asp            = At; % x * GlobalParams.time array of asp concentration in uM
        
        SimResults(step).peakLocs = peakLocs;
        
        toc(t1)
        
        % save current progress
        disp('Saving current progress.')
        save(GlobalParams.simName, 'SimParams', 'SimResults', 'GlobalParams', '-v6')
        matObj = matfile(GlobalParams.simName);
        vars = who(matObj);
        
        % check whether file save. if not, use v7.3
        if all(cellfun(@isempty, regexp(vars, 'SimResults')))
            save(GlobalParams.simName, 'SimParams', 'SimResults', 'GlobalParams', '-v7.3')
        end
        
        
        % set up next phase
        disp('Beginning next simulation block.')
        [SimParams,GlobalParams,flag] = setupXT(GlobalParams,SimParams,SimResults);
        
        if flag==0
            break
            
        else
            step = step+1; % block number
            tt = 0; % storage index
            t = 1; % time step index
            dt = SimParams(end).dt;
            dx = SimParams(end).dx;
            nextOutTime = SimParams(end-1).time(end) + SimParams(end).outDt;
            peakLocs = [];
            
            x = SimParams(end).x;
            L = SimParams(end).L;
            
            MFTmodel = initializeMFTmodel(SimParams(end).F0, GlobalParams.environment, GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
            Mu = MFTmodel.mu(:)*ones(1,L); % cell diffusion coefficients; matrix of size (nPhen x L)
            Chi = MFTmodel.chi(:)*ones(1,L); % chemotactic diffusion coefficients; matrix of size (nPhen x L)
            
            if GlobalParams.growthRate > 0
                growthMatrix = computeGrowthMatrix(GlobalParams,SimParams);
            end
            
            nPhen = length(SimParams(end).F0);
            
            % set up initial condition to continue from
            rho = zeros(nPhen,L+2);
            for i = 1:nPhen
                rho(i,2:end-1,:) = SimParams(end).rho0(i,:);
            end
            rho(:,1,:)=rho(:,2,:);
            rho(:,end,:)=rho(:,end-1,:);
            
            rhot = zeros(nPhen, L, 1);
            
            A = zeros(1,L+2);
            A(2:end-1) = SimParams(end).asp0;
            A(1)=A(2);
            A(end)=A(end-1);
            At = zeros(L, 1);
            
            
            backInd = find(A(2:end-1)<1e-3,1,'last');
            if isempty(backInd)
                backInd = 1;
            end
            rho_tot = squeeze(sum(rho(:,2:end-1,:),1));
            [pks,~] = findpeaks(rho_tot(:,backInd:end,:),'MinPeakHeight',0.1,'MinPeakProminence',0.01);%0.5
            if ~isempty(pks)
                startingPeakRho = pks(end);
            else
                startingPeakRho = max(rho_tot(:,1));
            end
            
            backInd = find(A<1e-3,1,'last'); % important to make sure this is consistent across functions. <...'last or >...'first'
            if isempty(backInd)
                backInd = 1;
            end
            NF0 = squeeze(sum(rho(:,backInd:end-1,1),2));
            PF0_prev = NF0/sum(NF0,1);
            
        end
    end
    
    
end


warning('on','signal:findpeaks:largeMinPeakHeight');




%% functions
    function fS = receptorSensing(L, Ki, Ka) % convert to perceived signal
        fS = log((1+L/Ki)./(1+L/Ka));
    end

    function [drhodt, dAdt] = integrateStep(rho1, A1) % atomic integration step
        % initialize
        drhodt = zeros(size(rho1));
        dAdt = zeros(size(A1));
        
        % cell density dynamics
        
        % calculate perceived signal
        f = receptorSensing(A1, GlobalParams.KiA, GlobalParams.KaA);
        
        % calculate chemotactic drift
        g = zeros(size(rho1)); % initialize drift term as temporary variable
        g(:,2:end-1,1) = Chi .* repmat((f(3:end) - f(1:end-2))./(2*dx),size(Chi,1),1); % chemotactic drift
        
        
        % chemotactic flux
        g = g .* rho1;
        g(A1==0) = 0;
        
        % divergence of all chemotactic and diffusive fluxes
        drhodt(:,2:end-1,:) = Mu .* (rho1(:,1:end-2,:) + rho1(:,3:end,:) - 2*rho1(:,2:end-1,:))./(dx).^2 ... % diffusion coefficient will copy into the third dimension
            - (g(:,3:end,:) - g(:,1:end-2,:))./(2*dx); % diffusion term
        
        % modify divergence in 2D
        if GlobalParams.dim == 2 % radius will copy from the second dimension into other dimensions
            drhodt(:,2:end-1,:) = drhodt(:,2:end-1,:) ...
                + 1./(x+dx/2) .* Mu(:,2:end-1) .* (rho1(:,3:end,:) - rho1(:,1:end-2,:))./(2*dx) ...
                - 1./(x+dx/2) .* g(:,2:end-1,:);
        end
        
        % add cell growth
        if GlobalParams.growthRate > 0
            drhodt(:,2:end-1,:) = drhodt(:,2:end-1,:) + (1-repmat(sum(rho1(:,2:end-1,:),1),nPhen,1)/GlobalParams.cap) .* (growthMatrix * rho1(:,2:end-1,:)); % divisions remove phenotype and produce two identical cells; and can produce self
        end
        
        % aspartate concentration dynamics
        maximalConsumption = squeeze(sum(rho1,1)) * kA;
        tempA = A1./ (A1+KmA);
        tempA(A1==0) = 0; % allows KmA = 0
        
        GammaA = maximalConsumption(:)' .* tempA; % uM/s, aspartate consumption term
        dAdt(2:end-1) = DA * (A1(1:end-2) + A1(3:end) - 2*A1(2:end-1))./(dx).^2 - GammaA(2:end-1);
        if GlobalParams.dim == 2
            dAdt(2:end-1) = dAdt(2:end-1) + 1./x * DA .* (A1(3:end) - A1(1:end-2))./(2*dx);
        end
        
        
    end


end



