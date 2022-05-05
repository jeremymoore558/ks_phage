function [SimParams,GlobalParams,flag] = setupXT(GlobalParams,SimParams,varargin)
% optional input: varargin{1} = SimResults, the output of a simulation

flag = 1;
x_discretize_factor = 50;
if isfield(GlobalParams,'x_discretize_factor')
    x_discretize_factor = GlobalParams.x_discretize_factor;
end
t_discretize_factor = 4;


if nargin==2
    % set up initial state
    
    % first run starts from initial condition. 
    
    % compute MFTmodel
    MFTmodel = initializeMFTmodel(SimParams.F0, GlobalParams.environment, GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
    
    % average value of chi-mu
    SimParams.chiScale = sum(GlobalParams.PF0(:).*(MFTmodel.chi(:)-MFTmodel.mu(:))); % um^2/s % (alpha-1)*<mu>
    
    % check for largest diffusivity among cells/chemicals
    muMax = max(MFTmodel.mu);
    muMax = max(muMax,GlobalParams.DA);
    if isfield(GlobalParams,'DN')
        muMax = max(muMax,GlobalParams.DN);
    end
    
    
    % scale for cell density first set by using the initial condition
    if isfield(SimParams(end),'rho0')
        SimParams.rhoScale = max(sum(SimParams(end).rho0,1));
    else
        switch GlobalParams.IC
            case 'sigmoid'
                SimParams.rhoScale = GlobalParams.OD0;
            case 'step'
                SimParams.rhoScale = GlobalParams.OD0;
        end        
    end
    SimParams.tScale = GlobalParams.asp/(GlobalParams.kA*SimParams.rhoScale); % s, time scale of consumption
    SimParams.xScale = sqrt(SimParams.chiScale*SimParams.tScale);
    
    
    SimParams.dx = round(SimParams.xScale/x_discretize_factor,10); % %/40 /25
    SimParams.dt = round(SimParams.dx^2/muMax/t_discretize_factor,10); % % / 5
    
    SimParams.outDt = round(5*SimParams.tScale,3);% max(SimParams.outDt,SimParams.dt);
    SimParams.outDt = min(SimParams.outDt,3*60); % max 3 min outDt
    
    switch GlobalParams.IC
        case 'sigmoid'
            SimParams.segmentLength = 10*SimParams.xScale; % to start %%%%%%%%%%%%%%%%%%%%%
        case 'step'
            SimParams.segmentLength = 6*SimParams.xScale; % to start %%%%%%%%%%%%%%%%%%%%%
    end
    
    SimParams.time    = 0;
    SimParams.x       = 0:SimParams.dx:SimParams.segmentLength; % um, position along channel
    SimParams.L       = length(SimParams.x);
    
    
elseif nargin>=3 && ~isempty(varargin{1})
    % continuing from an existing simulation
    
    firstInEnv = false; % assume that this is continuing an ongoing simulation
    if nargin>=4 && ~isempty(varargin{2})
        % but could be initializing a new simulation in a new environment.
        % this checks.
        firstInEnv = varargin{2};
    end
    
    % need to recompute chiScale, Ncells, scales, using values from the wave
    % use values from the last cell element of SimResults and SimParams and start the
    % next cell
    
    % get last results
    SimResults = varargin{1};
    
    step = length(SimParams);
    
    % new SimParams entry
    % SimParams and SimResults should be the same length
    f = fields(SimParams);
    for i = 1:length(f)
        val = getfield(SimParams,{step},f{i});
        SimParams = setfield(SimParams,{step+1},f{i},val);
    end
    
    % compute MFTmodel
    MFTmodel = initializeMFTmodel(SimParams(step+1).F0, GlobalParams.environment, GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
    
    % 
    totLengthFactorAgar = 1/3.8;

    % need to know whether the previous block was the same environment or
    % not (same simulation or not)
    % if switching environments, set new environment
    if firstInEnv
        switch GlobalParams.run_case
            case 'switch_environments'
                % verify the case
                
                % compute MFTmodel in the previous environment
                switch GlobalParams.environment
                    case 'liquid'
                        MFTmodel_prev = initializeMFTmodel(SimParams(step).F0, 'agar', GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
                        if GlobalParams.k>1
                            GlobalParams.totLength = 1/totLengthFactorAgar * GlobalParams.totLength;
                        end

                    case 'agar'
                        MFTmodel_prev = initializeMFTmodel(SimParams(step).F0, 'liquid', GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
                        if GlobalParams.k>1
                            GlobalParams.totLength = totLengthFactorAgar * GlobalParams.totLength;
                        end
                end
        end
        
        
    else
        % alternatively, ~firstInEnv, even if run_case =
        % 'switch_environments', then we are just preparing to simulate the
        % next block of the same environment. still need MFTmodel_prev:
        MFTmodel_prev = initializeMFTmodel(SimParams(step).F0, GlobalParams.environment, GlobalParams.biophysicalChi, GlobalParams.alpha, GlobalParams.pore_size);
    end
    
    % find the back of the wave: where asp is tiny
    asp_prev = SimResults(end).asp(:,end);
    rho_prev = squeeze(SimResults(end).rho(:,:,end));
    
    % check whether there is a wave...
    [~,locs] = findpeaks(squeeze(sum(rho_prev,1)),'MinPeakHeight',0.1,'MinPeakProminence',0.01);
    
    if isempty(locs)
        disp('No wave found (yet).')
    end
    
    
    % find the back of the wave
    backInd = find(asp_prev<1e-3,1,'last');
    if isempty(backInd)
        backInd = 1;
    end
    
    % sum up cells per cross-sectional area in wave
    Ncells_i = squeeze(sum(rho_prev(:,backInd:end),2))*SimParams(step).dx;
    SimParams(step+1).Ncells = sum(Ncells_i);
    
    
    % compute scales
    SimParams(step+1).chiScale = sum(Ncells_i(:).*(MFTmodel.chi(:)-MFTmodel.mu(:)))/sum(Ncells_i(:)); % um^2/s % (alpha-1)*<mu>
    
    % check for largest mu among represented cells (more than 1/100th of
    % the traveling population)
    muMax = max(MFTmodel.mu); %max(MFTmodel.mu(Ncells_i(:)>0.01*max(Ncells_i(:))));
    muMax = max(muMax,GlobalParams.DA);
    if isfield(GlobalParams,'DN')
        muMax = max(muMax,GlobalParams.DN);
    end
    
    % expected speed
    % this corrects the (steady state...) speed for new environments, improving choices of dx and dt
    % minimal effect when environment doesn't change
    speed_factor = min(1,sqrt(sum(Ncells_i(:).*MFTmodel.chi(:))./sum(Ncells_i(:).*MFTmodel_prev.chi(:))));
    
    SimParams(step+1).c = speed_factor*SimParams(step+1).Ncells*GlobalParams.kA/GlobalParams.asp; % cells/um^2 * uM/OD/s / uM * (OD/(cell/mL)) * um^3/mL = um/s;
    
    % set rhoScale using the simulation profile. use peak density and Ncells to
    % choose xScale.
    rhotot_prev = sum(rho_prev,1);
    [pk,locs] = findpeaks(rhotot_prev(:,backInd:end),'MinPeakProminence',0.01);
    if ~isempty(locs)
        [~,mxi] = max(locs);
        SimParams(step+1).rhoScale = pk(mxi);
        ind = backInd-1+locs(mxi);
        xScale3 = 10*1/sqrt(abs((rhotot_prev(ind+1)+rhotot_prev(ind-1)-2*rhotot_prev(ind))/SimParams(step).dx^2));
        % use curvature around the peak to resolve it
    else
        SimParams(step+1).rhoScale = max(rhotot_prev(:,backInd:end));
        xScale3 = Inf;
    end

    SimParams(step+1).tScale = GlobalParams.asp/(GlobalParams.kA*SimParams(step+1).rhoScale);
    
    xScale1=speed_factor*SimParams(step+1).Ncells/SimParams(step+1).rhoScale;
    xScale2=sqrt(SimParams(step+1).chiScale*SimParams(step+1).tScale);
    
    SimParams(step+1).xScale = min(xScale1,xScale2);
    SimParams(step+1).xScale = min(SimParams(step+1).xScale,xScale3);
    
    % new grid sizes
    SimParams(step+1).dx = round(SimParams(step+1).xScale/x_discretize_factor,10); 
    SimParams(step+1).dt = round(SimParams(step+1).dx^2/muMax/t_discretize_factor,10);
    
    SimParams(step+1).outDt = round(5*SimParams(step+1).tScale);
    SimParams(step+1).outDt = min(SimParams(step+1).outDt,3*60); % max 3 min outDt
    
    SimParams(step+1).time = [];
    
    % indices of positions from the previous simulation block to keep in
    % the next block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % only start chopping off the back once the peak is sufficiently far
    % from the initial point
    if ~isempty(locs) && SimParams(step).x(ind)>20*SimParams(1).xScale
        if strcmp(GlobalParams.run_case,'single_environment')==1  || GlobalParams.DA==0
            keepInds = find(SimParams(step).x > (SimParams(step).x(backInd) - 10*SimParams(step+1).xScale)); % - SimParams(step).xScale
        elseif GlobalParams.DA~=0
            keepInds = find(SimParams(step).x > (SimParams(step).x(backInd) - 10*SimParams(step+1).xScale - SimParams(step).xScale)); % this seemed to work in a previous version of the code
        end
    else
        keepInds = 1:length(SimParams(step).x);
    end
    
    
    % define x that overlaps with previous x and extend it, with new dx
    x1 = SimParams(step).x(keepInds(1));
    xasp = SimParams(step).x(find(asp_prev<(1-1e-6)*GlobalParams.asp,1,'last'));
    x2 = xasp+5*SimParams(step+1).xScale; % extend space as little as necessary.
    SimParams(step+1).x       = x1:SimParams(step+1).dx:x2; % um, position along channel
    
    SimParams(step+1).L       = length(SimParams(step+1).x);
    SimParams(step+1).segmentLength = SimParams(step+1).x(end)-SimParams(step+1).x(1);
    
    % interpolate last simulation time point onto new x grid
    ind1 = find(SimParams(step+1).x>SimParams(step).x(keepInds(end)),1,'first'); % need this because interpolation will fail ahead of the wave. no observations there, so the interpolation function doesn't know what to do. fill in with initial condition.
    nPhen = size(rho_prev,1);
    SimParams(step+1).rho0 = nan(nPhen,SimParams(step+1).L);
    for i = 1:nPhen
        SimParams(step+1).rho0(i,:) = interp1(SimParams(step).x(keepInds),rho_prev(i,keepInds),SimParams(step+1).x,'pchip');
        SimParams(step+1).rho0(i,ind1:end) = 0;
    end
    SimParams(step+1).asp0 = interp1(SimParams(step).x(keepInds),asp_prev(keepInds),SimParams(step+1).x,'pchip');
    SimParams(step+1).asp0(ind1:end) = GlobalParams.asp;
    
    
    % only for case with diversity and with growth, add or remove phenotypes at the edges of the range
    if size(rho_prev,1)>1 && GlobalParams.growthRate>0
        F0 = SimParams(step).F0;
        PF0 = Ncells_i/nansum(Ncells_i);
        
        PF0_prevprev = [];
        
        % needed quantities
        if length(SimParams(step).time)>1 % there will be cases where this won't hold, and you'll have to check back a step...
            asp_prevprev = SimResults(end).asp(:,end-1);
            rho_prevprev = squeeze(SimResults(end).rho(:,:,end-1));
            backInd_prevprev = find(asp_prevprev<1e-3,1,'last');
            NF0_prevprev = sum(rho_prevprev(:,backInd_prevprev:end),2);
            PF0_prevprev = NF0_prevprev/sum(NF0_prevprev); % same F0
            
        elseif step>1
            
            rho_stepminus1 = squeeze(SimResults(step-1).rho(:,:,end));
            asp_stepminus1 = SimResults(step-1).asp(:,end);
            F0_stepminus1 = SimParams(step-1).F0;
            backInd_stepminus1 = find(asp_stepminus1<1e-3,1,'last');
            NF0_stepminus1 = sum(rho_stepminus1(:,backInd_stepminus1:end),2);
            PF0_stepminus1 = NF0_stepminus1/sum(NF0_stepminus1); % same dF0, so sum is ok
            
            % only need to check that P(F0) for F0's currently being
            % simulated are decreasing
            PF0_prevprev = interp1(F0_stepminus1,PF0_stepminus1,F0);
            PF0_prevprev(isnan(PF0_prevprev)) = 0;
            
        end
        
        
        if ~isempty(PF0_prevprev)
            % remove phenotypes
            phenRemoveInds = find(PF0/max(PF0)<=1e-5 & (PF0-PF0_prevprev)<0);
            if ~isempty(phenRemoveInds)
                if phenRemoveInds(1)<length(F0)/2
                    phenRemoveInds = [(1:(phenRemoveInds(1)))';phenRemoveInds];
                end
                if phenRemoveInds(end)>length(F0)/2
                    phenRemoveInds = [phenRemoveInds;((phenRemoveInds(end)):length(F0))'];
                end
            end
            phenRemoveInds = unique(phenRemoveInds);
                        
            SimParams(step+1).rho0(phenRemoveInds,:,:) = [];
            SimParams(step+1).F0(phenRemoveInds,:,:) = [];
            
            F0(phenRemoveInds) = [];
            PF0(phenRemoveInds) = [];
            PF0 = PF0/sum(PF0);
            
            % add phenotypes
            offspringStd = GlobalParams.stdF0*sqrt(1-GlobalParams.phi.^2);
           
            phenAddInds1 = PF0(1)/max(PF0)>9e-5 & (PF0(1)-PF0_prevprev(1))>0;
            if phenAddInds1
                % add gridpoints for new phenotypes
                
                F0 = SimParams(step+1).F0;
                dF0 = min(diff(F0));
                newpts = (F0(1)-1/2*offspringStd):dF0:(F0(1)-dF0);
                F0 = [newpts(:); F0];
                SimParams(step+1).F0 = F0;
                SimParams(step+1).rho0 = cat(1,zeros(length(newpts),size(SimParams(step+1).rho0,2)),SimParams(step+1).rho0);
            end
            
            phenAddInds2 = PF0(end)/max(PF0)>9e-5 & (PF0(end)-PF0_prevprev(end))>0;
            if phenAddInds2
                % add gridpoints for new phenotypes
                
                F0 = SimParams(step+1).F0;
                dF0 = min(diff(F0));
                newpts = (F0(end)+dF0):dF0:(F0(end)+1/2*offspringStd);
                F0 = [F0; newpts(:)];
                SimParams(step+1).F0 = F0;
                SimParams(step+1).rho0 = cat(1,SimParams(step+1).rho0,zeros(length(newpts),size(SimParams(step+1).rho0,2)));
            end
        end
        
    end
    
end

end

