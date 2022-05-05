function varargout = mainFunction(GlobalParams,SimParams,varargin)
% GlobalParams is a structure containing parameter values for the
% simulation
% SimParams can be an empty variable
% varargin: optional input file name (like 'results.mat') to check whether the simulation results file exists
% if GlobalParams and SimParams are empty variables, this will run a
% simulation with default parameter values

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end

[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

GlobalParams.k = 1;

if length(SimParams)<=1 && ~isfield(SimParams,'time')
    [GlobalParams, SimParams] = initializeSimulationStructures(GlobalParams, SimParams);
end

if nargout>1
    varargout{2} = GlobalParams;
end

% check whether simulation has been run or not
if ~isempty(varargin)
    filename = varargin{1};
    exists = exist([GlobalParams.dataDir, filename],'file')==2;
    if exists
        varargout{1} = 1;
    else
        %   simulation has not been run
        varargout{1} = 0;
    end
    return
end

switch GlobalParams.run_case
    case 'single_environment' 
        
        % check for existing simulation results
        if exist(GlobalParams.simName, 'file')==0 || GlobalParams.reSimulate==1
            % simplest case first: no simulation file exists yet, or
            % resimulate is set to 1
            
            % make save directory if it doesn't exist
            if ~exist(GlobalParams.dataDir, 'dir')
                mkdir(GlobalParams.dataDir)
            end
            
            % run simulation
            tic
            disp(['Running simulation ', GlobalParams.simName, '...'])
            [SimResults, SimParams] = simulateWave(GlobalParams, SimParams);
            toc
            
            % save
            disp(['Saving data file ', GlobalParams.simName, '...'])
            
            save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v6')
            matObj = matfile(GlobalParams.simName);
            vars = who(matObj);
            
            % check whether file save. if not, use v7.3
            if all(cellfun(@isempty, regexp(vars, 'SimResults')))
                save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v7.3')
            end
            
        elseif exist(GlobalParams.simName, 'file')==2 && GlobalParams.reSimulate==0
            % if a simulation file with the same parameters exists, and
            % option is set to NOT resimulate...
            
            disp(['Loading simulation ', GlobalParams.simName, '...'])
            
            % get input total simulation time
            totTime = GlobalParams.totTime;
            
            % load previous simulation data
            s = load(GlobalParams.simName);
            GlobalParams_last = s.GlobalParams;
            SimParams_last = s.SimParams;
            
            
            if abs(SimParams_last(end).time(end)-totTime)<=SimParams_last(end).outDt
                disp('Simulation already run to completion.')
                return
            end
            
            % continue simulation where it left off
            GlobalParams_last.totTime = totTime;
            GlobalParams = GlobalParams_last;
            SimParams = SimParams_last;
            SimResults = s.SimResults;
            
            [GlobalParams, SimParams, flag] = prepExtendSimulationSS(GlobalParams, SimParams, SimResults);
            
            if flag==1
                tic
                disp(['Extending simulation ', GlobalParams.simName, '...'])
                % may double count the last time point of the previous
                    % simulation here?
                [SimResults, SimParams] = simulateWave(GlobalParams, SimParams, SimResults);
                
                toc
            end
            
            % save
            disp(['Saving data file ', GlobalParams.simName, '...'])
            
            save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v6')
            matObj = matfile(GlobalParams.simName);
            vars = who(matObj);
            
            % check whether file save. if not, use v7.3
            if all(cellfun(@isempty, regexp(vars, 'SimResults')))
                disp('Saving variables in .mat -v7.3 format.')
                save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v7.3')
            end
            
            
        end
        
        
    case 'switch_environments'
        % switching environments
        nEnv = GlobalParams.nEnv;
        
        % check for data files
        files = dir(GlobalParams.dataDir);
        files = files(~cellfun(@isempty, regexp({files(:).name}, '.mat')));
        files = files(~cellfun(@isempty, regexp({files(:).name}, 'results_env_')));
        
        % also need to check that the last existing file ran to completion....
        
        if ~isempty(files) && GlobalParams.reSimulate==0
            % if files exist, load where it left off last
            disp('Simulation files exist.')
            K = nan(size(files));
            for i = 1:length(files)
                [~,ind1] = regexp(files(i).name, 'results_env_');
                ind2 = regexp(files(i).name, '.mat');
                K(i) = str2double(files(i).name((ind1+1):(ind2-1)));
            end
            
            [lastK,mxi] = max(K);
            disp(['Checking simulation ' GlobalParams.dataDir, filesep, files(mxi).name])
            s=load([GlobalParams.dataDir, filesep, files(mxi).name],'SimParams','SimResults');
            
            if s.SimParams(end).time(end)<GlobalParams.totTime
                disp('Last simulation has not yet been run to completion.')
                
                s=load([GlobalParams.dataDir, filesep, files(mxi).name]);
                s.GlobalParams.totTime = GlobalParams.totTime;
                s.GlobalParams.nEnv = GlobalParams.nEnv;
                GlobalParams = s.GlobalParams;
                SimParams = s.SimParams;
                SimResults = s.SimResults;
                
                [GlobalParams, SimParams, flag] = prepExtendSimulationSS(GlobalParams, SimParams, SimResults);
                
                if flag==1
                    tic
                    disp(['Extending simulation ', GlobalParams.simName, '...'])
                    % may double count the last time point of the previous
                    % simulation here?
                    [SimResults, SimParams] = simulateWave(GlobalParams, SimParams, SimResults);
                    
                    toc
                end
                
                % save
                disp(['Saving data file ', GlobalParams.simName, '...'])
                
                save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v6')
                matObj = matfile(GlobalParams.simName);
                vars = who(matObj);
                
                % check whether file save. if not, use v7.3
                if all(cellfun(@isempty, regexp(vars, 'SimResults')))
                    disp('Saving variables in .mat -v7.3 format.')
                    save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v7.3')
                end
                
                
                
            end

            
            
            if lastK == nEnv
                disp('All simulations have been run!')
                return
            else
                disp('Continuing from where the last simulation left off...')
                k0 = lastK+1;
            end
        else
            % if no environments have been run yet, run the first one, k = 1
            
            % make save directory if it doesn't exist
            if ~exist(GlobalParams.dataDir, 'dir')
                mkdir(GlobalParams.dataDir)
            end
            
            % run simulation
            tic
            disp(['Running simulation ', GlobalParams.simName, '...'])
            [SimResults, SimParams] = simulateWave(GlobalParams, SimParams);
            toc
            
            % save
            disp(['Saving data file ', GlobalParams.simName, '...'])
            save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v6')
            matObj = matfile(GlobalParams.simName);
            vars = who(matObj);
            
            % check whether file save. if not, use v7.3
            if all(cellfun(@isempty, regexp(vars, 'SimResults')))
                disp('Saving variables in .mat -v7.3 format.')
                save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v7.3')
            end
            
            k0 = 2;
        end
        
        
        for k = k0:nEnv
        
            % load the previous data file
            disp(['Loading simulation ' GlobalParams.dataDir, filesep, sprintf('results_env_%03d',k-1), '.mat'])
            load([GlobalParams.dataDir, filesep, sprintf('results_env_%03d',k-1), '.mat'])

            GlobalParams.nEnv = nEnv;
            
            % prepare next simulation
            GlobalParams.k = k;
            [GlobalParams, SimParams, flag] = prepExtendSimulationEnvironments(GlobalParams, SimParams, SimResults);
            % now GlobalParams corresponds to simulation k
            
            % run simulation with options
            if flag==1 && (exist(GlobalParams.simName,'file')==0 || GlobalParams.reSimulate==1)
                % run simulation
                tic
                disp(['Running simulation ', GlobalParams.simName, '...'])
                [SimResults, SimParams] = simulateWave(GlobalParams, SimParams);
                toc
                
                %% save
                disp(['Saving data file ', GlobalParams.simName, '...'])
                save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v6')
                matObj = matfile(GlobalParams.simName);
                vars = who(matObj);
                
                % check whether file save. if not, use v7.3
                if all(cellfun(@isempty, regexp(vars, 'SimResults')))
                    disp('Saving variables in .mat -v7.3 format.')
                    save(GlobalParams.simName, 'GlobalParams', 'SimParams', 'SimResults', '-v7.3')
                end
                
%                 only keep last two iterations of EACH environment
                if k>=5
                    delete([GlobalParams.dataDir, filesep, sprintf('results_env_%03d',k-4), '.mat']);
                end
                  
            else
                disp('Error. Simulation will not be run.')
                return
                
            end
            
        end
        
end


end