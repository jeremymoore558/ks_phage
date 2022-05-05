function plotSwitchingEnvQuantities(stdTB,phis,totTimes,biophysicalChi,DA,xmesh)

% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineMarkerSize',5);
set(0,'DefaultLineColor',[0 0 0]);
set(0,'DefaultAxesBox','off');
set(0,'DefaultAxesColor','none');
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultTextInterpreter','tex');
set(0,'DefaultAxesTickLabelInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex');
whitebg('white');
close

%% global
close all

doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.environment = 'liquid';
GlobalParams.peakTB = 0.140;
GlobalParams.run_case = 'switch_environments';
GlobalParams.nEnv = 30;
SimParams = [];

GlobalParams.stdTB = stdTB;
GlobalParams.biophysicalChi = biophysicalChi;
GlobalParams.DA = DA;

if ~isempty(xmesh)
    GlobalParams.x_discretize_factor = xmesh;
end

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(GlobalParams.stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks',filesep]);
ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns',filesep,'figures'];

dense_factor = 10;

r = GlobalParams.growthRate;

%%
t_gt = cell(length(phis),length(totTimes));
c_gt = cell(length(phis),length(totTimes));
env_gt = cell(length(phis),length(totTimes)); % 1 for liquid, 2 for agar
cavg_gt = nan(length(phis),length(totTimes));
colors = jet(GlobalParams.nEnv);
for tt = 1:length(totTimes)
    disp(['Env duration r*T = ' num2str(totTimes(tt)*r)])
    for g = 1:length(phis)
        tic
        
        Dx = nan(GlobalParams.nEnv,1);
        Dt = nan(GlobalParams.nEnv,1);
        
        SimParams_gt = [];
        GlobalParams_gt = GlobalParams;
        GlobalParams_gt.phi = phis(g);
        GlobalParams_gt.totTime = totTimes(tt);
        [GlobalParams_gt,~] = initializeSimulationStructures(GlobalParams_gt,SimParams_gt);
        
        files = dir(GlobalParams_gt.dataDir);
        files = files(cellfun(@(x) startsWith(x,'results_env'), {files.name}));
        
        simK = nan(length(files),1);
        for i = 1:length(files)
            simK(i) = str2double(files(i).name(end-6:end-4));
        end
        
        maxK = max(simK)-1;
        if mod(maxK,2)~=0
            maxK = maxK-1;
        end
        
        if tt == 1 || tt==length(totTimes)
            figure;hold on;htb=gca;ftb=gcf;
            title(['\tau = ' num2str(phis(g)) ', rT = ' num2str(r*totTimes(tt))]) % could add back
        end
        for k = (maxK-1):maxK
            % put this inside the loop so that names come out right
            SimParams_gt = [];
            GlobalParams_gt = GlobalParams;
            GlobalParams_gt.phi = phis(g);
            GlobalParams_gt.totTime = totTimes(tt);
            [GlobalParams_gt,~] = initializeSimulationStructures(GlobalParams_gt,SimParams_gt);
            
            [~,sind2] = regexp(GlobalParams_gt.simName,'results_env_');
            GlobalParams_gt.simName = [GlobalParams_gt.simName(1:sind2), sprintf('%03d', k), '.mat'];
            
            s=load(GlobalParams_gt.simName);
            SimParams_gt = s.SimParams;
            GlobalParams_gt = s.GlobalParams;
            SimResults_gt = s.SimResults;
            
            x = SimParams_gt(end).x;
            rho_t = SimResults_gt(end).rho;
            asp_t = SimResults_gt(end).asp;
            x_dense = linspace(min(x),max(x),length(x)*dense_factor);
            dx_dense = SimParams_gt(end).dx/dense_factor;
            asp_i = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');
            rho_i = zeros(size(rho_t,1),length(x_dense));
            for kk = 1:size(rho_t,1)
                rho_i(kk,:) = interp1(x,squeeze(rho_t(kk,:,end)),x_dense,'pchip');
            end
            [pks,locs] = findpeaks(squeeze(sum(rho_i,1)));
            
            if ~isempty(locs)
                SimResults_gt(end).peakLocs(end) = x_dense(locs(end));
            end
            
            Dx(k) = SimResults_gt(end).peakLocs(end)-SimResults_gt(1).peakLocs(1); %%%%
            if k==1
                Dx(k) = SimResults_gt(end).peakLocs(end);
            end
            Dt(k) = SimParams_gt(end).time(end);
            
            
            if tt == 1 || tt==length(totTimes)
                backInd = find(asp_i<sqrt(GlobalParams_gt.KmA*GlobalParams_gt.KiA),1,'last');
                
                if isempty(backInd)
                    backInd = 1;
                end
                
                N = dx_dense*sum(rho_i(:,backInd:end),2);
                PF0 = N/sum(N,1);
                [TB,PTB] = convertPF0ToPTB(SimParams_gt(end).F0,PF0);
                plot(htb,TB,PTB)
                xlabel(htb,'Tumble bias, TB')
                ylabel(htb,'P(TB)')
                if k==maxK
                    hleg = legend('Liquid','Porous','Location','best');
                    hleg.Box = 'off';
                    print(ftb,[save_dir, filesep, 'PTBend_tau=' num2str(phis(g)) '_rT=' num2str(r*totTimes(tt)) '_', case_name, '.png'],fmt,rez)
                    savefig(ftb,[save_dir, filesep, 'PTBend_tau=' num2str(phis(g)) '_rT=' num2str(r*totTimes(tt)) '_', case_name, '.fig'])
                end
                
            end
        end
        
        cavg_gt(g,tt) = (Dx(maxK-1)+Dx(maxK)) / (Dt(maxK-1)+Dt(maxK));
        toc
        drawnow
    end
    close all
end

if length(phis)>1
    
    save([save_dir, filesep, 'outputs_plotSwitchingEnvQuantities_' case_name '.mat'],'Taus','T','r','cavg_gt');
    
else
    
    save([save_dir, filesep, 'outputs_plotSwitchingEnvQuantities_' case_name '.mat'],'r','cavg_gt');
    
end


end