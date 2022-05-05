function plotSteadyStateQuantities(stdTB,phis,biophysicalChi,DA,varargin)

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

GlobalParams = [];
doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.stdTB = stdTB;
GlobalParams.reSimulate = 0;
GlobalParams.biophysicalChi = biophysicalChi;
GlobalParams.DA = DA;
SimParams = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];

if nargin>=5
    if ~isempty(varargin{1})
        GlobalParams.growthRate = varargin{1};
        case_name = [case_name '_r=' num2str(round(3600*GlobalParams.growthRate,3)) 'perHr'];
    end
end

if nargin>=6
    if ~isempty(varargin{2})
        GlobalParams.asp = varargin{2};
        case_name = [case_name '_Asp=' num2str(GlobalParams.asp) 'uM'];
    end
end

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


TB_global = linspace(0,0.55,1000);
F0_global = linspace(0,5,1000);


%% liquid
GlobalParams.environment = 'liquid';

css_liq_g = nan(size(phis));
chiBack_liq_g = nan(size(phis));
muBack_liq_g = nan(size(phis));
chi_liq_g = nan(size(phis));
mu_liq_g = nan(size(phis));
t_liq_g = cell(size(phis));
c_liq_g = cell(size(phis));
PF0t_liq_g = cell(size(phis));

TB_liq = cell(size(phis));
PTBss_liq = cell(size(phis));
for g = length(phis):-1:1
    tic
    SimParams_g = [];
    GlobalParams_g = GlobalParams;
    GlobalParams_g.phi = phis(g);
    
    [GlobalParams_g,~] = initializeSimulationStructures(GlobalParams_g,SimParams_g);
    s=load(GlobalParams_g.simName);
    SimParams_g = s.SimParams;
    GlobalParams_g = s.GlobalParams;
    SimResults_g = s.SimResults;
    
    F0 = SimParams_g(end).F0;
    TB = 1./(1+exp(F0));
    x = SimParams_g(end).x;
    rho_t = SimResults_g(end).rho;
    asp_t = squeeze(SimResults_g(end).asp);
    
    x_dense = linspace(min(x),max(x),length(x)*dense_factor);
    dx_dense = SimParams_g(end).dx/dense_factor;
    asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');
    
    rho_end = zeros(size(rho_t,1),length(x_dense));
    for k = 1:size(rho_t,1)
        rho_end(k,:) = interp1(x,squeeze(rho_t(k,:,end)),x_dense,'pchip');
    end
    
    backInd = find(asp_end<sqrt(GlobalParams_g.KmA*GlobalParams_g.KiA),1,'last');
    N = dx_dense*sum(rho_end(:,backInd:end),2);
    PF0_ss = N/trapz(F0,N);
    dTBdF0 = TB.*(1-TB);
    PTB_ss = PF0_ss./dTBdF0;
    PTB_ss = PTB_ss/(-trapz(TB,PTB_ss));
    
    TB_liq{g} = TB;
    PTBss_liq{g} = PTB_ss;
    
    css_liq_g(g) = sum(N,1)*GlobalParams_g.kA/GlobalParams_g.asp;
    
    Nback = dx_dense*rho_end(:,backInd);
    PF0back = Nback/trapz(F0,Nback);
    MFT_g = initializeMFTmodel(F0,GlobalParams_g.environment,GlobalParams_g.biophysicalChi);
    chiBack_liq_g(g) = trapz(F0,MFT_g.chi.*PF0back);
    muBack_liq_g(g) = trapz(F0,MFT_g.mu.*PF0back);
    
    chi_liq_g(g) = trapz(F0,MFT_g.chi.*PF0_ss);
    mu_liq_g(g) = trapz(F0,MFT_g.mu.*PF0_ss);
    
    % get phenotype distributions on a common grid of F0
    
    for i = 1:length(SimParams_g)
        F0i = SimParams_g(i).F0;
        TBi = 1./(1+exp(F0i));
        x = SimParams_g(i).x;
        rho_t = SimResults_g(i).rho;
        asp_t = squeeze(SimResults_g(i).asp);
        
        t_liq_g{g} = [t_liq_g{g};SimParams_g(i).time(:)];
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        dx_dense = SimParams_g(i).dx/dense_factor;
        
        for j = 1:size(rho_t,3)
            asp_i = interp1(x(:),asp_t(:,j),x_dense(:),'pchip');
            rho_i = zeros(size(rho_t,1),length(x_dense));
            
            
            for k = 1:size(rho_t,1)
                rho_i(k,:) = interp1(x,squeeze(rho_t(k,:,j)),x_dense,'pchip');
            end
            
            backInd = find(asp_i<sqrt(GlobalParams_g.KmA*GlobalParams_g.KiA),1,'last');
            if isempty(backInd)
                backInd = 1;
            end
            N = dx_dense*sum(rho_i(:,backInd:end),2);
            
            PF0j = N/trapz(F0i,N);
            PF0t_j = interp1(F0i,PF0j,F0_global);
            PF0t_liq_g{g} = [PF0t_liq_g{g},PF0t_j(:)];
            
            c_i = sum(N,1)*GlobalParams_g.kA/GlobalParams_g.asp;
            c_liq_g{g} = [c_liq_g{g}(:); c_i];
            
        end
    end
    toc
end


%% agar
GlobalParams.environment = 'agar';

css_agar_g = nan(size(phis));
chiBack_agar_g = nan(size(phis));
muBack_agar_g = nan(size(phis));
chi_agar_g = nan(size(phis));
mu_agar_g = nan(size(phis));
t_agar_g = cell(size(phis));
c_agar_g = cell(size(phis));
PF0t_agar_g = cell(size(phis));

TB_ag = cell(size(phis));
PTBss_ag = cell(size(phis));
for g = 1:length(phis)
    tic
    SimParams_g = [];
    GlobalParams_g = GlobalParams;
    GlobalParams_g.phi = phis(g);
    
    [GlobalParams_g,~] = initializeSimulationStructures(GlobalParams_g,SimParams_g);
    s=load(GlobalParams_g.simName);
    SimParams_g = s.SimParams;
    GlobalParams_g = s.GlobalParams;
    SimResults_g = s.SimResults;
    
    F0 = SimParams_g(end).F0;
    TB = 1./(1+exp(F0));
    x = SimParams_g(end).x;
    rho_t = SimResults_g(end).rho;
    asp_t = squeeze(SimResults_g(end).asp);
    
    x_dense = linspace(min(x),max(x),length(x)*dense_factor);
    dx_dense = SimParams_g(end).dx/dense_factor;
    asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');
    
    rho_end = zeros(size(rho_t,1),length(x_dense));
    for k = 1:size(rho_t,1)
        rho_end(k,:) = interp1(x,squeeze(rho_t(k,:,end)),x_dense,'pchip');
    end
    
    backInd = find(asp_end<sqrt(GlobalParams_g.KmA*GlobalParams_g.KiA),1,'last');
    N = dx_dense*sum(rho_end(:,backInd:end),2);
    PF0_ss = N/trapz(F0,N);
    dTBdF0 = TB.*(1-TB);
    PTB_ss = PF0_ss./dTBdF0;
    PTB_ss = PTB_ss/(-trapz(TB,PTB_ss));
    
    TB_ag{g} = TB;
    PTBss_ag{g} = PTB_ss;
    
    css_agar_g(g) = sum(N,1)*GlobalParams_g.kA/GlobalParams_g.asp;
    
    Nback = dx_dense*rho_end(:,backInd);
    PF0back = Nback/trapz(F0,Nback);
    MFT_g = initializeMFTmodel(F0,GlobalParams_g.environment,GlobalParams_g.biophysicalChi);
    chiBack_agar_g(g) = trapz(F0,MFT_g.chi.*PF0back);
    muBack_agar_g(g) = trapz(F0,MFT_g.mu.*PF0back);
    
    chi_agar_g(g) = trapz(F0,MFT_g.chi.*PF0_ss);
    mu_agar_g(g) = trapz(F0,MFT_g.mu.*PF0_ss);
    
    
    % get phenotype distributions on a common grid of F0
    for i = 1:length(SimParams_g)
        F0i = SimParams_g(i).F0;
        TBi = 1./(1+exp(F0i));
        x = SimParams_g(i).x;
        rho_t = SimResults_g(i).rho;
        asp_t = squeeze(SimResults_g(i).asp);
        
        t_agar_g{g} = [t_agar_g{g};SimParams_g(i).time(:)];
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        dx_dense = SimParams_g(i).dx/dense_factor;
        
        for j = 1:size(rho_t,3)
            asp_i = interp1(x(:),asp_t(:,j),x_dense(:),'pchip');
            rho_i = zeros(size(rho_t,1),length(x_dense));
            
            for k = 1:size(rho_t,1)
                rho_i(k,:) = interp1(x,squeeze(rho_t(k,:,j)),x_dense,'pchip');
            end
            
            %             backInd = find(asp_i<1e-3,1,'last');
            backInd = find(asp_i<sqrt(GlobalParams_g.KmA*GlobalParams_g.KiA),1,'last');
            N = dx_dense*sum(rho_i(:,backInd:end),2);
            PF0j = N/trapz(F0i,N);
            PF0t_j = interp1(F0i,PF0j,F0_global);
            PF0t_agar_g{g} = [PF0t_agar_g{g},PF0t_j(:)];
            
            c_i = sum(N,1)*GlobalParams_g.kA/GlobalParams_g.asp;
            c_agar_g{g} = [c_agar_g{g}(:); c_i];
        end
        
    end
    toc
    
end


%% simulations without diversity



% generalist in liquid
GlobalParams_gen0 = GlobalParams;
GlobalParams_gen0.peakTB = generalistTB;
GlobalParams_gen0.stdTB = 0;
GlobalParams_gen0.phi = 0;
GlobalParams_gen0.environment = 'liquid';
SimParams_gen = [];

[GlobalParams_gen,~] = initializeSimulationStructures(GlobalParams_gen0,SimParams_gen);
s=load(GlobalParams_gen.simName);
SimParams_gen = s.SimParams;
GlobalParams_gen = s.GlobalParams;
SimResults_gen = s.SimResults;

F0 = SimParams_gen(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_gen(end).x;
rho_t = SimResults_gen(end).rho;
asp_t = squeeze(SimResults_gen(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_gen(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_gen.KmA*GlobalParams_gen.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_liq_gen = sum(N,1)*GlobalParams_gen.kA/GlobalParams_gen.asp;

% generalist in agar
GlobalParams_gen0.environment = 'agar';
SimParams_gen = [];

[GlobalParams_gen,~] = initializeSimulationStructures(GlobalParams_gen0,SimParams_gen);
s=load(GlobalParams_gen.simName);
SimParams_gen = s.SimParams;
GlobalParams_gen = s.GlobalParams;
SimResults_gen = s.SimResults;

F0 = SimParams_gen(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_gen(end).x;
rho_t = SimResults_gen(end).rho;
asp_t = squeeze(SimResults_gen(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_gen(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_gen.KmA*GlobalParams_gen.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_ag_gen = sum(N,1)*GlobalParams_gen.kA/GlobalParams_gen.asp;


% liquid specialist in liquid
GlobalParams_liq0 = GlobalParams_gen0;
GlobalParams_liq0.peakTB = liquid_specialist_TB;

GlobalParams_liq0.environment = 'liquid';
SimParams_liq = [];

[GlobalParams_liq,~] = initializeSimulationStructures(GlobalParams_liq0,SimParams_liq);
s=load(GlobalParams_liq.simName);
SimParams_liq = s.SimParams;
GlobalParams_liq = s.GlobalParams;
SimResults_liq = s.SimResults;

F0 = SimParams_liq(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_liq(end).x;
rho_t = SimResults_liq(end).rho;
asp_t = squeeze(SimResults_liq(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_liq(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_liq.KmA*GlobalParams_liq.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_liq_liqspec = sum(N,1)*GlobalParams_liq.kA/GlobalParams_liq.asp;

% liquid specialist in agar
GlobalParams_liq0.environment = 'agar';
SimParams_liq = [];

[GlobalParams_liq,~] = initializeSimulationStructures(GlobalParams_liq0,SimParams_liq);
s=load(GlobalParams_liq.simName);
SimParams_liq = s.SimParams;
GlobalParams_liq = s.GlobalParams;
SimResults_liq = s.SimResults;

F0 = SimParams_liq(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_liq(end).x;
rho_t = SimResults_liq(end).rho;
asp_t = squeeze(SimResults_liq(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_liq(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_liq.KmA*GlobalParams_liq.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_ag_liqspec = sum(N,1)*GlobalParams_liq.kA/GlobalParams_liq.asp;


% agar specialist in liquid
GlobalParams_ag0 = GlobalParams_liq0;
GlobalParams_ag0.peakTB = agar_specialist_TB;

GlobalParams_ag0.environment = 'liquid';
SimParams_ag = [];

[GlobalParams_ag,~] = initializeSimulationStructures(GlobalParams_ag0,SimParams_ag);
s=load(GlobalParams_ag.simName);
SimParams_ag = s.SimParams;
GlobalParams_ag = s.GlobalParams;
SimResults_ag = s.SimResults;

F0 = SimParams_ag(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_ag(end).x;
rho_t = SimResults_ag(end).rho;
asp_t = squeeze(SimResults_ag(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_ag(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_ag.KmA*GlobalParams_ag.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_liq_agspec = sum(N,1)*GlobalParams_ag.kA/GlobalParams_ag.asp;

% agar specialist in agar
GlobalParams_ag0.environment = 'agar';
SimParams_ag = [];

[GlobalParams_ag,~] = initializeSimulationStructures(GlobalParams_ag0,SimParams_ag);
s=load(GlobalParams_ag.simName);
SimParams_ag = s.SimParams;
GlobalParams_ag = s.GlobalParams;
SimResults_ag = s.SimResults;

F0 = SimParams_ag(end).F0;
TB = 1./(1+exp(F0));
x = SimParams_ag(end).x;
rho_t = SimResults_ag(end).rho;
asp_t = squeeze(SimResults_ag(end).asp);

x_dense = linspace(min(x),max(x),length(x)*dense_factor);
dx_dense = SimParams_ag(end).dx/dense_factor;
asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');

rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');

backInd = find(asp_end<sqrt(GlobalParams_ag.KmA*GlobalParams_ag.KiA),1,'last');
N = dx_dense*sum(rho_end(:,backInd:end),2);

css_ag_agspec = sum(N,1)*GlobalParams_ag.kA/GlobalParams_ag.asp;


% varying TB in liquid
TBs = logspace(-2,log10(0.4),50);
css_liq_TB = zeros(size(TBs));

GlobalParams_TB0 = GlobalParams_liq0;
GlobalParams_TB0.environment = 'liquid';

for tb = 1:length(TBs)
    
    GlobalParams_TB0.peakTB = TBs(tb);
    SimParams_TB = [];
    
    [GlobalParams_TB,~] = initializeSimulationStructures(GlobalParams_TB0,SimParams_TB);
    s=load(GlobalParams_TB.simName);
    SimParams_TB = s.SimParams;
    GlobalParams_TB = s.GlobalParams;
    SimResults_TB = s.SimResults;
    
    F0 = SimParams_TB(end).F0;
    TB = 1./(1+exp(F0));
    x = SimParams_TB(end).x;
    rho_t = SimResults_TB(end).rho;
    asp_t = squeeze(SimResults_TB(end).asp);
    
    x_dense = linspace(min(x),max(x),length(x)*dense_factor);
    dx_dense = SimParams_TB(end).dx/dense_factor;
    
    asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');
    rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');
    
    backInd = find(asp_end<sqrt(GlobalParams_TB.KmA*GlobalParams_TB.KiA),1,'last');
    N = dx_dense*sum(rho_end(:,backInd:end),2);
    
    css_liq_TB(tb) = sum(N,1)*GlobalParams_TB.kA/GlobalParams_TB.asp;
end

% in agar
GlobalParams_TB0.environment = 'agar';
css_ag_TB = zeros(size(TBs));

for tb = 1:length(TBs)
    
    if any(tb==1:4)
        css_ag_TB(tb) = 0;
    else
        GlobalParams_TB0.peakTB = TBs(tb);
        SimParams_TB = [];
        
        [GlobalParams_TB,~] = initializeSimulationStructures(GlobalParams_TB0,SimParams_TB);
        s=load(GlobalParams_TB.simName);
        SimParams_TB = s.SimParams;
        GlobalParams_TB = s.GlobalParams;
        SimResults_TB = s.SimResults;
        
        F0 = SimParams_TB(end).F0;
        TB = 1./(1+exp(F0));
        x = SimParams_TB(end).x;
        rho_t = SimResults_TB(end).rho;
        asp_t = squeeze(SimResults_TB(end).asp);
        
        x_dense = linspace(min(x),max(x),length(x)*dense_factor);
        dx_dense = SimParams_TB(end).dx/dense_factor;
        
        asp_end = interp1(x(:),asp_t(:,end),x_dense(:),'pchip');
        rho_end = interp1(x,squeeze(rho_t(1,:,end)),x_dense,'pchip');
        
        backInd = find(asp_end<sqrt(GlobalParams_TB.KmA*GlobalParams_TB.KiA),1,'last');
        N = dx_dense*sum(rho_end(:,backInd:end),2);
        
        css_ag_TB(tb) = sum(N,1)*GlobalParams_TB.kA/GlobalParams_TB.asp;
    end
end



%%
if exist([save_dir, filesep, 'outputs_plotSteadyStateQuantities_' case_name '.mat'],'file')==2
    save([save_dir, filesep, 'outputs_plotSteadyStateQuantities_' case_name '.mat'],'css_liq_gen','css_ag_gen','css_liq_liqspec','css_ag_liqspec','css_liq_agspec','css_ag_agspec','TBs','css_liq_TB','css_ag_TB','-append')
else
    save([save_dir, filesep, 'outputs_plotSteadyStateQuantities_' case_name '.mat'],'phis','css_liq_g','css_agar_g','chiBack_liq_g','chiBack_agar_g','chi_liq_g','chi_agar_g','muBack_liq_g','muBack_agar_g','mu_liq_g','mu_agar_g','c_agar_g','c_liq_g','t_agar_g','t_liq_g','TB_liq','PTBss_liq','TB_ag','PTBss_ag','css_liq_gen','css_ag_gen','css_liq_liqspec','css_ag_liqspec','css_liq_agspec','css_ag_agspec','css_liq_TB','css_ag_TB','TBs')
end
