%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultAxesFontName','Helvetica');
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
% close all

stdTB = 0.065;
biophysicalChi = 0;
DA = 800;

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

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

addpath([ks_dir, filesep, 'functions'])

save_dir = [ks_dir,filesep,'analysisFuns'];

dense_factor = 10;

load([save_dir, filesep, 'figures', filesep, 'outputs_plotSteadyStateQuantities_' case_name '.mat'])

plotInds = [1,6:2:length(taus)];
phis = exp(-1./(2*taus));

TBupper = 0.5;

fsize = [1.7,1.29];

fsize_slide = [4,3];


%% figure

% initial P(TB)
[GlobalParams,~] = initializeSimulationStructures(GlobalParams,[]);
TB = 1./(1+exp(GlobalParams.F0));
PF0_ss = GlobalParams.PF0;
dTBdF0 = TB.*(1-TB);
PTB_ss = PF0_ss./dTBdF0;
PTB_ss = PTB_ss/(-trapz(TB,PTB_ss));

ssTBfig_liq = figure; ssTBax_liq = gca;
hold on
xlim([0,TBupper])


green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];

liq_color_factors = logspace(log10(1/2),log10(1/max(colors(1,:))),length(taus));
ag_color_factors = logspace(log10(1/2),log10(1/max(colors(2,:))),length(taus));
liq_colors = repmat(colors(1,:),length(taus),1).*repmat(liq_color_factors(:),1,3);
ag_colors = repmat(colors(2,:),length(taus),1).*repmat(ag_color_factors(:),1,3);

TB_global = linspace(0,0.55,1000);
F0_global = linspace(0,5,1000);


%% liquid
GlobalParams.environment = 'liquid';

for ph = length(plotInds):-2:1
%     tic

    plot(ssTBax_liq,TB_liq{plotInds(ph)},PTBss_liq{plotInds(ph)},'Color',liq_colors(plotInds(ph),:)); 
    
    plot(ssTBax_liq_slide,TB_liq{plotInds(ph)},PTBss_liq{plotInds(ph)},'Color',liq_colors(plotInds(ph),:),'LineWidth',2); 
    
    plot(ssTBax_both_slide,TB_liq{plotInds(ph)},PTBss_liq{plotInds(ph)},'Color',liq_colors(plotInds(ph),:),'LineWidth',2); 
    
%     toc
end

%%
% original population
plot(ssTBax_liq,TB,PTB_ss,'k')


ssTBfig_ag = figure; ssTBax_ag = gca;
plot(ssTBax_ag,TB,PTB_ss,'k')
hold on
xlim([0,TBupper])



%% agar
GlobalParams.environment = 'agar';

for ph = 1:2:length(plotInds)
%     tic
    
    plot(ssTBax_ag,TB_ag{plotInds(ph)},PTBss_ag{plotInds(ph)},'Color',ag_colors(plotInds(ph),:)); 
    
    plot(ssTBax_ag_slide,TB_ag{plotInds(ph)},PTBss_ag{plotInds(ph)},'Color',ag_colors(plotInds(ph),:),'LineWidth',2); 
    plot(ssTBax_both_slide,TB_ag{plotInds(ph)},PTBss_ag{plotInds(ph)},'Color',ag_colors(plotInds(ph),:),'LineWidth',2); 
%     toc
end

%%
ssTBfig_liq.Units = 'inches';
ssTBfig_liq.Position(3:4) = fsize; 
ssTBax_liq.FontSize = 10;
ssTBax_liq.XTick = 0:0.2:TBupper;

ssTBax_liq.Box = 'off';
ssTBax_liq.YLim = [0,25];
ssTBax_liq.YTick = [];


ssTBfig_ag.Units = 'inches';
ssTBfig_ag.Position(3:4) = fsize; 
ssTBax_ag.FontSize = 10;
ssTBax_ag.XTick = 0:0.2:TBupper;
ssTBax_ag.Box = 'off';
ssTBax_ag.YTick = [];


%% save
print(ssTBfig_liq,'Fig2A_liq',fmt,rez)
print(ssTBfig_liq,'Fig2A_liq','-dsvg')

print(ssTBfig_ag,'Fig2B_ag',fmt,rez)
print(ssTBfig_ag,'Fig2B_ag','-dsvg')
