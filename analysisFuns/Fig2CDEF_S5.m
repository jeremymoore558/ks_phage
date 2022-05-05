%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultLineMarkerSize',3);
set(0,'DefaultLineColor',[0 0 0]);
set(0,'DefaultAxesBox','off');
set(0,'DefaultAxesColor','none');
set(0,'DefaultAxesLineWidth',1.0);
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultTextInterpreter','tex');
set(0,'DefaultAxesTickLabelInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex');
whitebg('white');
close


fsize = [1.7,1.29];
fsize_slide = [4,3];

%% global
% close all

stdTB = 0.065;
taus = [-1./(2*log(0:0.1:0.7)),1.8:0.6:6]; % 1.2
phis = exp(-1./(2*taus));

biophysicalChi = 0;
DA = 800;

GlobalParams = [];
doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.runTill = 'time';
GlobalParams.stdTB = stdTB;
GlobalParams.reSimulate = 0;
GlobalParams.biophysicalChi = biophysicalChi;
GlobalParams.DA = DA;
SimParams = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;

[GlobalParams, SimParams] = initializeSimulationStructures(GlobalParams, SimParams);

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

[ks_dir, filesep, 'functions'];
addpath([ks_dir, filesep, 'functions'])

save_dir = [ks_dir,filesep,'analysisFuns'];

dense_factor = 10;

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];

load([save_dir, filesep, 'figures', filesep, 'outputs_plotSteadyStateQuantities_' case_name '.mat'])


%% liquid relaxation time

t_adapt_liq = nan(1,length(taus));
t_adapt_se_liq = nan(1,length(taus));

for g = 1:length(taus)
    dt_global = min(diff(unique(t_liq_g{g})));
    dt_global = max(dt_global,60);
    t_global = min(t_liq_g{g}):dt_global:max(t_liq_g{g});
    [~,inds] = unique(t_liq_g{g});
    ct_global = interp1(t_liq_g{g}(inds),c_liq_g{g}(inds),t_global);
    
    figure;hold on
    plot(t_global,ct_global)
    
    t1 = find(abs(ct_global-ct_global(end))/abs(ct_global(end)-ct_global(1))>0.1,1,'last'); % 0.05
    
    f = fit(t_global(t1:end)',ct_global(t1:end)',@(b,c,d,x) b+c*exp(-x/d),'StartPoint',[mean(ct_global(end-10:end)),-5,-3e3],'Lower',[0;-20;1e2],'Upper',[20;0;1e5]);
    
    plot(t_global,f(t_global),'--')
    plot(t_global(t1:end),f(t_global(t1:end)),'--')
    
    ci = confint(f,0.68);
    
    t_adapt_liq(g) = f.d;
    t_adapt_se_liq(g) = 1/2*abs(diff(ci(:,2)));
    
    
end



%% agar relaxation time
t_adapt_agar = nan(1,length(taus));
t_adapt_se_agar = nan(1,length(taus));
for g = 1:length(taus)
    dt_global = min(diff(unique(t_agar_g{g})));
    t_global = min(t_agar_g{g}):dt_global:max(t_agar_g{g});
    [~,inds] = unique(t_agar_g{g});
    ct_global = interp1(t_agar_g{g}(inds),c_agar_g{g}(inds),t_global);

    [~,loc] = findpeaks(-c_agar_g{g},'MinPeakProminence',0.005);
    [~,t2] = min(abs(t_agar_g{g}-1.5e5));
    
    t1 = find(abs(ct_global-ct_global(end))/abs(ct_global(end)-ct_global(1))>0.1,1,'last'); % 0.05
    f = fit(t_global(t1:end)',ct_global(t1:end)',@(a,b,c,x) c+a*exp(-x/b),'StartPoint',[-5,-3e3,mean(ct_global(end-10:end))],'Lower',[-20;1e2;0],'Upper',[0;1e5;20]);
    
    
    figure;hold on
    
    plot(t_agar_g{g},c_agar_g{g})
    plot(t_global(t1:end),ct_global(t1:end))
    plot(t_global(t1:end),f(t_global(t1:end)),'--')
    ylim([0.6,1.15])
    
    ci = confint(f,0.68);
  
    t_adapt_agar(g) = f.b;
    t_adapt_se_agar(g) = 1/2*abs(diff(ci(:,2)));
    
end

%% some parameters and quantities
close all

ms = 2;

r = GlobalParams.growthRate;
ep = GlobalParams.KiA/GlobalParams.asp;

green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];

k_adapt_liq = 1./t_adapt_liq;
k_adapt_agar = 1./t_adapt_agar;

% predicted ss speed
css_liq_g_pred = sqrt((chiBack_liq_g-muBack_liq_g).^2./(chiBack_liq_g-muBack_liq_g+GlobalParams.DA));
css_agar_g_pred = sqrt((chiBack_agar_g-muBack_agar_g).^2./(chiBack_agar_g-muBack_agar_g+GlobalParams.DA));

% predicted speed using population avg chi-mu at each time point vs simulation speed
css_liq_g_pred_total_pop = sqrt((chi_liq_g-mu_liq_g).^2./(chi_liq_g-mu_liq_g+GlobalParams.DA));
css_agar_g_pred_total_pop = sqrt((chi_agar_g-mu_agar_g).^2./(chi_agar_g-mu_agar_g+GlobalParams.DA));

% trade-offs in non-diverse population speeds
TBs = logspace(-2,log10(0.4),50);

[~,mxi]=max(min(css_ag_TB/max(css_ag_TB),css_liq_TB/max(css_liq_TB)))
TBs(mxi)

[~,mxiA]=max(css_ag_TB)
TBs(mxiA)

[~,mxiL]=max(css_liq_TB)
TBs(mxiL)

css_liq_TB(1:4) = nan;
css_ag_TB(1:4) = nan;
TBs(1:4) = nan;


% batch culture composition relaxation time
krelax = (2*r*(1-phis));


%%

% plot predictions against simulations
figure;hold on
plot(css_liq_g/css_liq_g(1),css_liq_g_pred/css_liq_g_pred(1),'Color',colors(1,:),'Marker','o','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',ms)
plot(css_agar_g/css_agar_g(1),css_agar_g_pred/css_agar_g_pred(1),'Color',colors(2,:),'Marker','o','MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerSize',ms)

plot(css_liq_g/css_liq_g(1),css_liq_g_pred_total_pop/css_liq_g_pred_total_pop(1),'-','Color',colors(1,:))
plot(css_agar_g/css_agar_g(1),css_agar_g_pred_total_pop/css_agar_g_pred_total_pop(1),'-','Color',colors(2,:))

hline = line([1,2],[1,2]); hline.Color = 'k'; hline.LineStyle = '--';
hline.LineWidth = 1;
xlim([1,1.22])
ylim([1,1.22])
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h=gca;
h.FontSize = 10;
h.XTick = 1:.1:1.2;
h.YTick = 1:.1:1.2;
h.XTickLabelRotation = 0;
print('Fig2D',fmt,rez)
print('Fig2D','-dsvg')



% steady state speed vs mother-daughter correlation
figure;hold on

plot(phis,css_liq_g/css_liq_g(1),'Color',colors(1,:),'Marker','o','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',ms)
plot(phis,css_agar_g/css_agar_g(1),'Color',colors(2,:),'Marker','o','MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerSize',ms)

f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h=gca;
h.XTick = 0:0.2:1;
h.YTick = 1:0.1:1.2;
h.FontSize = 10;

print('Fig2C',fmt,rez)
print('Fig2C','-dsvg')


% trade-offs - Fig 2E
phis_plot = 0:0.2:0.8;
cssliqnorm = interp1(phis,css_liq_g/css_liq_liqspec,phis_plot);
cssagnorm = interp1(phis,css_agar_g/css_ag_agspec,phis_plot);

figure;hold on
hp = plot(css_liq_TB/max(css_liq_TB),css_ag_TB/max(css_ag_TB),'k');
hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([0,1])
ylim([0,1])

plot(css_liq_TB(mxiL)/max(css_liq_TB),css_ag_TB(mxiL)/max(css_ag_TB),'o','Color','k','MarkerFaceColor',colors(1,:),'MarkerSize',ms+2,'LineWidth',0.5)
plot(css_liq_TB(mxiA)/max(css_liq_TB),css_ag_TB(mxiA)/max(css_ag_TB),'o','Color','k','MarkerFaceColor',colors(2,:),'MarkerSize',ms+2,'LineWidth',0.5)
plot(css_liq_TB(mxi)/max(css_liq_TB),css_ag_TB(mxi)/max(css_ag_TB),'ko','MarkerFaceColor','k','MarkerSize',ms+1)

jetcolors=jet(2*length(phis_plot));
for i = 1:length(phis_plot)
    hp3 = plot(cssliqnorm(i),cssagnorm(i),'o','Color','k','MarkerFaceColor',jetcolors(2*i,:),'MarkerSize',ms+2,'LineWidth',0.5);
end


h=gca;
h.XTick = 0:0.2:1;
h.YTick = 0:0.2:1;
h.FontSize = 10;
xlim([0.4,1])

f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h=gca;
h.Units = 'inches';
pos = h.Position;
h.Position = pos;

print('Fig2E_v2',fmt,rez)
print('Fig2E_v2','-dsvg')


% speed/adaptation time trade-off
figure;hold on
plot((css_liq_g-css_liq_g(1))/range(css_liq_g),(k_adapt_liq-k_adapt_liq(end))/range(k_adapt_liq),'Color',colors(1,:),'Marker','o','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',ms)
plot((css_agar_g-css_agar_g(1))/range(css_agar_g),(k_adapt_agar-k_adapt_agar(end))/range(k_adapt_agar),'Color',colors(2,:),'Marker','o','MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerSize',ms)

h=gca;
h.XTick = 0:0.2:1;
h.YTick = 0:0.2:1;
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h.FontSize = 10;

print('Fig2F',fmt,rez)
print('Fig2F','-dsvg')



%% SI
set(0,'DefaultAxesLineWidth',1.5);

% adaptation time vs composition relaxation time -- SI figure
figure;hold on
plot(krelax/r,k_adapt_liq/r,'Color',colors(1,:))
plot(krelax/r,k_adapt_agar/r,'Color',colors(2,:))

h = gca;
hline = line([0,2],[0,2]); hline.Color = 'k'; hline.LineStyle = '--';
h.XLim = [0,2];

h.FontSize = 12;
f=gcf;
f.Units = 'inches';
f.Position(3:4) = [2.67,2];
h.XTick = [0:0.4:2];
h.YTick = [0:0.4:2];

print('FigS5',fmt,rez)
print('FigS5','-dsvg')

