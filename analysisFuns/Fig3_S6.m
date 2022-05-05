%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
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

%% load
taus = [-1./(2*log(0:0.1:0.7)),1.8:0.6:8.4];
phis = exp(-1./(2*taus));

doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.environment = 'liquid';
GlobalParams.run_case = 'switch_environments';
GlobalParams.nEnv = 30;

r = GlobalParams.growthRate;
totTimes = (1:0.7:(20*0.7+1))/r;
biophysicalChi = 0;
DA = 800;

GlobalParams.biophysicalChi = biophysicalChi; 
GlobalParams.DA = DA;

GlobalParams.x_discretize_factor = 25;

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;


ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end

[~,ind2] = regexp(ks_dir,[filesep,'ks']);

ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns'];

fsize_slide = [3,2];


%% without diversity
GlobalParams.stdTB = 0; 

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(GlobalParams.stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];

s=load([save_dir,filesep,'figures', filesep, 'outputs_plotSwitchingEnvQuantities_' case_name '.mat'],'cavg_gt');
cavg_gt_nodiv = s.cavg_gt;


%% with diversity
GlobalParams.stdTB = 0.065;

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(GlobalParams.stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];

s=load([save_dir,filesep,'figures', filesep, 'outputs_plotSwitchingEnvQuantities_' case_name '.mat'],'cavg_gt');
cavg_gt = s.cavg_gt;


%% plot

rT_thresh = 4;
mxs = nan(size(totTimes));
mxis = nan(size(totTimes));
for i = 1:length(totTimes)
    if r*totTimes(i)>=rT_thresh
        [mxs(i),mxis(i)] = max(cavg_gt(:,i)./cavg_gt_nodiv(i));
    end
end

mnis = nan(size(taus));
mns = nan(size(taus));
for i = 1:length(taus)
    [mns(i),mnis(i)] = min(cavg_gt(i,:)./cavg_gt_nodiv);
end
mnis(end)=4;
mnis(mnis==3)=4;

figure;
plot(r*totTimes,cavg_gt_nodiv)
xlabel('Environment duration, r T')
ylabel('Migration speed \langlec(t,\sigma=0)\rangle, (\mum/s)')
xlim([min(r*totTimes),max(r*totTimes)])

[T,Phis] = meshgrid(totTimes,phis);


figure;surf(r*T,1./(2*(1-Phis)),(cavg_gt./repmat(cavg_gt_nodiv,length(taus),1)),'LineStyle','none')
ylim([min(1./(2*(1-phis))),max(1./(2*(1-phis)))])
xlim([min(r*totTimes),max(r*totTimes)])
colormap(jet)
h=gca;
view(0,90)
colorbar
hold on
plot3(r*totTimes(r*totTimes>=rT_thresh),1./(2*(1-phis(mxis(r*totTimes>=rT_thresh)))),1.2*mxs(r*totTimes>=rT_thresh),'k--')
plot3(r*totTimes(mnis),1./(2*(1-phis)),1.2*mns,'w')
h.XTick = [1,5,10,15];
h.YTick = [1:2:8];
f=gcf;
f.Units = 'inches';
f.Position(3:4) = [3,2];
print('Fig3main',fmt,rez)
print('Fig3main','-dsvg')
f1=gcf;
h1=gca;
h1.FontSize = 12;


figure;surf(r*T,Phis,(cavg_gt./repmat(cavg_gt_nodiv,length(phis),1)),'LineStyle','none')
ylim([min(phis),max(phis)])
xlim([min(r*totTimes),max(r*totTimes)])
colormap(jet)
h=gca;
view(0,90)
colorbar
hold on
plot3(r*totTimes(r*totTimes>=rT_thresh),phis(mxis(r*totTimes>=rT_thresh)),1.2*mxs(r*totTimes>=rT_thresh),'k--')
plot3(r*totTimes(mnis),phis,1.2*mns,'w')
h.XTick = [1,5,10,15];
h.YTick = 0:0.2:0.8;
h.FontSize = 12;
f=gcf;
f.Units = 'inches';
f.Position(3:4) = [3,2];
print('FigS6C',fmt,rez)
print('FigS6C','-dsvg')


%% side panels
green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];

peakTB = generalistTB;
xupper = 0.45;
yupper = 26;

lw=1.5;
fsize = [1.6,1.2];

% bottom-left
uiopen([ks_dir '\analysisFuns\figures\PTBend_tau=0_rT=1_' case_name '.fig'],1)
h=gca;
f=gcf;
title(h,'');
h.Children(1).Color = colors(2,:);
h.Children(2).Color = colors(1,:);
h.Children(1).LineWidth = lw;
h.Children(2).LineWidth = lw;
h.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
h.Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h,'')
hline = line(peakTB*[1,1],[0,yupper]); hline.LineWidth = lw; hline.LineStyle = '--'; hline.Color = 'k'; hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
h.LineWidth = 1.0;
xlim([0,xupper])
ylim([0,yupper])
xlabel('')
ylabel('')
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h.XTick = 0:0.1:0.4;
h.YTick = [];
print(f,'Fig3botleft',fmt,rez)
print(f,'Fig3botleft','-dsvg')

% bottom right
uiopen([ks_dir '\analysisFuns\figures\PTBend_tau=0_rT=15_' case_name '.fig'],1)
h=gca;
title(h,'');
h.Children(1).Color = colors(2,:);
h.Children(2).Color = colors(1,:);
h.Children(1).LineWidth = lw;
h.Children(2).LineWidth = lw;
h.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
h.Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h,'')
hline = line(peakTB*[1,1],[0,yupper]); hline.LineWidth = lw; hline.LineStyle = '--'; hline.Color = 'k'; hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
h.LineWidth = 1.0;
xlim([0,xupper])
ylim([0,yupper])
xlabel('')
ylabel('')
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h.XTick = 0:0.1:0.4;
h.YTick = [];
print('Fig3botright',fmt,rez)
print('Fig3botright','-dsvg')

% top-left
uiopen([ks_dir '\analysisFuns\figures\PTBend_tau=8.4_rT=1_' case_name '.fig'],1)
h=gca;
f=gcf;
title(h,'');
h.Children(1).Color = colors(2,:);
h.Children(2).Color = colors(1,:);
h.Children(1).LineWidth = lw;
h.Children(2).LineWidth = lw;
h.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
h.Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h,'')
hline = line(peakTB*[1,1],[0,yupper]); hline.LineWidth = lw; hline.LineStyle = '--'; hline.Color = 'k'; hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
h.LineWidth = 1.0;
xlim([0,xupper])
ylim([0,yupper])
xlabel('')
ylabel('')
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h.XTick = 0:0.1:0.4;
h.YTick = [];
print('Fig3topleft',fmt,rez)
print('Fig3topleft','-dsvg')

% top-right
uiopen([ks_dir '\analysisFuns\figures\PTBend_tau=8.4_rT=15_' case_name '.fig'],1)
h=gca;
title(h,'');
h.Children(1).Color = colors(2,:);
h.Children(2).Color = colors(1,:);
h.Children(1).LineWidth = lw;
h.Children(2).LineWidth = lw;
h.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
h.Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h,'')
hline = line(peakTB*[1,1],[0,yupper]); hline.LineWidth = lw; hline.LineStyle = '--'; hline.Color = 'k'; hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
h.LineWidth = 1.0;
xlim([0,xupper])
ylim([0,yupper])
xlabel('')
ylabel('')
f=gcf;
f.Units = 'inches';
f.Position(3:4) = fsize;
h.XTick = 0:0.1:0.4;
h.YTick = [];
print('Fig3topright',fmt,rez)
print('Fig3topright','-dsvg')


%% double the stdTB for SI
% close all

GlobalParams = [];

GlobalParams.stdTB = 2*0.065;

GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.environment = 'liquid';
GlobalParams.run_case = 'switch_environments';
GlobalParams.nEnv = 30;

r = GlobalParams.growthRate;
totTimes = (1:0.7:(20*0.7+1))/r;
biophysicalChi = 0;
DA = 800;

GlobalParams.biophysicalChi = biophysicalChi; 
GlobalParams.DA = DA;

GlobalParams.x_discretize_factor = 25;


generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;

case_name = ['TB=' num2str(GlobalParams.peakTB) '_stdTB=' num2str(GlobalParams.stdTB) '_biophysChi=' num2str(GlobalParams.biophysicalChi) '_DA=' num2str(GlobalParams.DA)];


ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);

ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns'];

s=load([save_dir,filesep,'figures', filesep, 'outputs_plotSwitchingEnvQuantities_' case_name '.mat'],'cavg_gt');
cavg_gt = s.cavg_gt;


%% plot

rT_thresh = 5.2;
mxs = nan(size(totTimes));
mxis = nan(size(totTimes));
for i = 1:length(totTimes)
    if r*totTimes(i)>=rT_thresh
        [mxs(i),mxis(i)] = max(cavg_gt(:,i)./cavg_gt_nodiv(i));
    end
end

mnis = nan(size(taus));
mns = nan(size(taus));
for i = 1:length(taus)
    [mns(i),mnis(i)] = min(cavg_gt(i,:)./cavg_gt_nodiv);
end
mnis(end)=4;
mnis(mnis==3)=4;

figure;
plot(r*totTimes,cavg_gt_nodiv)
xlabel('Environment duration, r T')
ylabel('Migration speed \langlec(t,\sigma=0)\rangle, (\mum/s)')
xlim([min(r*totTimes),max(r*totTimes)])


figure;surf(r*T,1./(2*(1-Phis)),(cavg_gt./repmat(cavg_gt_nodiv,length(taus),1)),'LineStyle','none')
ylim([min(1./(2*(1-phis))),max(1./(2*(1-phis)))])
xlim([min(r*totTimes),max(r*totTimes)])
colormap(jet)
h=gca;
view(0,90)
colorbar
hold on
plot3(r*totTimes(r*totTimes>=rT_thresh),1./(2*(1-phis(mxis(r*totTimes>=rT_thresh)))),1.2*mxs(r*totTimes>=rT_thresh),'k--')
plot3(r*totTimes(mnis),1./(2*(1-phis)),1.2*mns,'w')
h.XTick = [1,5,10,15];
h.YTick = [1:2:8];
h.FontSize = 12;
f=gcf;
f.Units = 'inches';
f.Position(3:4) = [3,2];
print('FigS6A',fmt,rez)
print('FigS6A','-dsvg')


figure;surf(r*T,Phis,(cavg_gt./repmat(cavg_gt_nodiv,length(phis),1)),'LineStyle','none')
ylim([min(phis),max(phis)])
xlim([min(r*totTimes),max(r*totTimes)])
colormap(jet)
h=gca;
view(0,90)
colorbar
hold on
plot3(r*totTimes(r*totTimes>=rT_thresh),phis(mxis(r*totTimes>=rT_thresh)),1.2*mxs(r*totTimes>=rT_thresh),'k--')
plot3(r*totTimes(mnis),phis,1.2*mns,'w')
h.XTick = [1,5,10,15];
h.YTick = 0:0.2:0.8;
h.FontSize = 12;
f=gcf;
f.Units = 'inches';
f.Position(3:4) = [3,2];
print('FigS6D',fmt,rez)
print('FigS6D','-dsvg')

c = h.CLim;


h1.CLim = c;
print(f1,'FigS6B',fmt,rez)
print(f1,'FigS6B','-dsvg')
