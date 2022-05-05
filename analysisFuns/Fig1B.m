%% chi and mu in diff envs

%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
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

%%
ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns'];

biophysicalChi = 0;
GlobalParams.biophysicalChi = biophysicalChi;
case_name = ['biophysChi=' num2str(GlobalParams.biophysicalChi)];

green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];

fsize = [1.7,1.29];

%%
GlobalParams.growthRate = 0;%log(2)/50/60; % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.runTill = 'time';
GlobalParams.totTime = 5*60*60; % s
GlobalParams.environment = 'liquid'; % or 'agar'
GlobalParams.peakTB = 0.14; %.124
GlobalParams.stdTB = 0.065;
SimParams = [];
[GlobalParams,SimParams] = initializeSimulationStructures(GlobalParams,SimParams);

MFTliq = initializeMFTmodel(GlobalParams.F0,'liquid',biophysicalChi);
MFTagar = initializeMFTmodel(GlobalParams.F0,'agar',biophysicalChi);


%%
F0 = linspace(-10,20,1000);
MFTliq = initializeMFTmodel(F0,'liquid',biophysicalChi);
MFTagar = initializeMFTmodel(F0,'agar',biophysicalChi);
TB = 1./(1+exp(F0));


%% plot

figure;hold on
plot(TB,MFTliq.chi,'Color',colors(1,:),'LineWidth',1.5)
plot(TB,MFTagar.chi,'Color',colors(2,:),'LineWidth',1.5)
h=gca;
xlim([0,1])

h.YScale = 'log';
ylim([10,10^4]);
h.FontSize = 10;

% for Fig1A
% pos = h.Position
% ypos = h.YLabel.Position
% xpos = h.XLabel.Position


f=gcf;
f.Units = 'inches';
dx = 0.4;

f.Position(3:4) = fsize;

h.Units = 'inches';

hleg = legend('Liquid','Porous','Location','northeast');
hleg.Box = 'off';
hleg.FontSize = 9;

print('Fig1B',fmt,rez)
print('Fig1B','-dsvg')
