% density profiles, gradient

%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineMarkerSize',3);
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

fsize1 = [2,1.5];
fsize2 = [1.5,1.13];


%% 
dense_factor = 10;

stdTB = 0.065;
biophysicalChi = 0;
DA = 800;
OD0 = 6;
KiA = 3.5; % uM
KmA = 0.5;


GlobalParams_all.growthRate = 0;
GlobalParams_all.asp = 100; % uM
GlobalParams_all.stdTB = stdTB;
GlobalParams_all.reSimulate = 0;
GlobalParams_all.biophysicalChi = biophysicalChi; % this or the case that can be compared to theory?
GlobalParams_all.DA = DA;
GlobalParams_all.OD0 = OD0;
GlobalParams_all.KiA = KiA;
GlobalParams_all.KmA = KmA;
GlobalParams_all.IC = 'step'; % 'sigmoid','analytical'
SimParams_all = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir(1).path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns'];

ODumToCellsPercm2 = 8.4e8*(10^-2/10^-6); % OD * um * cells/cm^3/OD * cm/um

TBupper = 0.4;

case_name = ['_biophysChi=' num2str(GlobalParams_all.biophysicalChi) '_DA=' num2str(GlobalParams_all.DA) '_KiA=' num2str(GlobalParams_all.KiA)];

if ~isempty(KmA)
     case_name = [case_name, '_KmA=' num2str(GlobalParams_all.KmA) '_50x'];
end

green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];


%% liquid calculations
environment = 'liquid'; %'agar' % run both

GlobalParams_all.environment = environment;

[GlobalParams,SimParams] = initializeSimulationStructures(GlobalParams_all,SimParams_all);
ep = GlobalParams.KiA/GlobalParams.asp;

load(GlobalParams.simName)

% map indices within each block to indices if everything were all in one
% vector
blockInds = cell(1,length(SimParams));
for i = 1:length(SimParams)
    if i==1
        blockInds{i} = 1:length(SimParams(i).time);
    else
        blockInds{i} = blockInds{i-1}(end) + (1:length(SimParams(i).time));
    end
end

i = length(SimParams);
j = length(SimParams(i).time);

tInd = blockInds{i}(j);


F0 = SimParams(i).F0;
TB = 1./(1+exp(F0));

x = SimParams(i).x;
dx = SimParams(i).dx;

rhoF = squeeze(SimResults(i).rho(:,:,j));
asp = SimResults(i).asp(:,j);
dTBdF0 = TB.*(1-TB);

[~,mxInd] = max(sum(rhoF,1));
zsim = x-x(mxInd);

MFT = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);
chi = MFT.chi;
mu = MFT.mu;

backInd = find(asp<sqrt(KiA*KmA),1,'last');
NF0 = dx*sum(rhoF(:,backInd:end),2);

PF0 = NF0/trapz(F0,NF0,1);
PTB = PF0./dTBdF0;
PTB = PTB/(-trapz(TB,PTB));

meanChi = trapz(F0,PF0.*chi);
meanMu = trapz(F0,PF0.*mu);
c = sum(NF0,1)*GlobalParams.kA/GlobalParams.asp; % units are OD*um *uM/OD/s / uM = um/s...
lambda = c/(meanChi-meanMu);


% zmax(TB)
zmax_sim = nan(size(rhoF,1),1);
for k = 1:size(rhoF,1)
    [~,locs] = findpeaks(rhoF(k,backInd:end),'MinPeakHeight',1e-5);
    
    if ~isempty(locs)
        locs = locs(end)+backInd-1;
        zmax_sim(k) = zsim(locs);
    end
end


% <TB|z> in simulation
meanTB_z_sim = sum(repmat(TB(:),1,length(zsim)).*rhoF,1)./sum(rhoF,1);

% <\chi-\mu|z> in simulation
meanChiMu_z_sim = sum(repmat(chi(:)-mu(:),1,length(zsim)).*rhoF,1)./sum(rhoF,1);

% gradient
f_sim = log(1+asp/GlobalParams.KiA);
df_sim = (f_sim(3:end)-f_sim(1:end-2))/(2*dx);
[~,mxif] = max(df_sim);

% gradient at the back
dfback_sim = c./meanChiMu_z_sim;


%% plots

f = figure; hold on
plot(zsim/1000,sum(rhoF,1))
xlim([-11,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
print('FigS1A_rho_z_liq','-dsvg')
print('FigS1A_rho_z_liq',fmt,rez)

f = figure; hold on
plot(zsim/1000,asp/GlobalParams.asp)
xlim([-11,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
ylim([0,1])
print('FigS1B_s_z_liq','-dsvg')
print('FigS1B_s_z_liq',fmt,rez)

f = figure; hold on
plot(zsim(2:end-1)/1000,df_sim*1000)
xlim([-11,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
print('FigS1C_df_z_liq','-dsvg')
print('FigS1C_df_z_liq',fmt,rez)


f = figure; hold on
plot(chi,zmax_sim/1000)
xlim([3000,5600]);
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
print('FigS2A_zmax_chi_liq','-dsvg')
print('FigS2A_zmax_chi_liq',fmt,rez)


f = figure; hold on
plot(TB,zmax_sim/1000)
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
print('FigS2B_zmax_TB_liq','-dsvg')
print('FigS2B_zmax_TB_liq',fmt,rez)


f = figure; hold on
plot(zsim/1000,meanTB_z_sim)
xlim([-11,5])
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
h.Units = 'normalized';
h.Position(1)=.21;
h.YTick = 0.06:0.02:0.12;
print('FigS2C_meanTB_z_liq','-dsvg')
print('FigS2C_meanTB_z_liq',fmt,rez)

f = figure; hold on
plot(zsim/1000,meanChiMu_z_sim)
xlim([-11,5])
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
h.Units = 'normalized';
h.Position(1)=0.24;
print('FigS2D_meanChiMu_z_liq','-dsvg')
print('FigS2D_meanChiMu_z_liq',fmt,rez)



%% agar calculations
environment = 'agar'; %'agar' % run both

GlobalParams_all.environment = environment;

[GlobalParams,SimParams] = initializeSimulationStructures(GlobalParams_all,SimParams_all);
ep = GlobalParams.KiA/GlobalParams.asp;

load(GlobalParams.simName)


F0 = SimParams(end).F0;
TB = 1./(1+exp(F0));

x = SimParams(end).x;
dx = SimParams(end).dx;

rhoF = squeeze(SimResults(end).rho(:,:,end));
asp = SimResults(end).asp(:,end);
dTBdF0 = TB.*(1-TB);

ind1 = find(x/1000>10,1,'first');
[~,mxInd] = max(sum(rhoF(:,ind1:end),1));
mxInd = mxInd+ind1-1;
zsim = x-x(mxInd);

MFT = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);
chi = MFT.chi;
mu = MFT.mu;

backInd = find(asp<sqrt(KiA*KmA),1,'last');
NF0 = dx*sum(rhoF(:,backInd:end),2);

PF0 = NF0/trapz(F0,NF0,1);
PTB = PF0./dTBdF0;
PTB = PTB/(-trapz(TB,PTB));

meanChi = trapz(F0,PF0.*chi);
meanMu = trapz(F0,PF0.*mu);
c = sum(NF0,1)*GlobalParams.kA/GlobalParams.asp; % units are OD*um *uM/OD/s / uM = um/s...
lambda = c/(meanChi-meanMu);

% zmax(TB)
zmax_sim = nan(size(rhoF,1),1);
for k = 1:size(rhoF,1)
    [~,locs] = findpeaks(rhoF(k,backInd:end),'MinPeakHeight',1e-5);
    
    if ~isempty(locs)
        locs = locs(end)+backInd-1;
        zmax_sim(k) = zsim(locs);
    end
end

% <TB|z> in simulation
meanTB_z_sim = sum(repmat(TB(:),1,length(zsim)).*rhoF,1)./sum(rhoF,1);

% <\chi-\mu|z> in simulation
meanChiMu_z_sim = sum(repmat(chi(:)-mu(:),1,length(zsim)).*rhoF,1)./sum(rhoF,1);

% gradient
f_sim = log(1+asp/GlobalParams.KiA);
df_sim = (f_sim(3:end)-f_sim(1:end-2))/(2*dx);

% gradient at the back
dfback_sim = c./meanChiMu_z_sim;

%% plots

f = figure; hold on
plot(zsim/1000,sum(rhoF,1))
xlim([-5,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
h.Units = 'normalized';
h.Position(1)=0.17;
print('FigS1A_rho_z_ag','-dsvg')
print('FigS1A_rho_z_ag',fmt,rez)

f = figure; hold on
plot(zsim/1000,asp/GlobalParams.asp)
xlim([-5,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
ylim([0,1])
print('FigS1B_s_z_ag','-dsvg')
print('FigS1B_s_z_ag',fmt,rez)

f = figure; hold on
plot(zsim(2:end-1)/1000,df_sim*1000)
xlim([-5,5])
f.Units = 'inches';
f.Position(3:4) = fsize1;
h=gca;h.Box = 'off';
print('FigS1C_df_z_ag','-dsvg')
print('FigS1C_df_z_ag',fmt,rez)


f = figure; hold on
plot(chi,zmax_sim/1000)
xlim([275,475])
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
h.XTick = 275:100:475;
print('FigS2A_zmax_chi_ag','-dsvg')
print('FigS2A_zmax_chi_ag',fmt,rez)


f = figure; hold on
plot(TB,zmax_sim/1000)
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
h.XTick = 0:0.25:0.5;
print('FigS2B_zmax_TB_ag','-dsvg')
print('FigS2B_zmax_TB_ag',fmt,rez)


f = figure; hold on
plot(zsim/1000,meanTB_z_sim)
xlim([-5,5])
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
h.YTick = 0.18:0.02:0.24;
print('FigS2C_meanTB_z_ag','-dsvg')
print('FigS2C_meanTB_z_ag',fmt,rez)

f = figure; hold on
plot(zsim/1000,meanChiMu_z_sim)
xlim([-5,5])
f.Units = 'inches';
f.Position(3:4) = fsize2;
h=gca;h.Box = 'off';
print('FigS2D_meanChiMu_z_ag','-dsvg')
print('FigS2D_meanChiMu_z_ag',fmt,rez)