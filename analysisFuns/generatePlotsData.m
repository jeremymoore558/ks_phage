% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',20);
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

%%
ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir.path;
end

[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

addpath([ks_dir,filesep,'functions'])


%% steady state quantities
disp('Generating steady state plot data.')

stdTB = 0.065;
taus = [-1./(2*log(0:0.1:0.7)),1.8:0.6:6];
phis = exp(-1./(2*taus));
r = 1/50*log(2)/60;
plotSteadyStateQuantities(stdTB,phis,0,800); % stdTB, biophysChi, DA
plotSteadyStateQuantities(2*stdTB,phis,0,800); % stdTB, biophysChi, DA
% these run the no-diversity cases as well

%% switching environments, fixed durations
disp('Generating switching environments plot data.')

r = 1/50*log(2)/60;
taus = [-1./(2*log(0:0.1:0.7)),1.8:0.6:8.4];
phis = exp(-1./(2*taus));
totTimes = (1:0.7:(20*0.7+1))/r;

% diversity
disp('Diverse population...')
stdTB = 0.065;
plotSwitchingEnvQuantities(stdTB,phis,totTimes,0,800,25); % stdTB, biophysChi, DA = 800, x_discretization_factor
plotSwitchingEnvQuantities(2*stdTB,phis,totTimes,0,800,25); % stdTB, biophysChi, DA = 800, x_discretize_factor

% no diversity
disp('No diversity...')
stdTB = 0;
plotSwitchingEnvQuantities(stdTB,0,totTimes,0,800,25); % stdTB=0, biophysChi, DA = 800 um^2/s, x_discretize_factor

%%
disp('Done!.')
