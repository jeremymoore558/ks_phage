%% constant across simulations
clear GlobalParams

addpath(['.',filesep,'functions']);

doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.environment = 'liquid';
GlobalParams.reSimulate = 0;
GlobalParams.run_case = 'switch_environments';
GlobalParams.nEnv = 30;
GlobalParams.x_discretize_factor = 25;
SimParams = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.peakTB = generalistTB;

%% varying parameters
phis = exp(-1./(2*[-1./(2*log(0:0.1:0.7)),1.8:0.6:8.4]));
stdTBs = [0,0.065,2*0.065];
biophysChis = 0;
DAs = 800;
totTimes = (1:0.7:(20*0.7+1))/GlobalParams.growthRate;

[Phis,StdTBs,TotTimes,DAmat,BiophysChis] = ndgrid(phis,stdTBs,totTimes,DAs,biophysChis);

% remove redundant simulation cases
rmInds = find(StdTBs==0 & Phis>0);

Phis(rmInds) = [];
StdTBs(rmInds) = [];
TotTimes(rmInds) = [];
DAmat(rmInds) = [];
BiophysChis(rmInds) = [];

rmInds = find(DAmat==0 & BiophysChis==1);

Phis(rmInds) = [];
StdTBs(rmInds) = [];
TotTimes(rmInds) = [];
DAmat(rmInds) = [];
BiophysChis(rmInds) = [];

nSim = numel(StdTBs);
nSim

%% check whether simulation has been run to completion
% recommend running these in parallel on a cluster
for i = 1:nSim
  GlobalParams_i = GlobalParams;
  GlobalParams_i.phi = Phis(i);
  GlobalParams_i.stdTB = StdTBs(i);
  GlobalParams_i.totTime = TotTimes(i);
  GlobalParams_i.DA = DAmat(i);
  GlobalParams_i.biophysicalChi = BiophysChis(i);
 
  SimParams_i = SimParams;

  mainFunction(GlobalParams_i,SimParams_i);
end