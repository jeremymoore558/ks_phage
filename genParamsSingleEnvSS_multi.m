%% constant across simulations
clear GlobalParams

addpath(['.',filesep,'functions']);

doubling_time = 50; % min
GlobalParams.growthRate = 1/doubling_time*1/60*log(2); % 1/s
GlobalParams.asp = 100; % uM
GlobalParams.reSimulate = 0;

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

SimParams = [];
Envs = {'liquid','agar'}; 

%% varying parameters
envPlaceHolder = [1,2];
peakTBs = generalistTB;
phis = exp(-1./(2*[-1./(2*log(0:0.1:0.7)),1.8:0.6:6]));
stdTBs = [0,0.065,2*0.065];
biophysChis = 0;
DAs = 800; 

[PeakTBs,Phis,StdTBs,EnvPlaceHolder,DAmat,BiophysChis] = ndgrid(peakTBs,phis,stdTBs,envPlaceHolder,DAs,biophysChis);

% remove redundant simulation cases
rmInds = find(StdTBs==0 & Phis>0);

PeakTBs(rmInds) = [];
Phis(rmInds) = [];
StdTBs(rmInds) = [];
EnvPlaceHolder(rmInds) = [];
DAmat(rmInds) = [];
BiophysChis(rmInds) = [];

rmInds = find(DAmat==0 & BiophysChis==1);
PeakTBs(rmInds) = [];
Phis(rmInds) = [];
StdTBs(rmInds) = [];
EnvPlaceHolder(rmInds) = [];
DAmat(rmInds) = [];
BiophysChis(rmInds) = [];

TotTime = 15*60*60; % s
%
OD0s = nan(size(Phis));
OD0s(EnvPlaceHolder==1) = 3;
OD0s(EnvPlaceHolder==2) = 0.5; 

nSim = numel(StdTBs);
nSim

%% run simulations that have not been run -- will automatically skip those that are done
% recommend running these in parallel on a cluster
for i = 1:nSim
  GlobalParams_i = GlobalParams;
  GlobalParams_i.peakTB = PeakTBs(i);
  GlobalParams_i.phi = Phis(i);
  GlobalParams_i.stdTB = StdTBs(i);
  GlobalParams_i.totTime = TotTime;
  GlobalParams_i.DA = DAmat(i);
  GlobalParams_i.biophysicalChi = BiophysChis(i);
  GlobalParams_i.environment = Envs{EnvPlaceHolder(i)};
  GlobalParams_i.OD0 = OD0s(i);
 
  SimParams_i = SimParams;

  mainFunction(GlobalParams_i,SimParams_i);
end