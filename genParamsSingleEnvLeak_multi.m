%% constant across simulations
clear GlobalParams

addpath(['.',filesep,'functions']);

GlobalParams = [];
GlobalParams.growthRate = 0; 
GlobalParams.asp = 100; % uM
GlobalParams.totTime = 10*60*60; 
GlobalParams.environment = 'liquid';
GlobalParams.x_discretize_factor = 25;

GlobalParams.OD0 = 6;
GlobalParams.reSimulate = 0;
SimParams = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams.meanTB = generalistTB;

%% varying parameters
Envs = {'liquid','agar'};
envPlaceHolder = [1,2];
stdTBs = [0, 0.065]; 
biophysicalChis = 0;
Das = 800; % um^2/s

[Stds,BiophysChis,DAs,EnvPlaceHolder] = ndgrid(stdTBs,biophysicalChis,Das,envPlaceHolder);

nSim = numel(DAs);
nSim

%% run simulations that have not been run -- will automatically skip those that are done
% recommend running these in parallel on a cluster
runInds = ones(nSim,1);
for i = 1:nSim
   GlobalParams_i = GlobalParams;
   GlobalParams_i.stdTB = Stds(i);
   GlobalParams_i.environment = Envs{EnvPlaceHolder(i)};
   GlobalParams_i.DA = DAs(i);
   GlobalParams_i.biophysicalChi = BiophysChis(i);
   SimParams_i = SimParams;
   
   runInds(i) = mainFunction(GlobalParams_i,SimParams_i);
end
