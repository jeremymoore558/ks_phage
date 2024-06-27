%% Main function. Wrapper to run multiple simulations in parallel
%This script has no imputs. It just runs simulations by calling simulate
%wave, a function which takes in the struct SimParams. Various parameters
%can be edited after calling parameters(). 

%This wrapper is different from Main.m because it calls
%simulateWave_GaussPhage.m, which uses a different initial phage
%distribution (gaussian instead of a flat landscape.

clear 
close all

%% I vs. Chi
%List of parameter values to try
%Values to scale cell2 I and Chi by.
% cell2I = [0.01, 0.1, 1, 2, 5];
cell2I = [0.01, 0.1, 1, 10, 100];
cell2Chi = [0.18, 0.42, 1, 2.37, 5.62];
% cell2I = [1];
% cell2Chi = [1];

tic()
parfor ii = 1:length(cell2I)
    for jj = 1:length(cell2Chi)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.irate2 = SimParams.irate2 .* cell2I(ii);
        SimParams.Chi2 = SimParams.Chi2 .* cell2Chi(jj);
        SimParams.OutFolderName = './Outputs_WideGaussPhage/OutputsIChi_GaussPhage/';
        
        %Run simulation
        simulateWave_WideGaussPhage(SimParams)
    end
end
toc()

%% I vs. cA
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2I = [0.01, 0.1, 1, 10, 100];
cell2cA = [0.44, 0.67, 1, 1.5, 2.25];


tic()
parfor ii = 1:length(cell2I)
    for jj = 1:length(cell2cA)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.irate2 = SimParams.irate2 .* cell2I(ii);
        SimParams.cA2 = SimParams.cA2 .* cell2cA(jj);
        SimParams.OutFolderName = './Outputs_WideGaussPhage/OutputsIcA_GaussPhage/';

        %Run simulation
        simulateWave_WideGaussPhage(SimParams)
    end
end
toc()

%% I vs. cR
%List of parameter values to try
cell2I = [0.01, 0.1, 1, 10, 100];
%Values to scale cell2 I and Chi by.
cell2cR = [0.44, 0.67, 1, 1.5, 2.25];


tic()
parfor ii = 1:length(cell2I)
    for jj = 1:length(cell2cR)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.irate2 = SimParams.irate2 .* cell2I(ii);
        SimParams.cR2 = SimParams.cR2 .* cell2cR(jj);
        SimParams.OutFolderName = './Outputs_WideGaussPhage/OutputsIcR_GaussPhage/';

        %Run simulation
        simulateWave_WideGaussPhage(SimParams)
    end
end
toc()

%% I vs. Y
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2I = [0.01, 0.1, 1, 10, 100];
cell2Y = [0.71, 0.84, 1, 1.19, 1.41];


tic()
parfor ii = 1:length(cell2I)
    for jj = 1:length(cell2Y)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.irate2 = SimParams.irate2 .* cell2I(ii);
        SimParams.Y2 = SimParams.Y2 .* cell2Y(jj);
        SimParams.OutFolderName = './Outputs_WideGaussPhage/OutputsIY_GaussPhage/';

        %Run simulation
        simulateWave_WideGaussPhage(SimParams)
    end
end
toc()
print('All Experiments are Done')
