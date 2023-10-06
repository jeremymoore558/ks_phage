%% Main function. Wrapper to run multiple simulations in parallel
%This script has no imputs. It just runs simulations by calling simulate
%wave, a function which takes in the struct SimParams. Various parameters
%can be edited after calling parameters(). 

clear 
close all

%% b vs. Chi
%List of parameter values to try
%Values to scale cell2 b and Chi by.
cell2b = [10^(-1), 10^-0.5, 1, 10^0.5, 10];
cell2Chi = [0.18, 0.42, 1, 2.37, 5.62];

tic()
parfor ii = 1:length(cell2b)
    for jj = 1:length(cell2Chi)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.Chi2 = SimParams.Chi2 .* cell2Chi(ii);
        SimParams.b2 = SimParams.b2 .* cell2b(jj);
        SimParams.OutFolderName = './Outputs_b_NoPhage/OutputsbChi/';
        SimParams.PFU = 0;	
        
        %Run simulation
        simulateWave(SimParams)
    end
end
toc()

%% b vs. cA
%List of parameter values to try
%Values to scale cell2 b and Chi by.
cell2b = [10^-1, 10^-0.5, 1, 10^0.5, 10];
cell2cA = [0.44, 0.67, 1, 1.5, 2.25];


tic()
parfor ii = 1:length(cell2b)
    for jj = 1:length(cell2cA)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.b2 = SimParams.b2 .* cell2b(ii);
        SimParams.cA2 = SimParams.cA2 .* cell2cA(jj);
        SimParams.OutFolderName = './Outputs_b_NoPhage/OutputsbcA/';
        SimParams.PFU = 0;	

        %Run simulation
        simulateWave(SimParams)
    end
end
toc()

%% b vs. cR
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2b = [10^-1, 10^-0.5, 1, 10^0.5, 10];
%Values to scale cell2 I and Chi by.
cell2cR = [0.44, 0.67, 1, 1.5, 2.25];


tic()
parfor ii = 1:length(cell2b)
    for jj = 1:length(cell2cR)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.b2 = SimParams.b2 .* cell2b(ii);
        SimParams.cR2 = SimParams.cR2 .* cell2cR(jj);
        SimParams.OutFolderName = './Outputs_b_NoPhage/OutputsbcR/';
        SimParams.PFU = 0;	

        %Run simulation
        simulateWave(SimParams)
    end
end
toc()

%% b vs. Y
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2b = [10^-1, 10^-0.5, 1, 10^0.5, 10];
cell2Y = [0.71, 0.84, 1, 1.19, 1.41];


tic()
parfor ii = 1:length(cell2b)
    for jj = 1:length(cell2Y)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.b2 = SimParams.b2 .* cell2b(ii);
        SimParams.Y2 = SimParams.Y2 .* cell2Y(jj);
        SimParams.OutFolderName = './Outputs_b_NoPhage/OutputsbY/';
        SimParams.PFU = 0;	

        %Run simulation
        simulateWave(SimParams)
    end
end
toc()


%% I vs. b
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2b = [10^-1, 10^-0.5, 1, 10^0.5, 10];
cell2I = [0.01, 0.1, 1, 10, 100];


tic()
parfor ii = 1:length(cell2b)
    for jj = 1:length(cell2I)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.b2 = SimParams.b2 .* cell2b(ii);
        SimParams.irate2 = SimParams.irate2 .* cell2I(jj);
        SimParams.OutFolderName = './Outputs_b_NoPhage/OutputsIb/';
        SimParams.PFU = 0;	

        %Run simulation
        simulateWave(SimParams)
    end
end
toc()
