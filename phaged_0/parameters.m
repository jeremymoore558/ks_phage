%% Parameter file
function params = parameters()
    %Output folder name
    params.OutFolderName = './Outputs/';

    %Time-step
%     params.dt = 0.08; %s^-1
    params.dt = 0.08; %s^-1  
%     params.TFinal = 30000; %s
%     params.TFinal = 5000; %s
%     params.TFinal = 35000; %s
%    params.TFinal = 50000; %s
    params.TFinal = 70000; %s	
    %Initial phage concentration
    params.PFU = 2.5; %PFU/dx
	%params.PFU = 2.5 * 10^9; %PFU/dx for wide phage control
    

    %Mesh parameters
    params.dx = 20; %um
%     params.dx = 1000; %um. Here, each unit is a mm
    params.simulationSize = 1500; %Size of vectors, with meshstepsize spacing between each entry
    
    %Cell Phenotype parameters
    params.Chi = 315; %um^2/s Chemtactic coefficient
    params.Chi2 = 315; %um^2/s Chemtactic coefficient of cell 2
    params.Mu = 50; %um^2/s Diffusion coefficient
    
    %Phage parameters
    params.MuP = 0; %Diffusion coefficient
%     params.irate = 10^(-7); %Infections per cel/ml/pfu/ml/s
    params.irate = 10^(-14); %Infections per cel/ml/pfu/ml/s
    
    %Phage parameters for second cell phenotype
    params.irate2 = 10^(-14); %Infections per cel/ml/pfu/ml/s

    params.b = 30; %Burst size
%     params.b = 20; %Burst size
    params.b2 = 30;

    %Attractant parameters
    params.KiA = 18; %Receptor-binding
    params.KaA = 3000; %Receptor-binding
    
    params.KmA = 1; %Consumption K1/2
    params.DA = 800; %Diffusion coefficient
    params.cA = 4.17*10^(-16); %Consumption rate
%     params.cA = 4.17*10^(-11); %Consumption rate

    params.cA2 = params.cA; %Attractant consumption of phenotype 2

    %Nutrient parameters
    params.KmR = 50; %uM Michaelis constant for glycerol
    params.DR = 800;
    params.cR = 2.25.*10^(-14); %Nutrient consumption rate uM
%     params.cR = 2.25.*10^(-9); %Nutrient consumption rate
    params.Y = 1.229.*10^10; %Cell yield for mmol consumed
%     params.Y = 1.229.*10^9; %Cell yield for mmol consumed
    params.cR_eff = params.cR.*params.Y;    
    
    params.cR2 = params.cR; %Phenotype 2
    params.Y2 = params.Y; %Phenotype 2

    
end
