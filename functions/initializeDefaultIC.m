function SimParams = initializeDefaultIC(GlobalParams,SimParams)
% specifies initial OD of cells, rho0. also specifies other initial
% conditions
% SimParams should have length 1 if this is being used properly

% side note: for Nat Comm device dimensions, used 5000 cells initially
% OD2Cells = 8.4e8/1e12*600*14;

L = SimParams.L;

% initial cell density
% rho0 = zeros(nPhen, L);
xScale = SimParams.xScale;
dx = SimParams.dx;
x = 0:dx:((L-1)*dx);

IC = GlobalParams.IC;

switch IC
    case 'sigmoid'
        % cells
        rho0 = GlobalParams.OD0*GlobalParams.PF0(:)*(tanh(1-(x/2000).^2)+1)/(tanh(1)+1);
        
        % aspartate
        asp0 = ones(1,L)*GlobalParams.asp;
        
    case 'step'
       
        % like a step, but a very steep sigmoid instead
        rho0 = GlobalParams.OD0*GlobalParams.PF0(:)*(1./(1+exp((x-xScale)/(xScale/30))));

        % aspartate
        asp0 = ones(1,L)*GlobalParams.asp;
        
end

% store
SimParams.rho0 = rho0;
SimParams.asp0 = asp0;


end