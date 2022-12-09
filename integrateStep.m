function [drho_cell_dt, drho_cell2_dt, drho_phage_dt, dAdt, dRdt] = integrateStep(rho_cell_1, rho_cell2_1, rho_phage_1, A1, R1, SimParams)
    % initialize
    drho_cell_dt = zeros(size(rho_cell_1));
    drho_cell2_dt = zeros(size(rho_cell_1));
    drho_phage_dt = zeros(size(rho_phage_1));
    dAdt = zeros(size(A1));
    dRdt = zeros(size(R1));
    fS = @(L, Ki, Ka) log((1+L/Ki)./(1+L/Ka));
    
    dx = SimParams.dx;
    dt = SimParams.dt;
    
    %% Handle cell dynamics
    %Calculate perceived signal
    f = fS(A1, SimParams.KiA, SimParams.KaA);
  
    % calculate chemotactic drift
    g = zeros(size(rho_cell_1)); % initialize drift term as temporary variable
    g(:,2:end-1) = SimParams.Chi.*(f(3:end) - f(1:end-2))./(2*dx);
%     g(:,2:end-1) = SimParams.Chi.*(A1(3:end) - A1(1:end-2))./(2*dx);
    
    % chemotactic flux
    g = g .* rho_cell_1;
    g(A1==0) = 0;
    
    % divergence of all chemotactic and diffusive fluxes
    drho_cell_dt(:,2:end-1,:) = SimParams.Mu .* (rho_cell_1(:,1:end-2,:) + rho_cell_1(:,3:end,:) - 2*rho_cell_1(:,2:end-1,:))./(dx.^2) ... %Diffusion term
        - (g(:,3:end,:) - g(:,1:end-2,:))./(2*dx); %Chemotactic term
    
    %Growth
    drho_cell_dt(:,2:end-1,:) = drho_cell_dt(:,2:end-1,:) + SimParams.cR_eff.*rho_cell_1(:,2:end-1,:).*R1(2:end-1)./(R1(2:end-1) + SimParams.KmR);
    
    %Predation
    drho_cell_dt(:,2:end-1,:) = drho_cell_dt(:,2:end-1,:) - SimParams.irate .* rho_cell_1(:,2:end-1,:) .* rho_phage_1(:,2:end-1,:);
   
    %% Handle cell dynamics, second phenotype
    %Calculate perceived signal
    f = fS(A1, SimParams.KiA, SimParams.KaA);
  
    % calculate chemotactic drift
    g = zeros(size(rho_cell2_1)); % initialize drift term as temporary variable
    g(:,2:end-1) = SimParams.Chi2.*(f(3:end) - f(1:end-2))./(2*dx);
%     g(:,2:end-1) = SimParams.Chi.*(A1(3:end) - A1(1:end-2))./(2*dx);
    
    % chemotactic flux
    g = g .* rho_cell2_1;
    g(A1==0) = 0;
    
    % divergence of all chemotactic and diffusive fluxes
    drho_cell2_dt(:,2:end-1,:) = SimParams.Mu .* (rho_cell2_1(:,1:end-2,:) + rho_cell2_1(:,3:end,:) - 2*rho_cell2_1(:,2:end-1,:))./(dx.^2) ... %Diffusion term
        - (g(:,3:end,:) - g(:,1:end-2,:))./(2*dx); %Chemotactic term
    
    %Growth
    drho_cell2_dt(:,2:end-1,:) = drho_cell2_dt(:,2:end-1,:) + SimParams.cR2.*SimParams.Y2.*rho_cell2_1(:,2:end-1,:).*R1(2:end-1)./(R1(2:end-1) + SimParams.KmR);
    
    %Predation
    drho_cell2_dt(:,2:end-1,:) = drho_cell2_dt(:,2:end-1,:) - SimParams.irate2 .* rho_cell2_1(:,2:end-1,:) .* rho_phage_1(:,2:end-1,:);

    
    %% Handle Phage Dynamics   
    %Diffusion
    drho_phage_dt(:,2:end-1,:) = SimParams.MuP .* (rho_phage_1(:,1:end-2,:) + rho_phage_1(:,3:end,:) - 2*rho_phage_1(:,2:end-1,:))./(dx.^2);
    
    %Predation
    drho_phage_dt(:,2:end-1,:) = drho_phage_dt(:,2:end-1,:) + SimParams.irate.*(SimParams.b - 1).*rho_cell_1(:,2:end-1,:) .* rho_phage_1(:,2:end-1,:) + ...
        SimParams.irate2.*(SimParams.b - 1).*rho_cell2_1(:,2:end-1,:) .* rho_phage_1(:,2:end-1,:);
    
    
    %% Handle attractant dynamics
    maximalConsumption = squeeze(sum(rho_cell_1,1)) * SimParams.cA + squeeze(sum(rho_cell2_1,1)) * SimParams.cA2;
    tempA = A1./ (A1+SimParams.KmA);
    tempA(A1==0) = 0; % allows KmA = 0
    GammaA = maximalConsumption(:)' .* tempA; % uM/s, aspartate consumption term
    dAdt(2:end-1) = SimParams.DA * (A1(1:end-2) + A1(3:end) - 2*A1(2:end-1))./(dx.^2) ... %Diffusion
        - GammaA(2:end-1); %Consumption
    
    %% Handle Nutrient Dynamics
%     maximalConsumption = squeeze(sum(rho_cell_1,1)) * SimParams.cR;
    maximalConsumption = squeeze(sum(rho_cell_1,1)) * SimParams.cR + squeeze(sum(rho_cell2_1,1)).*SimParams.cR2;
    tempR = R1./ (R1+SimParams.KmR);
    tempR(R1==0) = 0; % allows KmR = 0
    GammaR = maximalConsumption(:)' .* tempR; % uM/s, aspartate consumption term
    dRdt(2:end-1) = SimParams.DR * (R1(1:end-2) + R1(3:end) - 2*R1(2:end-1))./(dx).^2 ... %Diffusion
        - GammaR(2:end-1); %Consumption


    %% Handle no-flux boundary conditions:
    %For this boundary condition, dc(x=2)/dt|x=0 = dc(x = 0)/dt and same
    %for the end

%     drho_cell_dt
    drho_cell_dt(:, 1, :) = drho_cell_dt(:, 2, :);
    drho_cell_dt(:, end, :) = drho_cell_dt(:, end-1,:);

%     drho_cell2_dt 
    drho_cell2_dt(:, 1, :) = drho_cell2_dt(:, 2, :);
    drho_cell2_dt(:, end, :) = drho_cell2_dt(:, end-1,:);

%     drho_phage_dt
    drho_phage_dt(:, 1, :) = drho_phage_dt(:, 2, :);
    drho_phage_dt(:, end, :) = drho_phage_dt(:, end-1, :);

%     dAdt
    dAdt(1) = dAdt(2);
    dAdt(end) = dAdt(end-1);

%     dRdt
    dRdt(1) = dRdt(2);
    dRdt(end) = dRdt(end-1);


end