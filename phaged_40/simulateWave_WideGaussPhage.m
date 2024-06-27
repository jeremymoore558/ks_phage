function [] = simulateWave(SimParams)
    %% Simultaneous simulation of phage and bacteria. Last update: 1/26/22 JM
    %Currently basing the code off of Henry's simulations: https://git.yale.edu/emonetlab/ks
    %Looking especially at the function: simulateWave.m
    close all
    
    %% Make filename 
    OutName = SimParams.OutFolderName;
    mkdir(OutName)
    SaveName = ['GaussPhage_', 'SimI2_', (num2str(SimParams.irate2)), 'Chi2_', (num2str(SimParams.Chi2)),...
        'cA2_',(num2str(SimParams.cA2)), 'cR2_',(num2str(SimParams.cR2)),...
        'Y_', (num2str(SimParams.Y2)), 'b_', (num2str(SimParams.b2))];
    SaveName = strrep(SaveName,'.','_');
    %% Initialization
    %There are 5 equations:
    %drho_N1/dt = f(N, N', N'', ...) is the bacteria concentration
    %drho_N2/dt = ...
    %drho_P/dt = f(P, P'', ...) is the phage concentration
    %dR/dt = f(N, X'') is the nutrient concetnration
    %dA/dt = f(N, A'') is the attractant concentration

%     SimParams = parameters();

    % grid sizes
    dx = SimParams.dx;
    dt = SimParams.dt;
    simSize = SimParams.simulationSize + 2; %Add two points for satisfying boundary conditions

    %Declare simulation vectors:
    rho_cell = zeros(1, simSize); %Cell density
    rho_cell2 = zeros(1, simSize); %Cell density of second phage resistance phenotype
    rho_phage = zeros(1, simSize); %Phage density
    R = zeros(1, simSize); %Nutrient Density
    A = zeros(1, simSize); %Attractant density

    %Set initial conditions
    rho_cell(1:simSize) = 1.1 * 10^10 .* normpdf((1:simSize)/1000.*dx, 0, 1) .* dx;
    rho_cell2 = rho_cell;
    %rho_phage(:) = 1.1 * 10^10 .* normpdf((1:simSize)/1000.*dx, 0, 1) .* dx;
    rho_phage(:) = 1.1 * 10^10 .* normpdf((1:simSize)/1000.*dx, 0, 1000) .* dx;
    A(:) = 2.*dx;
    R(:) = 50.*dx;

    % rho_cell(1:simSize) = 250*10^6 .* normpdf((1:simSize)./100, 0, 1);
    % rho_phage(:) = 250*10^6 .* normpdf((1:simSize)./100, 0, 1); %Inoculated with cells
    % rho_cell(1:200) = 2*10^11;
    % rho_phage(1:200) = 2*10^11;
    % A(:) = 2;
    % R(:) = 50;



    %% Plot initial condition
    fig = figure();
    subplot(4, 1, 1)
    plot(rho_cell, 'r')
    hold on
    plot(rho_cell2, 'b')
    ylabel("cell")

    subplot(4, 1, 2)
    plot(rho_phage)
    ylabel("phage")

    subplot(4, 1, 3)
    plot(A)
    ylabel("A")

    subplot(4, 1, 4)
    plot(R)
    ylabel("R")

    %% Simulation
%     nSteps = 5000;
    nSteps = round(SimParams.TFinal ./ dt);
    
    %Initialize storage vectors
    rho_phage_store = [];
    rho_cell_store = [];
    rho_cell2_store = [];
    t_store = [];
    A_store = [];
    R_store = [];
    store_ind = 1;
    
    for ii = 1:nSteps
        A(A < 0) = 0;
        R(R < 0) = 0;
        rho_phage(rho_phage < 0) = 0;
        rho_cell(rho_cell < 0) = 0;
        rho_cell2(rho_cell2 < 0) = 0;
        t = ii.*SimParams.dt;

        %Begin RK4
        % Intermediate variable 1
        rho_cell_1 = rho_cell;
        rho_cell2_1 = rho_cell2;
        rho_phage_1 = rho_phage;
        A1 = A;
        R1 = R;
        [drho_cell_dt_a, drho_cell2_dt_a, drho_phage_dt_a, dAdt_a, dRdt_a] = integrateStep(rho_cell_1, rho_cell2_1, rho_phage_1, A1, R1, SimParams);

        % Intermediate variable 2
        rho_cell_1 = rho_cell + drho_cell_dt_a*dt/2;
        rho_cell_1(rho_cell_1<0) = 0;
        rho_cell2_1 = rho_cell2 + drho_cell2_dt_a*dt/2;
        rho_cell2_1(rho_cell_1<0) = 0;
        rho_phage_1 = rho_phage + drho_phage_dt_a*dt/2;
        rho_phage_1(rho_phage_1<0) = 0;
        A1 = A + dAdt_a*dt/2;
        A1(A1<0) = 0;
        R1 = R + dRdt_a*dt/2;
        R1(R1<0) = 0;
        [drho_cell_dt_b, drho_cell2_dt_b, drho_phage_dt_b, dAdt_b, dRdt_b] = integrateStep(rho_cell_1, rho_cell2_1, rho_phage_1, A1, R1, SimParams);

        % Intermediate variable 3
        rho_cell_1 = rho_cell + drho_cell_dt_b*dt/2;
        rho_cell_1(rho_cell_1<0) = 0;
        rho_cell2_1 = rho_cell2 + drho_cell2_dt_b*dt/2;
        rho_cell2_1(rho_cell2_1<0) = 0;
        rho_phage_1 = rho_phage + drho_phage_dt_b*dt/2;
        rho_phage_1(rho_phage_1<0) = 0;
        A1 = A + dAdt_b*dt/2;
        A1(A1<0) = 0;
        R1 = R + dRdt_b*dt/2;
        R1(R1<0) = 0;
        [drho_cell_dt_c, drho_cell2_dt_c, drho_phage_dt_c, dAdt_c, dRdt_c] = integrateStep(rho_cell_1, rho_cell2_1,  rho_phage_1, A1, R1, SimParams);

        % Intermediate variable 4
        rho_cell_1 = rho_cell + drho_cell_dt_c*dt/2;
        rho_cell_1(rho_cell_1<0) = 0;
        rho_cell2_1 = rho_cell2 + drho_cell2_dt_c*dt/2;
        rho_cell2_1(rho_cell_1<0) = 0;
        rho_phage_1 = rho_phage + drho_phage_dt_c*dt/2;
        rho_phage_1(rho_phage_1<0) = 0;
        A1 = A + dAdt_c*dt/2;
        A1(A1<0) = 0;
        R1 = R + dRdt_c*dt/2;
        R1(R1<0) = 0;
        [drho_cell_dt_d, drho_cell2_dt_d, drho_phage_dt_d, dAdt_d, dRdt_d] = integrateStep(rho_cell_1, rho_cell2_1, rho_phage_1, A1, R1, SimParams);

        % Combine
        rho_cell = rho_cell + (drho_cell_dt_a + 2*drho_cell_dt_b + ...
                                2*drho_cell_dt_c + drho_cell_dt_d)*dt/6;
        rho_cell2 = rho_cell2 + (drho_cell2_dt_a + 2*drho_cell2_dt_b + ...
                                2*drho_cell2_dt_c + drho_cell2_dt_d)*dt/6;
        rho_phage = rho_phage + (drho_phage_dt_a + 2*drho_phage_dt_b + ...
                                2*drho_phage_dt_c + drho_phage_dt_d)*dt/6;
        A = A + (dAdt_a + 2*dAdt_b + 2*dAdt_c + dAdt_d)*dt/6;
        R = R + (dRdt_a + 2*dRdt_b + 2*dRdt_c + dRdt_d)*dt/6;

        if mod(ii, nSteps/5) == 0
            A(A < 0) = 0;
            R(R < 0) = 0;
            rho_phage(rho_phage < 0) = 0;
            rho_cell(rho_cell < 0) = 0;
            rho_cell2(rho_cell2 < 0) = 0;

            subplot(4, 1, 1)
            hold on
            plot(rho_cell, 'r')
            plot(rho_cell2, 'b')

            subplot(4, 1, 2)
            hold on
            plot(rho_phage)

            subplot(4, 1, 3)
            hold on
            plot(A)

            subplot(4, 1, 4)
            hold on
            plot(R)
            
            %Update storage vectors
            rho_phage_store(store_ind, :) = rho_phage;
            rho_cell_store(store_ind, :) = rho_cell;
            rho_cell2_store(store_ind, :) = rho_cell2;
            t_store(store_ind) = t;
            A_store(store_ind, :) = A;
            R_store(store_ind, :) = R;
            store_ind = store_ind + 1;
        end

    end

    %% Final storage vector update
    %Update storage vectors
    rho_phage_store(store_ind, :) = rho_phage;
    rho_cell_store(store_ind, :) = rho_cell;
    rho_cell2_store(store_ind, :) = rho_cell2;
    t_store(store_ind) = t;
    A_store(store_ind, :) = A;
    R_store(store_ind, :) = R;
    
    %% Plot
    setxlim = [0,simSize];

    subplot(4, 1, 1)
    hold on
    plot(rho_cell, 'r');
    plot(rho_cell2, 'b');
    ylabel("cell")
    xlim(setxlim)
    set(gca,'xticklabel',{[]}, 'YScale', 'linear')
    % set(gca,'xticklabel',{[]}, 'YScale', 'log')

    subplot(4, 1, 2)
    hold on
    plot(rho_phage)
    ylabel("phage")
    xlim(setxlim)
    %ylim([0, 10^10]) %Apr 18, 2022
    set(gca,'xticklabel',{[]})


    subplot(4, 1, 3)
    hold on
    plot(A)
    ylabel("A")
    xlim(setxlim)
    set(gca,'xticklabel',{[]})


    subplot(4, 1, 4)
    hold on
    plot(R)
    ylabel("R")
    xlim(setxlim)
    xt=arrayfun(@num2str,get(gca,'xtick')*dx/1000,'un',0);
    set(gca,'xticklabel',xt)
    xlabel("X Position (mm)")

    set(gcf, 'Position', [256.2,88.2,1101.6,678.4000000000001])
    saveas(fig, [OutName, SaveName, '.png'])
    
    save([OutName, SaveName], 'rho_phage_store', 'rho_cell_store', 'rho_cell2_store',...
        't_store', 'A_store', 'R_store', 'SimParams')

end






