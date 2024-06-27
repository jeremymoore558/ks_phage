%% Load Wave Data
clear all
close all

FileList = dir('../Outputs_NoPhageControl/OutputsIY/*.mat');

%Store measurements (PoplationCell 1, Pop. Cell 2, relative Chi,
%relative i)
Measurements = zeros(length(FileList), 6); 

for ii = 1:length(FileList)
    File = ['../Outputs_NoPhageControl/OutputsIY/', FileList(ii).name];
    load(File)

    %% Get wave at last timepoint
    rho_cell_end = rho_cell_store(end, :);
    rho_cell2_end = rho_cell2_store(end, :);

    % set(gca', 'yscale', 'log')

%     %% Get indices of the wave
%     %Find the peaks and start of wave 1
%     peaks = findpeaks(rho_cell_end);
%     peakInd = find(rho_cell_end == peaks(end));
%     WaveStartInd = max(find(rho_cell_end(1:peakInd) <= min(rho_cell_end(1:peakInd))));
% 
%     %Find peaks and start of wave 2
%     peaks2 = findpeaks(rho_cell2_end);
%     peakInd2 = find(rho_cell2_end == peaks2(end));
%     WaveStartInd2 = max(find(rho_cell2_end(1:peakInd2) <= min(rho_cell2_end(1:peakInd))));

%     figure();
%     plot(rho_cell_end)
%     hold on
%     xline(peakInd)
%     xline(WaveStartInd)
% 
%     plot(rho_cell2_end)
%     xline(peakInd2)
%     xline(WaveStartInd2)

    %% Sum cell density in the wave
    cell_population = sum(rho_cell_end(1:end));
    cell2_population = sum(rho_cell2_end(1:end));
    
    %% Get Parameters
    relativeY = SimParams.Y2 ./ SimParams.Y;
    RelativeI = SimParams.irate2 ./ SimParams.irate;
    
    %% Store measurements
    Measurements(ii, :) = [cell_population, cell2_population, relativeY, RelativeI, SimParams.Y, SimParams.irate];
end

%% Plot
relY = Measurements(:, 3);
relI = Measurements(:, 4);
relPopulation = Measurements(:, 2) ./ Measurements(:, 1);


figure()
scatter(relY, relI, 60, log(relPopulation), 'filled')
a = colorbar;
colormap Jet
a.Label.String = 'log N2/N1';
a.Label.Rotation = -90;
xlabel("Relative Y")
ylabel("Relative I")
set(gca,'xscale', 'log', 'yscale', 'log')
saveas(gcf, './I_vs_Y_NoPhage_tot.png')

%% Export Data
headers = ["Cell_population", "Cell2_population", "relativeY", "relativeI", "Y", "irate"];
textHeader = strjoin(headers, ',');
fid = fopen('./I_vs_Y_NoPhage_tot.csv', 'w');
fprintf(fid,'%s\n',textHeader)
fclose(fid)

dlmwrite('./I_vs_Y_NoPhage_tot.csv',Measurements,'-append');