%% Load Wave Data
clear all
close all

FileList = dir('../Outputs_b/OutputsIb/*.mat');

%Store measurements (PoplationCell 1, Pop. Cell 2, relative Chi,
%relative i)
Measurements = zeros(length(FileList), 6); 

for ii = 1:length(FileList)
    File = ['../Outputs_b/OutputsIb/', FileList(ii).name];
    load(File)

    %% Get wave at last timepoint
    rho_cell_end = rho_cell_store(end, :);
    rho_cell2_end = rho_cell2_store(end, :);

    % set(gca', 'yscale', 'log')

    %% Get indices of the wave
    %Find the peaks and start of wave 1
%     peaks = findpeaks(rho_cell_end);
%     if ~isempty(peaks)
%         peakInd = find(rho_cell_end == peaks(end));
%         WaveStartInd = max(find(rho_cell_end(1:peakInd) <= min(rho_cell_end(1:peakInd))));
%         cell_population = sum(rho_cell_end(WaveStartInd:end));
%     else
%         cell_population = sum(rho_cell_end(1:end));% Added to catch cases where all the cells have died
%         ii
%     end
% 
%     %Find peaks and start of wave 
%     peaks2 = findpeaks(rho_cell2_end);
%     if ~isempty(peaks2)
%         peakInd2 = find(rho_cell2_end == peaks2(end));
%         WaveStartInd2 = max(find(rho_cell2_end(1:peakInd2) <= min(rho_cell2_end(1:peakInd))));
%         cell2_population = sum(rho_cell2_end(WaveStartInd:end));
%     else
%         cell2_population = sum(rho_cell2_end(1:end));
%         ii
%     end
%         
% 

    %% Sum cell density in the wave
    cell_population = sum(rho_cell_end(1:end));
    cell2_population = sum(rho_cell2_end(1:end));
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
    
    
    
    %% Get Parameters
    relativeI = SimParams.irate2 ./ SimParams.irate;
    Relativeb = SimParams.b2 ./ SimParams.b;
    
    %% Store measurements
    Measurements(ii, :) = [cell_population, cell2_population, relativeI, Relativeb, SimParams.irate, SimParams.b];
end

%% Plot
relI = Measurements(:, 3);
relb = Measurements(:, 4);
relPopulation = Measurements(:, 2) ./ Measurements(:, 1);


figure()
scatter(relI, relb, 60, log(relPopulation), 'filled')
a = colorbar;
colormap Jet
a.Label.String = 'log N2/N1';
a.Label.Rotation = -90;
xlabel("Relative I")
ylabel("Relative b")
set(gca,'xscale', 'log', 'yscale', 'log')
saveas(gcf, './b_vs_I.png')

%% Export Data
headers = ["Cell_population", "Cell2_population", "relativeI", "relativeb", "irate", "b"];
textHeader = strjoin(headers, ',');
fid = fopen('./b_vs_I.csv', 'w');
fprintf(fid,'%s\n',textHeader)
fclose(fid)

dlmwrite('./b_vs_I.csv',Measurements,'-append');