%% Presets
fmt = '-dpng';
rez = '-r300';
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultLineMarkerSize',5);
set(0,'DefaultLineColor',[0 0 0]);
set(0,'DefaultAxesBox','off');
set(0,'DefaultAxesColor','none');
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultTextInterpreter','tex');
set(0,'DefaultAxesTickLabelInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex');
whitebg('white');
close

%% PANEL ARRANGEMENT MIGHT CHANGE
dense_factor = 10;

stdTB = 0.065;
biophysicalChi = 0; 
DA = 800;
OD0 = 6;
KiA = 3.5; % uM
KmA = 0.5;

GlobalParams_all.growthRate = 0;
GlobalParams_all.asp = 100; % uM
GlobalParams_all.stdTB = stdTB;
GlobalParams_all.reSimulate = 0;
GlobalParams_all.biophysicalChi = biophysicalChi; 
GlobalParams_all.DA = DA;
GlobalParams_all.OD0 = OD0;
GlobalParams_all.KiA = KiA;
GlobalParams_all.KmA = KmA;
GlobalParams_all.IC = 'step'; 
GlobalParams_all.x_discretize_factor = 25;
SimParams_all = [];

generalistTB = 0.140;
liquid_specialist_TB = 0.025;
agar_specialist_TB = 0.32;

GlobalParams_all.peakTB = generalistTB;

ks_dir = what('ks');
if isempty(ks_dir) || length(ks_dir)>1
    ks_dir = pwd;
else
    ks_dir = ks_dir(1).path;
end
[~,ind2] = regexp(ks_dir,[filesep,'ks']);
ks_dir = ks_dir(1:ind2);

save_dir = [ks_dir,filesep,'analysisFuns'];

ODumToCellsPercm2 = 8.4e8*(10^-2/10^-6); % OD * um * cells/cm^3/OD * cm/um

TBupper = 0.4;

case_name = ['_biophysChi=' num2str(GlobalParams_all.biophysicalChi) '_DA=' num2str(GlobalParams_all.DA) '_KiA=' num2str(GlobalParams_all.KiA)];

if ~isempty(KmA)
     case_name = [case_name, '_KmA=' num2str(GlobalParams_all.KmA) '_' num2str(x_discretize_factor)];
end


green = [126,186,86]/255;
orange = [241,160,105]/255;
colors = [green;orange];

tmin = 1; % hr
dtplot = [0,0.5,2,6];

liq_color_factors = logspace(log10(1/2),log10(1/max(colors(1,:))),length(dtplot)-1);
ag_color_factors = logspace(log10(1/2),log10(1/max(colors(2,:))),length(dtplot)-1);

liq_colors = repmat(colors(1,:),length(dtplot)-1,1).*repmat(liq_color_factors(:),1,3);
ag_colors = repmat(colors(2,:),length(dtplot)-1,1).*repmat(ag_color_factors(:),1,3);

fsize = [1.7,1.29];
fs = 10;
lw = 1.5;

fsize_slide = [4,3];


%% liquid calculations
environment = 'liquid'; 

GlobalParams_all.environment = environment;

[GlobalParams,SimParams] = initializeSimulationStructures(GlobalParams_all,SimParams_all);
ep = GlobalParams.KiA/GlobalParams.asp;

load(GlobalParams.simName)

% map indices within each block to indices if everything were all in one
% vector
blockInds = cell(1,length(SimParams));
for i = 1:length(SimParams)
    if i==1
        blockInds{i} = 1:length(SimParams(i).time);
    else
        blockInds{i} = blockInds{i-1}(end) + (1:length(SimParams(i).time));
    end
end

if exist(['Fig1D_logZ_liq' case_name '.mat'],'file')==2
    load(['Fig1D_logZ_liq' case_name '.mat']);
    loaded = 1;
else
    loaded = 0;
    logZ_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
end


% leakage over time
TB = [];
PTB_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
NF0_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
t = nan(1,blockInds{end}(end));
c_t = nan(1,blockInds{end}(end));
PF0back_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
chimuBack_t = nan(1,blockInds{end}(end));

for i = 1:length(SimParams)
    F0 = SimParams(i).F0;
    TB = 1./(1+exp(F0));

    t(blockInds{i}) = SimParams(i).time(:);
    
    x_i = SimParams(i).x;
    rho_t = SimResults(i).rho;
    asp_t = SimResults(i).asp;
    dTBdF0 = TB.*(1-TB);
    
    x_i_dense = linspace(min(x_i),max(x_i),length(x_i)*dense_factor);
    dx_i_dense = SimParams(i).dx/dense_factor;
    MFT = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);
    chi = MFT.chi;
    mu = MFT.mu;

    for j = 1:length(blockInds{i})
        tInd = blockInds{i}(j);
        
        asp_j = interp1(x_i(:),asp_t(:,j),x_i_dense(:),'pchip');
        
        rho_j = zeros(size(rho_t,1),length(x_i_dense));
        for k = 1:size(rho_t,1)
            rho_j(k,:) = interp1(x_i,squeeze(rho_t(k,:,j)),x_i_dense,'pchip');
        end
        
        backInd = find(asp_j<=sqrt(GlobalParams.KmA*GlobalParams.KiA),1,'last');
        
        if ~isempty(backInd)
            NF0_t(:,tInd) = dx_i_dense*sum(rho_j(:,backInd:end),2);
            NF0_t(NF0_t(:,tInd)<1e-4,tInd) = 0;
            
            if loaded==0
                if i==1 && j==1
                    logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                elseif i>1&&j==1 
                    if all(isnan(logZ_t(:,blockInds{i-1}(end))))
                        logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                    else
                        logZj = logZ_t(:,blockInds{i-1}(end));
                    end
                elseif j>1
                    if all(isnan(logZ_t(:,tInd-1)))
                        logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                    else
                        logZj = logZ_t(:,tInd-1);
                    end
                end
                
                if all(isnan(logZj))
                    logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                end
                
                logZj2 = findLogZ_relerr(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0,logZj); drawnow
                
            else
                logZj2 = logZ_t(:,tInd);
                
            end
            
            PF0 = NF0_t(:,tInd)/trapz(F0,NF0_t(:,tInd),1);
            PTB = PF0./dTBdF0;
            PTB = PTB/(-trapz(TB,PTB));
            
            meanChi = trapz(F0,PF0.*chi);
            meanMu = trapz(F0,PF0.*mu);
            c_t(tInd) = sum(NF0_t(:,tInd),1)*GlobalParams.kA/GlobalParams.asp; % units are OD*um *uM/OD/s / uM = um/s...
            lambda = c_t(tInd)/(meanChi-meanMu);
            
            z = -8/lambda:(SimParams(i).dx):3/lambda;
            
            [rhoKS,sKS] = KSpredict(F0,NF0_t(:,tInd),GlobalParams,SimParams(i).dx,logZj2,z);
            sKS(sKS<0) = 0;
            z0Ind = find(sKS==0,1,'last')+1;
            
            PF0back = rhoKS(:,z0Ind)/abs(trapz(F0,rhoKS(:,z0Ind),1));
            PF0back_t(:,tInd) = PF0back;
            
            chimuBack_t(tInd) = trapz(F0,PF0back.*(chi(:)-mu(:)),1);
            
            
            
            if loaded == 0
                logZ_t(:,tInd) = logZj2;
            end
            PTB_t(:,tInd) = PTB;
        end
        
    end

end

    
save(['Fig1D_logZ_liq' case_name '.mat'],'logZ_t'); % cut short at i = 11, j = 7

% compute predicted leakage
L_t_pred = -ep.*repmat(c_t.^2./chimuBack_t,length(TB),1).*PF0back_t.*(F0(2)-F0(1)).*sum(NF0_t,1);
    
% compute leakage in simulation
t = t(:);
L_ti = ((NF0_t(:,3:end))-(NF0_t(:,1:end-2)))./repmat((t(3:end)-t(1:end-2))',size(NF0_t,1),1); % convert to reasonable units, or divide by N, total cells/area
L_t = [nan(length(TB),1),L_ti,nan(length(TB),1)];


% other calculations
N_t = c_t/GlobalParams.kA*GlobalParams.asp;

L_t = L_t*ODumToCellsPercm2./repmat((TB.*(1-TB)),1,length(t));
L_t_pred = L_t_pred*ODumToCellsPercm2./repmat((TB.*(1-TB)),1,length(t));

Ltotal = -trapz(TB,L_t,1)';
Ltotal_pred = -trapz(TB,L_t_pred,1)';

inds = find(isnan(Ltotal));
for i = 1:length(inds)
    if inds(i)>2 && inds(i)<(length(Ltotal)-1)
        Ltotal(inds(i)) = 1/2*(Ltotal(inds(i)-2)+Ltotal(inds(i)+2));
        
        L_t(:,inds(i)) = 1/2*(L_t(:,inds(i)-2)+L_t(:,inds(i)+2));
    end
end

[t,ui] = unique(t);
N_t = N_t(ui);
L_t = L_t(:,ui);
L_t_pred = L_t_pred(:,ui);
Ltotal = Ltotal(ui);
Ltotal_pred = Ltotal_pred(ui);
PTB_t = PTB_t(:,ui);
c_t = c_t(ui);
logZ_t = logZ_t(:,ui);


% time points are not evently spaced -- interpolate
t_global = 3600*(tmin+dtplot);
nt = length(t_global);

% interpolate
Ltot_global = interp1(t,Ltotal,t_global);
L_global = zeros(size(L_t,1),nt);
L_pred_global = zeros(size(L_t,1),nt);
NTB_global = zeros(size(L_t,1),nt);
for i = 1:size(L_t,1)
    L_global(i,:) = interp1(t,L_t(i,:),t_global);
    L_pred_global(i,:) = interp1(t,L_t_pred(i,:),t_global);
    
    NTB_global(i,:) = interp1(t,PTB_t(i,:),t_global); 
end

c_t_global = interp1(t,c_t,t_global); 
N_t_global = c_t_global/GlobalParams.kA*GlobalParams.asp;
for j = 1:size(NTB_global,2)
    NTB_global(:,j) = N_t_global(j).*NTB_global(:,j)/abs(trapz(TB,NTB_global(:,j))); % each column should integrate to N_t(i)
end

%% liquid plots

% initial condition
[GlobalParams,~] = initializeSimulationStructures(GlobalParams,[]);
TB0 = 1./(1+exp(GlobalParams.F0));
PF0 = GlobalParams.PF0;
dTBdF0 = TB0.*(1-TB0);
PTB0 = PF0./dTBdF0;
PTB0 = PTB0/(-trapz(TB0,PTB0));

% colors = lines(2);
[~,t1] = min(abs(t-t_global(1)));

% total leakage over time
figure;hold on; hLtot=gca; fLtot=gcf;
plot((t(t1:end)-t(t1))/3600,Ltotal(t1:end)/Ltotal(t1),'Color',colors(1,:),'LineWidth',2)
hp=plot((t(t1:end)-t(t1))/3600,Ltotal_pred(t1:end)/Ltotal_pred(t1),'k--','LineWidth',1.5); hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
for i = 2:length(t_global)
    plot(hLtot,(t_global(i)-t_global(1))/3600,Ltot_global(i)/Ltot_global(1),'o','Color',liq_colors(i-1,:),'MarkerFaceColor',liq_colors(i-1,:),'MarkerSize',4);
end

hLtot.YLim(1)=0;
fLtot.Units = 'inches';
fLtot.Position(3:4) = [2.67,2];
hLtot.XLim(2) = 8;
hLtot.Units = 'inches';
hLtot.FontSize = 12;


% TB distribution over time
figure; hold on
hPTB = gca;
fPTB = gcf;
hPTB.FontSize = fs;
hPTB.LineWidth = lw;
hPTB.XTick = 0:0.1:0.4;
fPTB.Units = 'inches';
fPTB.Position(3:4) = fsize;

plot(TB,PTB0,'k','LineWidth',lw)
for i = 2:length(t_global)
    
    plot(TB,NTB_global(:,i)/abs(trapz(TB,NTB_global(:,i))),'Color',liq_colors(i-1,:),'LineWidth',lw);
    
end
xlim([0,TBupper])
hPTB.XTickLabelRotation=0;

print(fPTB,'Fig1D_PTB_t_liq',fmt,rez) % named by delta t
print(fPTB,'Fig1D_PTB_t_liq','-dsvg')


%% agar calculations
environment = 'agar'; 

GlobalParams_all.environment = environment;

[GlobalParams,SimParams] = initializeSimulationStructures(GlobalParams_all,SimParams_all);
ep = GlobalParams.KiA/GlobalParams.asp;

load(GlobalParams.simName)

% map indices within each block to indices if everything were all in one
% vector
blockInds = cell(1,length(SimParams));
for i = 1:length(SimParams)
    if i==1
        blockInds{i} = 1:length(SimParams(i).time);
    else
        blockInds{i} = blockInds{i-1}(end) + (1:length(SimParams(i).time));
    end
end

if exist(['Fig1D_logZ_ag' case_name '.mat'],'file')==2
    load(['Fig1D_logZ_ag' case_name '.mat']);
    loaded = 1;
else
    loaded = 0;
    logZ_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
end


% leakage over time
TB = [];
PTB_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
NF0_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
t = nan(1,blockInds{end}(end));
c_t = nan(1,blockInds{end}(end));
PF0back_t = nan(length(GlobalParams.F0(:)),blockInds{end}(end));
chimuBack_t = nan(1,blockInds{end}(end));

for i = 1:length(SimParams)
    F0 = SimParams(i).F0;
    TB = 1./(1+exp(F0));

    t(blockInds{i}) = SimParams(i).time(:);
    
    x_i = SimParams(i).x;
    rho_t = SimResults(i).rho;
    asp_t = SimResults(i).asp;
    dTBdF0 = TB.*(1-TB);
    
    x_i_dense = linspace(min(x_i),max(x_i),length(x_i)*dense_factor);
    dx_i_dense = SimParams(i).dx/dense_factor;
    MFT = initializeMFTmodel(F0,GlobalParams.environment,GlobalParams.biophysicalChi);
    chi = MFT.chi;
    mu = MFT.mu;

    for j = 1:length(blockInds{i})
        tInd = blockInds{i}(j);
        
        asp_j = interp1(x_i(:),asp_t(:,j),x_i_dense(:),'pchip');
        
        rho_j = zeros(size(rho_t,1),length(x_i_dense));
        for k = 1:size(rho_t,1)
            rho_j(k,:) = interp1(x_i,squeeze(rho_t(k,:,j)),x_i_dense,'pchip');
        end
        
        backInd = find(asp_j<=sqrt(GlobalParams.KmA*GlobalParams.KiA),1,'last');
        
        if ~isempty(backInd)
            NF0_t(:,tInd) = dx_i_dense*sum(rho_j(:,backInd:end),2);
            NF0_t(NF0_t(:,tInd)<1e-4,tInd) = 0;
            
            if loaded==0
                if i==1 && j==1
                    logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                elseif i>1&&j==1 
                    if all(isnan(logZ_t(:,blockInds{i-1}(end))))
                        logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                    else
                        logZj = logZ_t(:,blockInds{i-1}(end));
                    end
                elseif j>1
                    if all(isnan(logZ_t(:,tInd-1)))
                        logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                    else
                        logZj = logZ_t(:,tInd-1);
                    end
                end
                
                if all(isnan(logZj))
                    logZj = findLogZ(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0); drawnow
                end
                
                logZj2 = findLogZ_relerr(F0,NF0_t(:,tInd),SimParams(i).dx,GlobalParams,0,logZj); drawnow
                
            else
                logZj2 = logZ_t(:,tInd);
                
            end
            
            PF0 = NF0_t(:,tInd)/trapz(F0,NF0_t(:,tInd),1);
            PTB = PF0./dTBdF0;
            PTB = PTB/(-trapz(TB,PTB));
            
            meanChi = trapz(F0,PF0.*chi);
            meanMu = trapz(F0,PF0.*mu);
            c_t(tInd) = sum(NF0_t(:,tInd),1)*GlobalParams.kA/GlobalParams.asp; % units are OD*um *uM/OD/s / uM = um/s...
            lambda = c_t(tInd)/(meanChi-meanMu);
            
            z = -8/lambda:(SimParams(i).dx):3/lambda;
            
            [rhoKS,sKS] = KSpredict(F0,NF0_t(:,tInd),GlobalParams,SimParams(i).dx,logZj2,z);
            sKS(sKS<0) = 0;
            z0Ind = find(sKS==0,1,'last')+1;
            
            PF0back = rhoKS(:,z0Ind)/abs(trapz(F0,rhoKS(:,z0Ind),1));
            PF0back_t(:,tInd) = PF0back;
            
            chimuBack_t(tInd) = trapz(F0,PF0back.*(chi(:)-mu(:)),1);
            
            PF0_theory = sum(rhoKS(:,z0Ind:end),2);
            PF0_theory = PF0_theory/sum(PF0_theory,1);
            chimu_t(tInd) = sum(PF0_theory.*(chi(:)-mu(:)));
            
            if loaded == 0
                logZ_t(:,tInd) = logZj2;
            end
            PTB_t(:,tInd) = PTB;
        end
        
    end

end

    
save(['Fig1D_logZ_ag' case_name '.mat'],'logZ_t');


% compute predicted leakage
L_t_pred = -ep.*repmat(c_t.^2./chimuBack_t,length(TB),1).*PF0back_t.*(F0(2)-F0(1)).*sum(NF0_t,1);
    
% compute leakage in simulation
t = t(:);
L_ti = ((NF0_t(:,3:end))-(NF0_t(:,1:end-2)))./repmat((t(3:end)-t(1:end-2))',size(NF0_t,1),1); % convert to reasonable units, or divide by N, total cells/area
L_t = [nan(length(TB),1),L_ti,nan(length(TB),1)];


% other calculations
N_t = c_t/GlobalParams.kA*GlobalParams.asp;

L_t = L_t*ODumToCellsPercm2./repmat((TB.*(1-TB)),1,length(t));
L_t_pred = L_t_pred*ODumToCellsPercm2./repmat((TB.*(1-TB)),1,length(t));

Ltotal = -trapz(TB,L_t,1)';
Ltotal_pred = -trapz(TB,L_t_pred,1)';

inds = find(isnan(Ltotal));
for i = 1:length(inds)
    if inds(i)>2 && inds(i)<(length(Ltotal)-1)
        Ltotal(inds(i)) = 1/2*(Ltotal(inds(i)-2)+Ltotal(inds(i)+2));
        
        L_t(:,inds(i)) = 1/2*(L_t(:,inds(i)-2)+L_t(:,inds(i)+2));
    end
end

[t,ui] = unique(t);
N_t = N_t(ui);
L_t = L_t(:,ui);
L_t_pred = L_t_pred(:,ui);
Ltotal = Ltotal(ui);
Ltotal_pred = Ltotal_pred(ui);
PTB_t = PTB_t(:,ui);
c_t = c_t(ui);
logZ_t = logZ_t(:,ui);

% time points are not evently spaced
nt = length(t_global);

% interpolate
Ltot_global = interp1(t,Ltotal,t_global);
L_global = zeros(size(L_t,1),nt);
L_pred_global = zeros(size(L_t,1),nt);
NTB_global = zeros(size(L_t,1),nt);
for i = 1:size(L_t,1)
    L_global(i,:) = interp1(t,L_t(i,:),t_global);
    L_pred_global(i,:) = interp1(t,L_t_pred(i,:),t_global);
    
    NTB_global(i,:) = interp1(t,PTB_t(i,:),t_global); 
end

c_t_global = interp1(t,c_t,t_global); 
N_t_global = c_t_global/GlobalParams.kA*GlobalParams.asp;
for j = 1:size(NTB_global,2)
    NTB_global(:,j) = N_t_global(j).*NTB_global(:,j)/abs(trapz(TB,NTB_global(:,j))); % each column should integrate to N_t(i)
end


%% agar plots

figure; hold on
hPTB = gca;
fPTB = gcf;
hPTB.FontSize = fs;
hPTB.LineWidth = lw;
hPTB.XTick = 0:0.1:0.4;
fPTB.Units = 'inches';
fPTB.Position(3:4) = fsize;


% composition
plot(TB,PTB0,'k','LineWidth',lw)
for i = 2:length(t_global)
    
    plot(TB,NTB_global(:,i)/abs(trapz(TB,NTB_global(:,i))),'Color',ag_colors(i-1,:),'LineWidth',lw);

end
xlim([0,TBupper])
hPTB.XTickLabelRotation=0;
print(fPTB,'Fig1D_PTB_t_ag',fmt,rez) 
print(fPTB,'Fig1D_PTB_t_ag','-dsvg')


[~,t1] = min(abs(t-t_global(1)));

% total leakage
plot(hLtot,(t(t1:end)-t(t1))/3600,Ltotal(t1:end)/Ltotal(t1),'LineWidth',2,'Color',colors(2,:))
hp=plot(hLtot,(t(t1:end)-t(t1))/3600,Ltotal_pred(t1:end)/Ltotal_pred(t1),'k--','LineWidth',1.5); hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
for i = 2:length(t_global)
    plot(hLtot,(t_global(i)-t_global(1))/3600,Ltot_global(i)/Ltot_global(1),'o','Color',ag_colors(i-1,:),'MarkerFaceColor',ag_colors(i-1,:),'MarkerSize',4);
end


%% save figures

print(fLtot,'FigS3',fmt,rez)
print(fLtot,'FigS3','-dsvg')
