% This script produces the Figures S1-S3

% Plot from parameter sweep (currently only for VelActFmaxStiffSweep)
blankSlate
% load data
mainDir = 'TrackO1ParamSweep\VelActFmaxStiffSweep\Data_45to135_flat_T0p4_warmStart\';%'TrackO1ParamSweep/VelActFmaxStiffSweep/Data_45to135_flat_T0p4_ageOnly';%
fList = dir([mainDir,'/Data_*.mat']);
[A,S] = loadTrackO1Solutions(mainDir);
saveFig = false;
saveFigDir = []; % if empty, will save to mainDir

fprintf('Parameter range:\n----------------------\n')
disp(['c1: ',sprintf('\t%.2f',S.c1_range)])
disp(['act: ',sprintf('\t%.0f  ',S.act_range*1000),'  ms'])
disp(['deact: ',sprintf('\t%.0f  ',S.deact_range*1000),'  ms'])
disp(['Vmax: ',sprintf('\t%.2f  ',S.Vmax_range)])
disp(['Fmax: ',sprintf('\t%.0f  ',S.Fmax_range)])


%% load performance data

objType = 'sse';
recomputeRMSE = true;
switch lower(objType)
    case 'work'
        JField = 'J_sol';
        cbarTxt = 'Work (J)';
    case 'sse'
        JField = 'J_trackErr_sol';
        recomputeRMSE = true;
        cbarTxt = 'RMSE (l_0)';
end
% load cost data
[J,meanCoact] = deal(NaN(size(S.c1Mat)));
for i = 1:numel(S.c1Mat)
    if ~isempty(A(i).(JField)) % && isfield(A(i),'J_trackErr_sol')
    J(i) = A(i).(JField);

    act = A(i).x_t([3,6],:);
    coact = min(act);
    meanCoact(i) = mean(coact);
    end
end
J_SSE = J;
if recomputeRMSE
% recompute into RMSE
J = sqrt(A(1).Psim.dt/A(1).Psim.T*J)/A(1).Psim.lO/2; % get the root-mean-squared error for one muscle (it will be the same for the other).
end

%% compute time metrics

TimeToTarget = NaN(size(A));
TimeToStable = NaN(size(A));
HomingInErr = NaN(size(A));
meanCoactBallistic = NaN(size(A));
showfigs = false; % for debugging; if enabled, will show figures for every instance, and user must press a key to continue

for i = 1:length(A)
    t = A(i).t_track;
    x_t = A(i).x_t(1,:)/A(i).Psim.lO;
    v_t = A(i).x_t(2,:)/A(i).Psim.lO/A(i).Psim.Vmax;
    l_t = A(i).l1_track/A(i).Psim.lO;
    prox = x_t-l_t;
    
    iTTT = find(abs(prox) <= 0.005,1);
    if isempty(iTTT)
        TimeToTarget(i) = max(t);
    else
    TimeToTarget(i) = t(iTTT);
    n_h = length(l_t) - iTTT;
    HomingInErr(i) = sqrt(sum(prox((iTTT+1):end).^2)/n_h);
    n_b = iTTT;
    Total_SSE(i) = sum(prox.^2);
    Ballistic_SSE(i) = sum(prox(1:iTTT).^2);
    Ballistic_RMSE(i) = sqrt(Ballistic_SSE(i)/n_b);
    
    end

    StableCriterion = abs(v_t) <= 0.01 & abs(prox) <= 0.01;
    iSC = find(StableCriterion,1);
    
    if isempty(iSC)
        iSC = length(t);
        TimeToStable(i) = max(t);
    else
        TimeToStable(i) = t(iSC);
    end
    act = A(i).x_t([3,6],:);
    coact = min(act);
    meanCoactBallistic(i) = mean(coact(1:iSC));
    
    if showfigs
    close all
    subplot(2,1,1)
    plot(t,x_t)
    hold on
    plot(t,l_t,'b--')
    plot(TimeToTarget(i),x_t(iTTT),'r^');
    plot(TimeToStable(i),x_t(iSC),'rp')
    subplot(2,1,2)
    plot(t,v_t)
    pause
    end
    
end

%% Figure S1: Fraction of JRMSE due to ballistic error 

Jrel = sqrt(Ballistic_SSE(:)./Total_SSE(:));
close all
figure('color','w')
histogram(Jrel,'binwidth',0.00001,'normalization','probability')
%set(gca,'yscale','log','xscale','log')
xlim([0.999,1.0001])
prctile(Jrel,10)
xlabel('$J_{RMSE,rel}$','interpreter','latex')
ylabel(sprintf('Fraction of %i trials',numel(Jrel)),'interpreter','latex')
exportgraphics(gcf,'Figures/FigureS1.pdf')

%% Figure S2: compares Time to Target and Time to Stable point against RMSE

close all

figure('position',[683   313   454   475],'color','w');
tlo = tiledlayout('vertical','tilespacing','compact');

nexttile
plot(J(:),TimeToTarget,'.k')
ylabel('Time to Target (s)')
xlabel('RMSE (l_0)')

nexttile
plot(J(:),TimeToStable,'.k')
ylabel('Time to Stable Point (s)')
xlabel('RMSE (l_0)')
exportgraphics(gcf,'Figures/FigureS2.pdf')


%% Figure S3: Plot Coactivation as a predictor of performance and time to stable point / time to Target
close all
figure('color','w','Position',[680    100   560   660])
tlo = tiledlayout(3,1,'TileSpacing','compact');

xplot = meanCoact(:);

nexttile
plot(xplot(:),TimeToTarget,'.k')
ylabel('Time to Target (s)')
set(gca,'XTickLabel',[])

nexttile
plot(xplot(:),TimeToStable,'.k')
ylabel('Time to Stable Point (s)')
set(gca,'XTickLabel',[])

nexttile
plot(xplot(:),J(:),'.k')
ylabel('RMSE (l_0)')
xlabel('Mean Coactivation')

exportgraphics(gcf,'Figures/FigureS3.pdf')
