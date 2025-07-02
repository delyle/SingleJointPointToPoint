% Plot from parameter sweep (currently only for VelActFmaxStiffSweep)
blankSlate
% load data
mainDir = 'TrackO1ParamSweep\VelActFmaxStiffSweep\Data_45to135_flat_T0p4_baselineFineGridStiffDeact\';
fList = dir([mainDir,'/Data_*.mat']);
[A,S] = loadTrackO1Solutions(mainDir);
saveFig = true;
saveFigDir = 'Paper'; % if empty, will save to mainDir

fprintf('Parameter range:\n----------------------\n')
disp(['c1: ',sprintf('\t%.2f',S.c1_range)])
disp(['act: ',sprintf('\t%.0f  ',S.act_range*1000),'  ms'])
disp(['deact: ',sprintf('\t%.0f  ',S.deact_range*1000),'  ms'])
disp(['Vmax: ',sprintf('\t%.2f  ',S.Vmax_range)])
disp(['Fmax: ',sprintf('\t%.0f  ',S.Fmax_range)])

%% Load cost data

objType = 'sse';
recomputeRMSE = false;
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
J = NaN(size(S.c1Mat));
for i = 1:numel(S.c1Mat)
    if ~isempty(A(i).(JField)) % && isfield(A(i),'J_trackErr_sol')
    J(i) = A(i).(JField);
    end
end
if recomputeRMSE
% recompute into RMSE
J = sqrt(A(1).Psim.dt/A(1).Psim.T*J)/A(1).Psim.lO/2; % get the root-mean-squared error for one muscle (it will be the same for the other).
end
%% Want to plot performance for deact vs stiffness
close all; clc
idx.act = 1;
idx.Vmax = 1;
idx.Fmax = 1;
flipDeact = -1; % -1 to flip deact axis, 1 to leave it.

J_plot = permute(J(:,idx.act,idx.Vmax,idx.Fmax,:),[1,5,2,4,3]);

dJdK = diff(J_plot,1,2);
[mindJdK,imin] = min(abs(dJdK),[],1);
deactAtMin = [NaN,flipDeact*S.deact_range(imin)*1000]; 

singleDeactNearMin = mode(deactAtMin);

figure('color','w')

A = load('coolWhite.mat');
cmap = A.coolWhite;
colormap('parula')

contourf(S.c1_range,flipDeact*S.deact_range*1000,J_plot,12)
hold on
plot(S.c1_range,singleDeactNearMin*ones(1,length(S.c1_range)),'w-','LineWidth',1)

fprintf('Maximum dJdK at Deact = %.0f ms is %.3g\n',singleDeactNearMin*flipDeact,max(dJdK(mode(imin),:)))

xlabel('Stiffness (c_1)')
ylabel('Deactivation time constant [ms]')
s = strsplit(mainDir,'Data_');

titleTxt = [];%titleTxt = strrep(s{end},'_',' ');
subtitleTxt = sprintf('Act. $=%.0f$ ms, $v_{\\rm max}=%.2f \\,\\,\\,l_0/$s, $F_{\\rm max}=%.0f$ N',S.act_range(idx.act)*1e3,S.Vmax_range(idx.Vmax),S.Fmax_range(idx.Fmax));
title(titleTxt,subtitleTxt,'interpreter','latex','fontsize',9)
cbh = colorbar;
title(cbh,cbarTxt)

caseName  = sprintf('_Act%.0fms_Vmax%.2f_Fmax%.0f',S.act_range(idx.act)*1e3,S.Vmax_range(idx.Vmax),S.Fmax_range(idx.Fmax));

ax = gca;
ax.YTickLabel = strrep(ax.YTickLabel,'-','');

% plot additional deactivation lines
lit_deact = [40,48,60,65]; % OpenSim default; Kowakowski 2022; Thelan 2003; Murtola 2023
romnum = {'i','ii','iii','iv','v'};
tb = {'bottom','top','bottom','top'};
nudgeRight = 0.0005*[1,1,1,1];
for i = 1:length(lit_deact)
    c = interp1(-singleDeactNearMin*linspace(-0.4,0.5,size(cmap,1))',cmap,lit_deact(i)+singleDeactNearMin);
    plot(S.c1_range([1,end]),-lit_deact(i)*[1 1],'--','color',c,'linewidth',1)
    text(max(S.c1_range)+nudgeRight(i),-lit_deact(i)+0.5,['(',romnum{i},')'],'verticalalignment','middle','horizontalalignment','left')
end

if saveFig
    if isempty(saveFigDir)
        saveFigDir = mainDir;
    end
    saveName = ['Data_',s{end},caseName,'_DeactVsStiff.pdf'];
    if length(saveName) > 90
        saveName = [caseName(2:end),'_DeactVsStiff.pdf'];
    end
    exportgraphics(gcf,[saveFigDir,'/',saveName])
end
