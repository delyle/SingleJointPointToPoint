% Plot from parameter sweep (currently only for VelActFmaxStiffSweep)
blankSlate
% load data
mainDir = 'TrackO1ParamSweep\VelActFmaxStiffSweep\Data_45to135_flat_T0p4_warmStart\';
fList = dir([mainDir,'/Data_*.mat']);
A = load([mainDir,'/',fList(1).name]);
saveFig = true;
saveFigDir = 'Paper'; % if empty, will save to mainDir

fprintf('Parameter range:\n----------------------\n')
disp(['c1: ',sprintf('\t%.2f',A.c1_range)])
disp(['act: ',sprintf('\t%.0f  ',A.act_range*1000),'  ms'])
disp(['deact: ',sprintf('\t%.0f  ',A.deact_range*1000),'  ms'])
disp(['Vmax: ',sprintf('\t%.2f  ',A.Vmax_range)])
disp(['Fmax: ',sprintf('\t%.0f  ',A.Fmax_range)])

%% Want to plot active antagonist force, passive agonist force
%  and their difference across deact time constant ranges

idx.c1 = 5;
idx.act = 5;
idx.Vmax = 2;
idx.Fmax = 3;

for i = 1:length(A.deact_range)
    deact_val = A.deact_range(i);
    Psim.c = A.c1_range(idx.c1);
    Psim.alp = [A.act_range(idx.act),deact_val];
    Psim.Vmax = A.Vmax_range(idx.Vmax);
    Psim.Fmax = A.Fmax_range(idx.Fmax);
    % load the iteration that matches that case
    iterName = A.iterNameFun(Psim);
    solDir = A.adir(iterName);
    fileName = [solDir,'/sol_',iterName];
    B(i) = load(fileName);
end
%% Plot the force profile for each case
close all
figure('color','w')
meanForceTable = false;

tlh = tiledlayout(4,1,'tileSpacing','compact');

% determine if time is in a "t" or "t_track" field;
tField = 't';
lTargField = 'lTarget1';
lTargStart = 2;
if isfield(B(1),'t_track')
    tField = 't_track';
    lTargField = 'l1_track';
    lTargStart = 1;
end

tol = 0.05;

for i = 1:4
    ax(i) = nexttile(i);
    hold(ax(i),'on')
end
[ax.LineWidth] = deal(1);
[ax.Box] = deal('on');

legendTxt = compose('%.0f',A.deact_range*1000);



N = length(B);
colors = linspecer(N,'sequential');
meanTorque = NaN(1,N);
for i = 1:N
    
    iFinal = find(abs(B(i).x_t(1,:) - B(i).(lTargField)(lTargStart:end))/B(i).Psim.lO < tol,1);
    irange = 1:iFinal;
    t_plot = B(i).(tField)(irange);

    plot(ax(1),t_plot,B(i).F_tPassive2(irange),'-','color',colors(i,:));


    plot(ax(2),t_plot,B(i).F_tActive1(irange),'-','color',colors(i,:));

    netTorque = A.baselineMomentArm*(B(i).F_tPassive2(irange) - B(i).F_tActive1(irange));
    plot(ax(3),t_plot,netTorque,'-','color',colors(i,:));
    meanTorque(i) = 1/t_plot(end)*trapz(t_plot,netTorque(irange));

    plot(ax(4),t_plot,min(B(i).x_t(3,irange),B(i).x_t(6,irange)))
end
title(ax(1),'Agonist Passive Force')
title(ax(2),'Antagonist Active Force')
title(ax(3),'Net Torque about joint')
title(ax(4),'Coactivation')

ax(1).YLim = ax(2).YLim;
ax(1).XLim = ax(2).XLim;

lgh = legend(ax(1),legendTxt,'location','eastoutside');
lgh.Layout.Tile = 'east';
title(lgh,sprintf('Deact. rate (ms)'));
lgh.Position(2) = 0.86;

ylfs = 12;
ylabel(ax(1),'Force [N]','FontSize',ylfs)
ylh2 = ylabel(ax(2),'Force [N]','FontSize',ylfs);
ylh3 = ylabel(ax(3),'Torque [Nm]','FontSize',ylfs);
ylh3.Position(1) = ylh2.Position(1);
xlabel(tlh,'Time [s]')
sweepName = regexprep(fList.name,'Data_|.mat','');
sweepNameCorr = strrep(sweepName,'_',' ');
caseName = regexprep(A.iterNameFun(B(1).Psim),'alp(.*?)ms_','');
caseNameCorr = strrep(caseName,'_',' ');
titleTxt = [caseNameCorr,' ',sweepNameCorr];
%title(tlh,titleTxt,'fontsize',10)

ABC = 'ABC';
for i = 1:3
    text(ax(i),0,1.06,['(',ABC(i),')'],'units','normalized','FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','left')
end

if meanForceTable
    % plot mean force in legend next to last axis
    meanForceTxt = compose('%.1f',meanTorque);
    lgh_mf = legend(ax(3),meanForceTxt,'Location','east');
    title(lgh_mf,'Mean Force')
    lgh_mf.Position(1) = 0.86;
    lgh_mf.FontSize = 6;
end

if saveFig
    if isempty(saveFigDir)
        saveFigDir = mainDir;
    end
    saveName = ['Data_',sweepName,'_',caseName,'_PassiveForce.pdf'];
    exportgraphics(gcf,[saveFigDir,'/',saveName])
end