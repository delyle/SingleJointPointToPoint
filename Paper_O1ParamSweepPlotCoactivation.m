% Plot from parameter sweep (currently only for VelActFmaxStiffSweep)
blankSlate
% load data
mainDir = 'TrackO1ParamSweep/VelActFmaxStiffSweep/Data_45to135_flat_T0p4_warmStart';
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

%% Want to plot activation profile across deact time constant ranges

idx.c1 = 1;
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
%% Plot the activation profile for each case
%close all
figure('color','w','position',[680   374   560   504])
firstPlot = 'velocity'; % force or velocity

tlh = tiledlayout(4,1,'tileSpacing','compact');

for i = 1:4
ax(i) = nexttile(i);
hold(ax(i),'on')
end

[ax.LineWidth] = deal(1);
[ax.Box] = deal('on');

legendTxt = compose('%.0f',A.deact_range*1000);

% determine if time is in a "t" or "t_track" field;
tField = 't';
if isfield(B(1),'t_track')
    tField = 't_track';
end

N = length(B);
colors = linspecer(N,'sequential');

FS = 8;
AgonistFirstPlot = 1;
    switch lower(firstPlot)
        case 'force'
            firstPlotYlabel = 'Force   [F_{max}]';
        case 'jointtorque'
            firstPlotYlabel = 'Joint Torque [rF_{max}]';
            FS = 8;
            AgonistFirstPlot = 2;
        case 'velocity'
            firstPlotYlabel = 'Velocity  [v_{max}]';
        otherwise
            error('firstPlot should be force, jointtorque or velocity')
    end

for i = 1:N
    switch lower(firstPlot)
        case 'force'
            firstPlotY = B(i).F_tActive2/Psim.Fmax;
        case 'jointtorque'
            firstPlotY = (B(i).F_tActive2-B(i).F_tActive1)/Psim.Fmax;
        case 'velocity'
            firstPlotY = -B(i).x_t(5,:)/A.Psim.lO/Psim.Vmax;
    end
    plot(ax(1),B(i).(tField),firstPlotY,'-','color',colors(i,:));

    plot(ax(2),B(i).(tField),B(i).x_t(6,:),'-','color',colors(i,:));

    plot(ax(3),B(i).(tField),B(i).x_t(3,:),'-','color',colors(i,:));

    plot(ax(4),B(i).(tField),min([B(i).x_t(3,:);B(i).x_t(6,:)]),'-','color',colors(i,:));
end


title(ax(AgonistFirstPlot),'Agonist')
title(ax(3),'Antagonist')
title(ax(4),'Coactivation')
ax(3).YLim = ax(2).YLim;
ax(4).YLim = ax(2).YLim;

[ax(1:3).XTickLabel] = deal('');

ABC = 'abcd';

for i = 1:4
    text(ax(i),0.03,0.85,['(',ABC(i),')'],'units','normalized','FontSize',10)
end

lgh = legend(ax(1),legendTxt);
lgh.Layout.Tile = 'east';
title(lgh,sprintf('Deact. time (ms)'));

ylabel(ax(1),firstPlotYlabel,'Fontsize',FS)
ylabel(tlh,'Activation')
xlabel(tlh,'Time [s]')
sweepName = regexprep(fList.name,'Data_|.mat','');
sweepNameCorr = strrep(sweepName,'_',' ');
caseName = regexprep(A.iterNameFun(B(1).Psim),'alp(.*?)ms_','');
caseNameCorr = strrep(caseName,'_',' ');
titleTxt = [caseNameCorr,' ',sweepNameCorr];
%title(tlh,titleTxt,'fontsize',10)

if saveFig
    if isempty(saveFigDir)
        saveFigDir = mainDir;
    end
    saveName = ['Data_',sweepName,'_',caseName,'_Coactivation_',firstPlot,'.pdf'];
    exportgraphics(gcf,[saveFigDir,'/',saveName])
end