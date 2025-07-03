% plotsFigure 7 from the manuscript

blankSlate
% load data
mainDir = 'TrackO1ParamSweep/VelActFmaxStiffSweep/Data_45to135_flat_T0p4_warmStart'; %'TrackO1ParamSweep\DeactSweepSmallerMovement\Data_85to95_flat_T0p2_warmStart';%
fList = dir([mainDir,'/Data_*.mat']);
A = load([mainDir,'/',fList(1).name]);
saveFig = true;
saveFigDir = 'Figures'; % if empty, will save to mainDir

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
figure('color','w','position',[680   219   560   659])

tlh = tiledlayout(5,1,'tileSpacing','compact');

for i = 1:5
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


for i = 1:N
    plot(ax(1),B(i).(tField),B(i).x_t(4,:)/B(i).Psim.lO,'-','color',colors(i,:));

    plot(ax(2),B(i).(tField),B(i).x_t(5,:)/B(i).Psim.lO/B(i).Psim.Vmax,'-','color',colors(i,:));

    plot(ax(3),B(i).(tField),B(i).x_t(6,:),'-','color',colors(i,:));
    plot(ax(3),B(i).(tField),-B(i).x_t(3,:),'-','color',colors(i,:));

    plot(ax(4),B(i).(tField),B(i).F_tActive2/B(i).Psim.Fmax,'-','color',colors(i,:));
    plot(ax(4),B(i).(tField),-B(i).F_tActive1/B(i).Psim.Fmax,'-','color',colors(i,:));

    plot(ax(5),B(i).(tField),min([B(i).F_tActive1;B(i).F_tActive2])/B(i).Psim.Fmax,'-','color',colors(i,:));
    

end


[ax(1:4).XTickLabel] = deal('');

ABC = 'abcde';

for i = 1:5
    text(ax(i),0.03,0.85,['(',ABC(i),')'],'units','normalized','FontSize',10)
end

lgh = legend(ax(1),legendTxt,'location','northeast','numcolumns',2);
title(lgh,sprintf('Deact. time (ms)'));
lgh.Box = 'off';
lgh.Position = [0.6520    0.7973    0.2339    0.1206];

ylabel(ax(1),'Strain','Fontsize',FS)
ylabel(ax(2),'Strain rate','Fontsize',FS)
ylabel(ax(3),'Activation','FontSize',FS)
ylabel(ax(4),'Muscle Force (F_{max})','FontSize',FS)
ylabel(ax(5),'Cocontraction (F_{max})','FontSize',FS)


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
    saveName = ['Data_',sweepName,'_',caseName,'_States_NoStiff.pdf'];
    exportgraphics(gcf,[saveFigDir,'/Figure7.pdf'])
end