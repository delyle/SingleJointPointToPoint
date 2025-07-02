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

%% want to get three conditions

% "youngest"
idx(1).c1 = 1;
idx(1).act = 1;
idx(1).deact = 1;
idx(1).Vmax = 4;
idx(1).Fmax = 4;

% baseline values
idx(2).c1 = 3;
idx(2).act = 5;
idx(2).deact = 5;
idx(2).Vmax = 2;
idx(2).Fmax = 3;

% "oldest"
idx(3).c1 = 5;
idx(3).act = 7;
idx(3).deact = 7;
idx(3).Vmax = 1;
idx(3).Fmax = 1;

for i = 1:3
    deact_val = A.deact_range(idx(i).deact);
    Psim.c = A.c1_range(idx(i).c1);
    Psim.alp = [A.act_range(idx(i).act),deact_val];
    Psim.Vmax = A.Vmax_range(idx(i).Vmax);
    Psim.Fmax = A.Fmax_range(idx(i).Fmax);
    % load the iteration that matches that case
    iterName = A.iterNameFun(Psim);
    solDir = A.adir(iterName);
    fileName = [solDir,'/sol_',iterName];
    B(i) = load(fileName);
end

%%

close all
figure('position',[506   104   508   677],'color','w')
colors = linspecer(3);
relAntColor = 0.40;
LW = 1;
labelFS = 8;
agAntSpec = 'within'; % or within

for i = 1:3
t = B(i).t_track;
x_t = B(i).x_t;
l = x_t(4,:);
v = x_t(5,:);
a = x_t([3,6],:);
F1 = B(i).F_t1;
F2 = B(i).F_t2;


m = 4; n = 1;
sp(1) = subplot(m,n,1);
hold on
if i == 1; plot(t,B(1).l2_track,'k:'); end
plot(t,l,'color',colors(i,:),'Linewidth',LW)
ylabel('Length (m)')
hold on

sp(2) = subplot(m,n,2);
if i == 1;plot(t([1,end]),[0 0],'k-'); end
hold on
plot(t,v,'color',colors(i,:),'LineWidth',LW)
ylabel('Velocity (m/s)')

% find point of peak velocity
[vpeak,ipeak] = max(-v);
Vmax = -B(i).Psim.Vmax*B(i).Psim.lO;
plot(t(ipeak:end),Vmax*ones(size(ipeak:length(t))),'--','color',colors(i,:))
plot(t(ipeak),-vpeak,'o','markersize',2,'color',colors(i,:),'markerfacecolor',colors(i,:))

sp(3) = subplot(m,n,3);
if i == 1;plot(t([1,end]),[0 0],'k-'); end
hold on
plot(t,a(2,:),'color',colors(i,:),'Linewidth',LW)
plot(t,-a(1,:),'color',colors(i,:)+relAntColor*([1,1,1]-colors(i,:)),'Linewidth',LW)
ylabel('Activation')
sp(4) = subplot(m,n,4);
if i == 1;plot(t([1,end]),[0 0],'k-'); end
hold on
plot(t,F2,'color',colors(i,:),'Linewidth',LW)
plot(t,-F1,'color',colors(i,:)+relAntColor*([1,1,1]-colors(i,:)),'Linewidth',LW)
ylabel('Force (N)')

end
xlabel('Time (s)')

subplot(sp(1))
box on
legend(sp(1).Children(end-3:end-1),{'"oldest"','baseline','"youngest"'})
text(0.025,B(1).l2_track(1),'Target','verticalalignment','bottom','fontsize',labelFS)
ylim( [0.28 0.3601]);

subplot(sp(2))
lh = sp(2).Children(end-3:end-2);

text(sp(2),0.2,-B(1).Psim.Vmax*B(1).Psim.lO,'$v_{\rm max}$','interpreter','latex','verticalalignment','bottom','FontSize',labelFS+2)

sp(2).YLim = [-1.1 0.2];

switch lower(agAntSpec)
    case 'side'
text(sp(3),1.01,0.95,'Agonist','rotation',270,'units','normalized','verticalalignment','bottom','horizontalalignment','left','fontsize',labelFS)
text(sp(3),1.01,0.05,'Antagon.','rotation',270,'units','normalized','verticalalignment','bottom','horizontalalignment','right','fontsize',labelFS)

text(sp(4),1.01,0.95,'Agonist','rotation',270,'units','normalized','verticalalignment','bottom','horizontalalignment','left','fontsize',labelFS)
text(sp(4),1.01,0.05,'Antagon.','rotation',270,'units','normalized','verticalalignment','bottom','horizontalalignment','right','fontsize',labelFS)
    case 'within'
        text(sp(3),0.99,0.97,'Agonist','units','normalized','verticalalignment','top','horizontalalignment','right','fontsize',labelFS)
text(sp(3),0.99,0.03,'Antagonist','units','normalized','verticalalignment','bottom','horizontalalignment','right','fontsize',labelFS)

text(sp(4),0.1875,0.97,'Agonist','units','normalized','verticalalignment','top','horizontalalignment','center','fontsize',labelFS)
text(sp(4),0.1875,0.03,'Antagonist','units','normalized','verticalalignment','bottom','horizontalalignment','center','fontsize',labelFS)
end

sp(4).YLim = [-1500 1200];

txt = 'ABCD';
for i = 1:4
    text(sp(i),0,1.03,['(',txt(i),')'],'units','normalized','FontSize',labelFS+2,'verticalalignment','bottom','HorizontalAlignment','left')
end

exportgraphics(gcf, 'D:\Delyle\SimpleAnatagonistParameterSweeps\Paper\TimeSeries.pdf', 'ContentType', 'vector')

%print(gcf,'D:\Delyle\SimpleAnatagonistParameterSweeps\Paper\TimeSeries.png','-dpng','-r150')

