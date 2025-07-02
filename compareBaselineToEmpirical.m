% Plot single simulated elbow pointing task against empirical data
% Generates Figure 2 in the paper

blankSlate

A = load('60to120Comp.mat'); % Note: this can 

fh = figure('color','w');
fh.Position(2) = fh.Position(2)-200;
lw = 1;
t= A.t_track*1000; % put in ms
l = A.x_t(1,:);
t_adjust = [-80,t];

momentArm = 0.04;
restingTheta = pi/2;

theta2l = @(o) (o-restingTheta)*momentArm + A.Psim.lO;
l2theta = @(l) (l-A.Psim.lO)/(momentArm)+restingTheta; 

thetaDeg = l2theta(l)*180/pi-60;

load('EmpiricalData.mat','Data');



plot(t_adjust,[0,A.u_inp(2,:)],'--','color',[1 0 0],'linewidth',lw)
hold on
lh(1) = plot(t_adjust,[0,-A.u_inp(1,:)],'--','color',[1 0.5 0],'linewidth',lw);

%plot(t,A.x_t(6,:),'k:')
%plot(t,-A.x_t(3,:),'k:')

ylim([-1.5 1.5])
plot(t_adjust([1,end]),[0 0],'-','color',0*[1 1 1])
ylabel('Excitation / Maximum')
lh(:,2) = plot(Data.ago_act_MVC(:,1),Data.ago_act_MVC(:,2)/max(Data.ago_act_MVC(:,2)),'-','color',[0.5 0 0],'linewidth',lw);
plot(Data.ant_act_MVC(:,1),-Data.ant_act_MVC(:,2)/min(Data.ant_act_MVC(:,2)),'-','color',[0.5 0.25 0],'linewidth',lw)
xlim([t_adjust(1),Data.ant_act_MVC(end,1)])

text(-20,0.5,['Agonist',newline,'Excitation'],'HorizontalAlignment','right','verticalalignment','bottom')
text(75,-0.5,['Antagonist',newline,'Excitation'],'HorizontalAlignment','right','verticalalignment','top')

yyaxis(gca,'right');
plot(t,thetaDeg,'--','color',[1 0.25 0])
plot(Data.deg(:,1),Data.deg(:,2),'-','color',[0.5 0.125 0],'linewidth',lw)

ylim([-80, 80])
set(gca,'YColor','k')
ylabel('Relative Displacement (^\circ)')
xlabel('Time (ms)')

text(250,65,'Elbow extension','HorizontalAlignment','right','verticalalignment','bottom')
ytick = get(gca,'ytick');
n = numel(ytick);
ytick(1:floor(n/2)) = [];
set(gca,'ytick',ytick,'linewidth',lw)

legend(lh,'Simulation','Empirical','location','northoutside','NumColumns',3)



%exportgraphics(gcf,'baselineForCompToEmp.emf')
exportgraphics(gcf,'Figures/Figure2.pdf')

