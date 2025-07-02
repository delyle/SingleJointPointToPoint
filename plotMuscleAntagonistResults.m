function ax = plotMuscleAntagonistResults(t,x_tMat,uMat,lTargetMat,F_tMat,Psim,revInputs)
% this function takes two muscles from the simple antagonist simulations
% and plots them on one figure in two columns
%
% x_t is a matrix as [x_t1;x_t2]. Same idea for u, lTarget, F_t.

if nargin < 7
    revInputs = false;
end

figure('color','w','position',[544         145        1007         718])

Nrow = 5;
Ncol = 2;
cmap = colororder;
ParentOpts = struct('TileSpacing','Tight');
ChildOpts = struct('TileSpacing','compact');
[ax,rowh,Tlh] = nestedTiledLayout([1 Ncol],[Nrow 1],ChildOpts,ParentOpts);

[ax.LineWidth] = deal(1);
for ii = 1:2

idx = ii;
if revInputs
    idx = -ii + 3; 
end
istart = (idx-1)*5+1;
x_t = x_tMat(istart+(0:4),:);
u = uMat(idx,:);
lTarget = lTargetMat(idx,:);
nF = size(F_tMat,1);
istart = (idx-1)*(nF/2-1)+idx;
F_t = F_tMat(istart:(nF/2*idx),:);


plot(ax(ii),t,x_t(1,:)/Psim.lO,'linewidth',2)
hold(ax(ii),'on');
plot(ax(ii),t,lTarget(2:end)/Psim.lO,'b--')

plot(ax(Ncol+ii),t,x_t(2,:)/Psim.lO,'linewidth',2,'color',cmap(2,:));
hold(ax(Ncol+ii),'on');
plot(ax(Ncol+ii),t([1 end]),-Psim.Vmax*[1 1],'--','color',cmap(2,:))


plot(ax(2*Ncol+ii),t,x_t(5,:),'linewidth',2,'color',cmap(3,:));
hold(ax(2*Ncol+ii),'on');
stairs(ax(2*Ncol+ii),t,u,'linewidth',2,'color',cmap(4,:));

if size(F_t,1) == 2
    % interpret the second dimension as a passive component, and the first
    % as an active component
    plot(ax(3*Ncol+ii),t,F_t(1,:)/Psim.Fmax,'linewidth',2,'color',cmap(5,:))
    hold(ax(3*Ncol+ii),'on');
    plot(ax(3*Ncol+ii),t,F_t(2,:)/Psim.Fmax,'--','linewidth',2,'color',cmap(5,:))
    lgdTxt = {'Active force','Passive force'};
else
    plot(ax(3*Ncol+ii),t,F_t/Psim.Fmax,'linewidth',2,'color',cmap(5,:))
    hold(ax(3*Ncol+ii),'on');
    plot(ax(3*Ncol+ii),t([1 end]),F_load/Psim.Fmax,'k--')
    lgdTxt = {'Contractile force'};
end

fl = Otten_fl(x_t(1,:)/Psim.lO,Psim.b);
fv = Otten_fv_nc(-x_t(2,:)/Psim.lO/Psim.Vmax,Psim.d);
fp = Murtola_fp_nc(x_t(1,:)/Psim.lO,Psim.c);
colororderTmp = cmap([1,2,5,3],:);
plot(ax(4*Ncol+ii),t,[fl;fv;fp;x_t(5,:)],'linewidth',2)
set(ax(4*Ncol+ii),'ColorOrder',colororderTmp)

end

[ax(1:end-2).XTickLabels] = deal([]);

legOrient = 'vertical';
legLoc = 'east';
lh = legend(rowh(1).Children(1),{'muscle length','target length'},'orientation',legOrient);
lh.Layout.Tile = legLoc;
ylabel(rowh(1),'Length [l_0]')

lh(2) = legend(rowh(2).Children(1),{['Rate of Muscle',newline,'length change'],'Vmax'},'orientation',legOrient );
ylabel(rowh(2),'Velocity [l_0/s]')
lh(2).Layout.Tile = legLoc;

lh(3) = legend(rowh(3).Children(1),{'activation','excitation'},'orientation',legOrient );
lh(3).Layout.Tile = legLoc;

lh(4) = legend(rowh(4).Children(1),lgdTxt,'orientation',legOrient );
ylabel(rowh(4),'Force [F_{max}]')
lh(4).Layout.Tile = legLoc;

lh(5) = legend(rowh(5).Children(1),{'force-length','force-velocity','passive','activation'},'orientation',legOrient );
lh(5).Layout.Tile = legLoc;


%[ax(2:2:end).YAxisLocation] = deal('right');
xlabel(Tlh,'Time (s)')