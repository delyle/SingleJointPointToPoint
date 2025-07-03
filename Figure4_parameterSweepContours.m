% This script reproduces Figure 4 in the manuscript

blankSlate % clear the workspace

mDir = 'TrackO1ParamSweep/VelActFmaxStiffSweep/';
dataDir = 'Data_45to135_flat_T0p4_warmStart';
adir = @(x) [mDir,'/',dataDir,'/',x];
printFigs = @(x,n)  exportgraphics(figure(n),sprintf('%s_%i.pdf',x,n), 'ContentType', 'vector');
iterNameFun = @(Psim) strrep(sprintf('c1_%.2f_alp_%.2f_%.2fms_Vmax%.2f_Fmax%.0f',Psim.c(1),Psim.alp*1e3,Psim.Vmax,Psim.Fmax),'.','p');

%% load data
load(adir(dataDir))

%% plot range of surfaces

close all
stiffRange = [1,length(c1_range)];
plotTime = false;

irange = length(Vmax_range)+[0, -3];
iirange = length(Fmax_range)+[0, -3];

J=J_rmse_l0;
figure('color','w','position',[ 209   100   776   840])
rangeJ = [min(J(:)),max(J(:))];
floor1em3 = @(x) floor(x*1000)/1000;
contourLevels = floor1em3(rangeJ(1)):2.5e-3:floor1em3(rangeJ(2)+1e-3);
cmax = contourLevels([1,end]);

flipDeact = -1; % For plotting deactivation. -1 to flip, 1 to keep the same.

RomNum = ["i","ii","iii","iv","v","vi","vii","viii"];
j = 0;
StiffText_i = 1;
for iii = stiffRange
    for i = irange
        for ii = iirange
            j = j+1;
            s(j) = subplot(length(iirange)*length(stiffRange),length(irange),j);
            contourf(act_range*1e3,flipDeact*deact_range*1e3,J(:,:,i,ii,iii),contourLevels)
            clim(cmax);
            if i == 1 
                xlabel('Activation (ms)')
            else
                set(gca,'xticklabels',[])
            end
            if ii == 4
                ylabel('Deactivation (ms)')
                text(-0.25,0.2,sprintf('$\\mathbf{v_{max} =%.1f \\,\\,\\, l_0/s}$',Vmax_range(i)),'interpreter','latex','units','normalized','rotation',90)
            end
            
            regUsedHere = regUsed(:,:,i,ii,iii);
            actUsedHere = actMat(:,:,i,ii,iii)*1e3;
            deactUsedHere = deactMat(:,:,i,ii,iii)*1e3;
            hold on
            plot(actUsedHere(regUsedHere(:)),flipDeact*deactUsedHere(regUsedHere(:)),'.','color',0.75*[1 1 1])
            if i == 4 
                title(sprintf('$\\mathbf{F_{max} = %.0f}$ \\bf{N}',Fmax_range(ii)),'interpreter','latex');
            end

            xlim(act_range([1,end])*1e3)
            ylim(sort(flipDeact*deact_range([1,end])*1e3))
            YTL = s(j).YTickLabel;
            s(j).YTickLabel = strrep(YTL,'-','');


            if j == 1 || j == 5
                AB = 'AB';
                StiffText(StiffText_i) = text(-0.3,1.25,sprintf('\\textbf{(%s)} \\quad Parallel Stiffness $\\mathbf{c_1=%.2f}$',AB(StiffText_i), c1_range(iii)),'units','normalized','Fontweight','bold','HorizontalAlignment','left','Interpreter','latex');
                StiffText_i = StiffText_i+1;
            end
            text(0.045,0.955,sprintf('(%s)',RomNum(j)),'units','normalized','verticalalignment','top','horizontalalignment','left','color','w','fontweight','bold')

        end
    end

    [optJ,iOptJ] = min(J(:,:,:,:,iii));
    [iopt,jopt,kopt,iiopt] = ind2sub(size(actMat),iOptJ);

end

cbh = colorbar(s(2),'northoutside','TickLength',0.02);
s(2).Position(3) = s(1).Position(3);
%cbth = title(cbh,'RMSE (l_0)');
cbh.Position([1,2,3]) = [0.362,0.965,0.3351-0.1];
cbh.AxisLocation = 'in';
cbth = title(cbh,'Root-mean-squared position error (l_0)');
cbth.Position = [68.4141   10 0];
cbh.Label.FontSize = 9;

% inset showing young --> old
axIn = axes(gcf,'Position',[0.4211    0.495    0.1200    0.0600],'box','on','xtick','','ytick','');
xlim(axIn,[0 1]);
ylim(axIn,[0 1]);

text(axIn,0.05,1,'Younger','Fontsize',8,'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold')
text(axIn,0.95,0,'Older','Fontsize',8,'HorizontalAlignment','right','VerticalAlignment','bottom','FontWeight','bold')
annotation('arrow',[0.471,0.502],[0.5236,0.5]+0.012,'headsize',8)

newXright = s(2).Position(1)-0.075;

% adjust positions
for i = 2:2:8
    s(i).Position(1) = newXright;
    set(s(i),'YAxisLocation','right')
end

yadjust = -0.035;

for i = 3:8
    s(i).Position(2) = s(i).Position(2) - yadjust*(1+floor((i-3)/2)-2.5*(i>4));
end

%% save the figure

exportgraphics(figure(1),'Figures/Figure4.pdf', 'ContentType', 'vector');
