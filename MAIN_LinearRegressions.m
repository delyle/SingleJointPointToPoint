% This script extracts the parameter sweep data, and performs linear
% regressions

blankSlate % clear the workspace

% load data
mainDir = 'TrackO1ParamSweep/VelActFmaxStiffSweep/Data_45to135_flat_T0p4_warmstart';%
fList = dir([mainDir,'/Data_*.mat']);
[A,S] = loadTrackO1Solutions(mainDir);
saveFig = true;
saveFigDir = []; % if empty, will save to mainDir

fprintf('Parameter range:\n----------------------\n')
disp(['c1: ',sprintf('\t%.2f',S.c1_range)])
disp(['act: ',sprintf('\t%.0f  ',S.act_range*1000),'  ms'])
disp(['deact: ',sprintf('\t%.0f  ',S.deact_range*1000),'  ms'])
disp(['Vmax: ',sprintf('\t%.2f  ',S.Vmax_range)])
disp(['Fmax: ',sprintf('\t%.0f  ',S.Fmax_range)])

% Extract performance data

[meanCoact,J] = deal(NaN(size(S.c1Mat)));
for i = 1:numel(S.c1Mat)
    J(i) = A(i).J_trackErr_sol;
    act = A(i).x_t([3,6],:);
    coact = min(act);
    meanCoact(i) = mean(coact);
end

% recompute into RMSE
J = sqrt(A(1).Psim.dt/A(1).Psim.T*J)/A(1).Psim.lO/2; % get the root-mean-squared error for one muscle (it will be the same for the other).


% Normalise data
Nfun = @(v) (v-min(v))./range(v); % normalizes a matrix's columns to be on [0 1];

Xraw = [S.deactMat(:),S.c1Mat(:),S.actMat(:),S.VmaxMat(:),S.FmaxMat(:)];

X = Nfun(Xraw); 
y = Nfun(J(:));

%% Test linear model with all terms

% terms matrix
Tlmfull = zeros(16,6); % need extra column of zeros for response variable
Tlmfull(2:6,1:5) = eye(5); % single terms; start at 2 to include intercept
Tlmfull(7:end,1:5) = [
    1 1 0 0 0
    1 0 1 0 0
    1 0 0 1 0
    1 0 0 0 1
    0 1 1 0 0
    0 1 0 1 0
    0 1 0 0 1
    0 0 1 1 0 
    0 0 1 0 1
    0 0 0 1 1]; %interaction terms

Xtable = array2table([X,y],'VariableNames',{'Deact','Stiff','Act','Vmax','Fmax','RMSE'});

lmFull = fitlm(Xtable,Tlmfull)

%% test linear model; only first order terms
% terms matrix
TlmFO = zeros(16,6); % need extra column of zeros for response variable
TlmFO(4:6,3:5) = eye(3); % single terms; start at 4 to avoid deact, stiff
TlmFO(7,1:5) = [
    1 1 0 0 0]; %interaction terms

lmFirstOrder = fitlm(Xtable,TlmFO)


%% Test linear model, but only do Fmax, vmax and activation and no interactive effects
XtableAFV = Xtable(:,3:end)

TlmAFV = zeros(7,4); % need extra column of zeros for response variable
TlmAFV(2:4,1:3) = eye(3); % only activation, fmax and vmax terms are included

lmFirstOrderNonInteract = fitlm(XtableAFV,TlmAFV)

%% fit performance against coactivation

lmCoact = fitlm(meanCoact(:),J(:)/max(J(:)))

min(meanCoact(:))
max(meanCoact(:))

close all
figure('color','w')
plot(meanCoact(:),J(:),'k.')

xlabel('Mean Coactivation')
ylabel('J_{RMSE}')