blankSlate
% In this script we sweep through various parameters in a simple
% point-to-point joint movement task

mDir = 'TrackO1ParamSweep/VelActFmaxStiffSweep'; % main directory
dataDir = 'Data_45to135_flat_T0p4_warmStart'; % specific directory to save the data, within main directory

% make the directories and define some convinient functions
addpath(mDir)
adir = @(x) [mDir,'/',dataDir,'/',x];
mkdir(adir(''))
printFigs = @(x,n)  exportgraphics(figure(n),sprintf('%s_%i.pdf',x,n), 'ContentType', 'vector');
iterNameFun = @(Psim) strrep(sprintf('c1_%.0fe-3_alp_%.2f_%.2fms_Vmax%.2f_Fmax%.0f',Psim.c(1)*1000,Psim.alp*1e3,Psim.Vmax,Psim.Fmax),'.','p');

%% common parameters.

% define baseline parameters
baselineFmax = 1300; % peak tension in N
baselineVmax = 1.6; % peak strain rate in lengths / s
baselineMomentArm = 0.04; % in m
baselineRestingTheta = pi/2; % position of the arm (from maximum flexion) where muscle length is at resting; in radians
baselineThetaRange = [45,135]*pi/180; % movement excursion in radians; see figure 1 for definition of rotational excursion
baseline_c = [0.05,5,1]; % parameters for stiffness model
baseline_d = [4, 1.4, 30.24]; % parameters for force-length characteristics
baseline_act = 94.5e-3; % activation time constant in s
baseline_deact = 65e-3; % deactivation time constant in s

% Define parameters to test
Vmax_range = [0.75,1,1.5,2]*baselineVmax;
act_range = [0.15,0.25,0.5,0.75,1,1.5,2]*baseline_act;
deact_range = [0.15,0.25,0.5,0.75,1,1.5,2]*baseline_deact;
c1_range = [0, 0.03:0.02:0.09];
Fmax_range = [0.5,0.8,1,1.5]*baselineFmax;

scale_d2 = true; % if true, maximum eccentric gain is scaled by Fmax so that the maximum eccentric force matches the baseline case
useBaselineSol = false; % if true, uses the baseline solution as a guess for the parameter sweep.
recomputeSolutions = false; % if true, will recompute the solution even if one is already located in the folder
showFigs = false; % whether to show figures after every solution is computed. Setting to true increases runtime.

Psim.T = 0.4; % final time, in s
Psim.dt = 5e-3; % time discretisation, in s
W.act = 0;
W.exc = 0;
W.pow = 0;
W.endTrack = 0.01; % a weight adding a penalty for error at end or tracking period


% simulation parameters
Psim.M = 43;
Psim.lO = 0.32;
Psim.actOrder = 1;


% parameters for hill-type model, based on Murtola 2022 and Otten 1985
Psim.b = [2, -1.3, 0.53]; % force-length NOTE: b(1) CANNOT be less than 2, or the hessian is undefined!
Psim.d = baseline_d; % force-velocity characteristics
Psim.s = 500; % smoothing parameter for discontinuity in O(1) activation model

theta2l = @(o) (o-baselineRestingTheta)*baselineMomentArm;

startEnd_l1 = theta2l(baselineThetaRange)+Psim.lO;

D = 2*Psim.lO;
Start_ldl = [startEnd_l1(1);0;D-startEnd_l1(1);0];

trackCase = startEnd_l1(2)*[1,1];

% create the parameter sweep matrices
[deactMat,actMat,VmaxMat,FmaxMat,c1Mat] = ndgrid(deact_range,act_range,Vmax_range,Fmax_range,c1_range);
N = numel(actMat);

Psim.c = baseline_c; % passive parallel elastic properties (baseline)
Psim.Fmax = baselineFmax;
Psim.Vmax = baselineVmax;
Psim.alp = [baseline_act,baseline_deact]; % activation time constants for O(1) model

% check if all the iteration names are unique
iterNames = strings(N,1);
for i = 1:N
    Psim.alp = [actMat(i),deactMat(i)];
    Psim.Vmax = VmaxMat(i);
    Psim.Fmax = FmaxMat(i);
    Psim.c(1) = c1Mat(i);
    iterNames(i) = iterNameFun(Psim);
end

if numel(iterNames) ~= numel(unique(iterNames))
    error('iterNameFun produces duplicate iterNames with these parameters')
end

clearvars iterNames
%% Get initial guess; baseline


Psim.c = baseline_c; % passive stretch (baseline)
Psim.Fmax = baselineFmax;
Psim.Vmax = baselineVmax;
Psim.alp = [baseline_act,baseline_deact]; % activation time constants for O(1) model

Psim_baseline = Psim;

saveNameBaseline1 = adir('sol1_baseline');
linSolver = 'ma57';
close all

hticOverall = tic;
optimizeSimpleAntagonistTrackO1(Psim_baseline,Start_ldl,trackCase,W,saveNameBaseline1,'',...
    'EndForceDiffTol',0/Psim.Fmax,'linSolver',linSolver)

if showFigs
    for i = 1:2
        printFigs(saveNameBaseline1,i)
    end
end
htoc = toc;
%% Determine how many cases need to be recomputed

nSols = 0;
if ~recomputeSolutions
    listSols = dir(adir('c1_*/sol_*.mat'));
    nSols = length(listSols);
    clearvars listSols;
end

%% run sweep

guessFile = ['',saveNameBaseline1(useBaselineSol==1,:)]; % will be empty if useBaselineSol==false
wbh = waitbar(0,'Running optimisations');
ratebaseline = htoc;
tic;
iDiscount = 0;
for i = 1:N
    figh = findall(0,'type','figure');% get figure handles
    close(setdiff(figh,wbh));% close all except waitbar
    Psim.alp = [actMat(i),deactMat(i)];
    Psim.Vmax = VmaxMat(i);
    Psim.Fmax = FmaxMat(i);
    Psim.c(1) = c1Mat(i);
    if scale_d2
        Psim.d(2) = baseline_d(2)*baselineFmax/Psim.Fmax;
    end
    iterName = iterNameFun(Psim);
    rate = (ratebaseline+htoc)/(i+1-iDiscount);
    iLeft = (N-i+1-nSols);
    minLeft = round(rate*iLeft/60);
    wbTxt = ['Current case: ',sprintf('[%.3f,%.0f,%.0f,%.1f,%.0f]',Psim.c(1),Psim.alp*1e3,Psim.Vmax,Psim.Fmax),sprintf(' est. %i min remaining',minLeft)];
    waitbar(i/N,wbh,wbTxt);
    saveDir = adir(iterName);
    mkdir(saveDir)
    saveName = [saveDir,'/sol_',iterName];
    if recomputeSolutions || (~isfile([saveName,'.mat']) && ~isfile([saveName,'_excReg.mat']))
        optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,W,saveName,guessFile,...
            'EndForceDiffTol',10/Psim.Fmax,'linSolver',linSolver,'showFigs',showFigs)
        if showFigs
            for ii = 1:2
                printFigs(saveName,ii)
            end
        end
    else
        nSols = nSols - 1;
        iDiscount = iDiscount+1;
    end
    htoc = toc;
end
close(wbh)

htocOverall = toc(hticOverall) %#ok<NOPTS>

%% Load data
regUsed = false(size(actMat));
recomputeReg = false; % NOTE: if set to true, this will find cases that 
% failed to produce a solution, and attempt to recompute using a
% regularisation term. No regularisation was used in the study
for i = 1:N
    Psim.alp = [actMat(i),deactMat(i)];
    Psim.Vmax = VmaxMat(i);
    Psim.Fmax = FmaxMat(i);
    Psim.c(1) = c1Mat(i);

    iterName = iterNameFun(Psim);
    saveName = adir([iterName,'/sol_',iterName]);
    try
        A(i) = load([saveName,'.mat']);
    catch
        saveNameReg = [saveName,'_excReg'];
        if recomputeReg
            % recompute with regularisation
            Wtmp = W;
            Wtmp.exc = 1e-5;

            sReg = optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,Wtmp,saveNameReg,'',...
                'EndForceDiffTol',0/Psim.Fmax,'linSolver',linSolver,'showFigs',showFigs);
            % use the regularisation to as guess
            sFromReg = optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,W,saveName,saveNameReg,...
                'EndForceDiffTol',0/Psim.Fmax,'linSolver',linSolver,'showFigs',showFigs);
        end
        if sFromReg
            A(i) = load([saveName,'.mat']);
        elseif sReg
            A(i) = load([saveNameReg,'.mat']);
            regUsed(i) = true;
        end
    end
end

%% Get Cost

J_sse = NaN(size(actMat));
TimeToTarget = NaN(size(actMat));
TimeToStop = NaN(size(actMat));
l_tol = 0.0005;
zci = @(v) find(diff(sign(v))); % finds zero crossings of vector; from https://uk.mathworks.com/matlabcentral/answers/267222-easy-way-of-finding-zero-crossing-of-a-function#answer_209072

for i = 1:length(A(:))
    if ~isempty(A(i).J_trackErr_sol)
        J_sse(i) = A(i).J_trackErr_sol;
        j = find(abs(A(i).x_t(1,:)-A(i).l1_track) < l_tol,1);
        if ~isempty(j)
            TimeToTarget(i) = A(i).t_track(j);
        end
        izero =  zci(A(i).x_t(2,:));
        if isempty(izero)
            TimeToStop(i) = A(i).t_track(end);
        else
            TimeToStop(i) = mean(A(i).t_track(izero(1)+[0,1]));
        end
    end
end

J_rmse_l0 = sqrt(Psim.dt/Psim.T*J_sse)/Psim.lO/2; % get the root-mean-squared error for one muscle (it will be the same for the other).

%% Save data
save(adir(dataDir))