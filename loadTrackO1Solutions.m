function [A,S] = loadTrackO1Solutions(dirName,recomputeReg,irange)

% find main data file in listed directory
fList = dir([dirName,'/Data_*.mat']);
S = load([dirName,'/',fList(1).name]);

if nargin < 3 || isempty(irange)
    irange = 1:S.N;
    if nargin < 2
        recomputeReg = false;
    end
end

% get baseline Psim
if isfield(S,'Psim_baseline')
    Psim = S.Psim_baseline;
else
    % load the baseline case
    Bsol = load([dirName,'/sol1_baseline.mat']);
    Psim = Bsol.Psim;
end
wbh = waitbar(0,'Loading data...');
j = 0;
N = length(irange);
for i = irange
    Psim.alp = [S.actMat(i),S.deactMat(i)];
    Psim.Vmax = S.VmaxMat(i);
    Psim.Fmax = S.FmaxMat(i);
    Psim.c(1) = S.c1Mat(i);

    iterName = S.iterNameFun(Psim);
    if i == 1
        A = load(S.adir([iterName,'/sol_',iterName,'.mat']));
    else
        saveName = S.adir([iterName,'/sol_',iterName]);
        try
            A(i) = load([saveName,'.mat']);
        catch
            saveNameReg = [saveName,'_excReg'];
            if recomputeReg
                % recompute with regularisation
                Wtmp = S.W;
                Wtmp.exc = 0.01;

                optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,Wtmp,saveNameReg,'',...
                    'EndForceDiffTol',0/Psim.Fmax,'linSolver',linSolver,'showFigs',false)
            end
            try
                A(i) = load([saveNameReg,'.mat']);
                regUsed(i) = true;
            catch
            end
        end
    end
    j = j+1;
    waitbar(j/N,wbh,'Loading data...')
end
close(wbh)
