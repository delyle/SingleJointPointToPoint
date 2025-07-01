blankSlate

import casadi.*
           
guessFile = '';
saveName = 'tmp';
showFigs = true;
lO = 0.32;
Psim.T = 0.4; % final time
Psim.dt = 5e-3;

linSolver = 'ma57';

% simulation parameters
Psim.M = 43;
Psim.Fmax = 1300;
Psim.Vmax = 1.6;
Psim.lO = lO;
Psim.actOrder = 1;


% parameters for hill-type model, based on Martula 2022 and Otten 1985
Psim.b = [2, -1.3, 0.53]; % force-length NOTE: b(1) CANNOT be less than 2, or the hessian is undefined!
Psim.c = [0.05,5,1]; % passive stretch
Psim.d = [4, 1.8, 30]; % force-velocity
Psim.alp = [94.5, 65]*1e-3; % activation time constants for O(1) model
Psim.s = 500; % smoothing parameter for discontinuity in O(1) model

thetaRange = [60,120]*pi/180;
restingTheta = pi/2;
momentArm = 0.04;
theta2l = @(o) (o-restingTheta)*momentArm + Psim.lO;

Start_l1 = theta2l(thetaRange(1));
D = 2*Psim.lO;
trackCase = theta2l(thetaRange(2))*[1 1];  % could also be a file name, vector, or 'sigmoid', 'halfsigmoid', 'flat'

Start_dl1 = 0;

Start_ldl = [Start_l1;Start_dl1;D-Start_l1;-Start_dl1];


End_act = zeros(2,1);
End_act_tol = 1e-2*ones(2,1); % if set to zero, will not converge as activation can only approach zero once excited
End_dact_tol = 1e-1*ones(2,1);
setEndAct = false; % if true, will set end activation to set amount, with end amount set to within some tolerance if close to zero or 1
                    % if false, will set end activation to be steady within some tolerance.
EndForceDiffTol = 10/Psim.Fmax;
constrainEndTrack = false; % if true, sets the endpoint to exactly match the tracking value there.

intgType = 'intg'; % 'intg' or ExpEuler
W.exc = 1e-6;
W.act = 0;
W.pow = 0;
W.endTrack = 0.01; % a weight adding a penalty for error at end or tracking period
smoothPower = 500;
smoothType = 'tanh';

optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,W,saveName,guessFile,...
    'EndForceDiffTol',EndForceDiffTol,'showFigs',showFigs,'intgType',intgType)
