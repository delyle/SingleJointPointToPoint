blankSlate

%%% vvvvv  USER INPUTS  vvvvv         
guessFile = ''; % can be blank, or a solution file
saveName = 'example'; % name to save the results. If blank, no data will be saved.
showFigs = true; % whether to show the results when the simulation ends
linSolver = 'ma57'; % 'mumps' or 'ma57'

% Physical model parameters
Psim.M = 43; % Mass equivalent in kg: Inertia / r^2
Psim.Fmax = 1300; % Peak tension in N
Psim.Vmax = 1.6; % maximum strain rate in lengths per second
Psim.lO = 0.32; % muscle fibre length in m
Psim.actOrder = 1; % order of the activation model; set to 1 in this paper

% parameters for hill-type model, based on Martula 2022 and Otten 1985
Psim.b = [2, -1.3, 0.53]; % force-length NOTE: b(1) CANNOT be less than 2, or the hessian is undefined!
Psim.c = [0.05,5,1]; % passive stretch
Psim.d = [4, 1.8, 30]; % force-velocity
Psim.alp = [94.5, 65]*1e-3; % activation time constants for O(1) activation model
Psim.s = 500; % smoothing parameter for discontinuity in O(1) activation model

% parameters for tracking task
Psim.T = 0.4; % final simulation time in s
Psim.dt = 5e-3; % Simulation time step duration in s
thetaRange = [45,135]*pi/180; % start and end arm position in radians
restingTheta = pi/2; % position where muscle is at rest
momentArm = 0.04; % muscle moment arm in meters

theta2l = @(o) (o-restingTheta)*momentArm + Psim.lO;
Start_l1 = theta2l(thetaRange(1)); % determine starting muscle length
Start_dl1 = 0; % starting muscle speed

trackCase = theta2l(thetaRange(2))*[1 1];  % What the simulation should track. Could be a file name, vector, or 'sigmoid', 'halfsigmoid', 'flat'

EndForceDiffTol = 10/Psim.Fmax; % tolerance for difference in force between antagonist and agonist at final timestep, relative to Fmax

W.exc = 0; % Weighting to penalise excitations; set to 0 for this study
W.act = 0; % Weighting to penalise activation; set to 0 for this study
W.pow = 0; % Weighting to penalise muscle work; set to 0 for this study
W.endTrack = 0.01; % a weight adding a penalty for error at end of tracking period

%%% ^^^^^ USER INPUTS ^^^^^^

Start_ldl = [Start_l1;Start_dl1;2*Psim.lO-Start_l1;-Start_dl1]; 

optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,W,saveName,guessFile,...
    'EndForceDiffTol',EndForceDiffTol,'showFigs',showFigs)
