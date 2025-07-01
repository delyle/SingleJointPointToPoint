function success = optimizeSimpleAntagonistTrackO1(Psim,Start_ldl,trackCase,W,varargin)
% Finds the excitations minimizing the sum-squared difference between
% muscle lengths and target lengths over a set time, starting from rest.
% Eventual functionality will allow starting from an arbitrary initial state
%
% trackCase can be simulation data from optimizeSimpleAntagonist, a vector
% of muscle lengths to be interpolated, or the strings 'sigmoid',
% 'halfsigmoid', or 'flat'


p = inputParser;

addRequired(p,'Psim',@isstruct);
addRequired(p,'Start_ldl',@(x) isnumeric(x) && isvector(x));
addRequired(p,'trackCase',@(x) ischar(x) || isstring(x) || isnumeric(x))
addRequired(p,'W',@isstruct)
addOptional(p,'saveName','',@(x) isstring(x) || ischar(x));
addOptional(p,'guessFile','',@(x) isstring(x) || ischar(x) || isempty(x))
addParameter(p,'End_act_tol',1e-2*ones(6,1),@(x) isnumeric(x) && length(x) == 6 && all(x > 0))
addParameter(p,'setEndAct',false,@isscalar) % if a two-element vector, will set end activation to set amount within some tolerance if close to zero or 1
% if false, will set end activation to be steady within some tolerance.
addParameter(p,'End_dact_tol',1e-1*ones(3,1),@(x) isnumeric(x) && length(x) == 3);
addParameter(p,'EndForceDiffTol',1e-2,@(x) isscalar(x) && x >= 0)
addParameter(p,'linSolver','mumps',@ischar);% options are: 'mumps', or (with COIN_HSL) 'MA27','MA57','MA77','MA86', or 'MA97'
addParameter(p,'intgType','intg',@(x) any(strcmp(x,{'intg','ExpEuler'})));
addParameter(p,'smoothPower',100,@(x) isnumeric(x) && isscalar(x))
addParameter(p,'smoothType','tanh',@ischar);
addParameter(p,'showFigs',true,@islogical)
addParameter(p,'constrainEndTrack',false)
addParameter(p,'setInitialAct',false,@(x) isnumeric(x) && length(x) == 2)
% will override the steady force rule at the beginning and use the
% input value [ag, ant]

parse(p,Psim,Start_ldl,trackCase,W,varargin{:});

user_input = p.Results;

saveName = p.Results.saveName;
guessFile = p.Results.guessFile;
End_act_tol = p.Results.End_act_tol; % if set to zero, will not converge as activation can only approach zero once excited
setEndAct = p.Results.setEndAct;
End_dact_tol = p.Results.End_dact_tol;
EndForceDiffTol = p.Results.EndForceDiffTol;
linSolver = p.Results.linSolver;
intgType = p.Results.intgType;
smoothPower = p.Results.smoothPower;
smoothType = p.Results.smoothType;
showFigs = p.Results.showFigs;
setInitialAct = p.Results.setInitialAct;

import casadi.*


if numel(Start_ldl) == 2
Start_l1 = Start_ldl(1);
Start_dl1 = Start_ldl(2);

Start_ldl = [Start_l1;Start_dl1;Psim.lO*2-Start_l1;-Start_dl1];
end


% if false, will set end activation to be steady within some tolerance.
constrainEndTrack = false; % if true, sets the endpoint to exactly match the tracking value there.


if (ischar(trackCase) || isstring(trackCase)) && isfile(trackCase)

        A = load(trackCase);
        X_track = A.x_t;
        l1_track = X_track(1,:);
        l2_track = X_track(4,:);
        t_track = A.t;
        N = length(t_track)-1;
        T = t_track(end);
        dt = diff(t_track(1:2));
else
        T = Psim.T;
        dt = Psim.dt;
        t_track = 0:dt:T; % initialize time vector (for states)
        N = length(t_track)-1; % initialize number of timesteps (these are number of control iterations)
        s = 9;
        sigmoid = @(x,s) tanh(s*x)/2+1/2;
        if ischar(trackCase)
        switch trackCase
            case 'sigmoid'
                l1_track = (sigmoid(t_track-T/2,s)*0.2+0.9)*Psim.lO;
            case 'halfsigmoid'
                l1_track = (sigmoid(t_track,s)*0.2+0.9)*Psim.lO;
            case 'flat'
                l1_track = 1.1*ones(size(t_track))*Psim.lO;
            otherwise
                error('trackCase not recognized')
        end
        else
      % interpolate a numeric vector over the time span
        trackCaseTime = ((1:length(trackCase))-1)/(length(trackCase)-1)*T;
        l1_track = interp1(trackCaseTime,trackCase,t_track);              
        end
        D = 2*Psim.lO;
        l2_track = D-l1_track;
end


% Create Guess
if isempty(guessFile)
    X0 = [Start_ldl(1:2);0;Start_ldl(3:4);0]; % initial state 2*(muscle length, velocity, a12, activation)
    XF = [l1_track(end);0;0;l2_track(end);0;0];
    % interpolate between start and end positions
    X_guess = interp1([1,N+1],[X0,XF]',1:N+1)';
    U_guess = zeros(2,N+1);
else
    A = load(guessFile);
    if isfield(A,'T_sol')
        T_guess = A.T_sol;
    else
        T_guess = A.T;
    end
    N_inp = A.N;
    t_guess = linspace(0,T_guess,N+1);
    if isfield(A,'t_track')
        A.t = A.t_track;
    end
    X_guess = interp1(A.t,A.x_t',t_guess)';
    U_guess = interp1(A.t,A.u_inp',t_guess,'linear','extrap')';
end
Psim.dt = dt;

% generate casadi simulation and force function
switch intgType
    case 'intg'
        [forward_dt,F_fun,F_funActive] = SimpleAntagonist(Psim);
    case 'ExpEuler'
        [forward_dt,F_fun] = SimpleAntagonistVarTimeExpEuler(Psim);
end

if isnumeric(setInitialAct)
Start_act = setInitialAct(:);
else

% determine initial activation based on initial force

F0(1) = full(F_fun(Start_ldl(1:2),0));
F0(2) = full(F_fun(Start_ldl(3:4),0));

% determine which is bigger
[~,i] = max(F0);

% find the activation needed such that the other muscle provides equal
% force

% create a casadi rootfinder object
a_root = MX.sym('a_root',1);
i_start = [3,4;1,2];
F0vs_a = Function('F0vs_a',{a_root},{F_fun(Start_ldl(i_start(i,:)),a_root)-F0(i)});
F0root = rootfinder('F0root','newton',F0vs_a);
% find the required activation
a_root = full(F0root(0.5));
disp(['Left-right difference in initial muscle force is ',num2str(abs(full(F0vs_a(a_root))))]) % check, should be zero

Start_act = a_root*[1;1];
Start_act(i) = 0;
end
%% Multiple Shooting optimal control problem

opti = casadi.Opti();

% desion variable for states
X = opti.variable(6,N+1);

% name states for readability
l1 = X(1,:);
dl1 = X(2,:);
a1 = X(3,:);

l2 = X(4,:);
dl2 = X(5,:);
a2 = X(6,:);
act = [a1;a2];

% Get forces (useful)

F1 = F_fun([l1;dl1],a1);
F2 = F_fun([l2;dl2],a2);

% Decision variables for control vector (excitations)
U = opti.variable(2,N+1);

% Gap closing shooting constraints

switch intgType
    case 'intg'
        for k = 1:N
            opti.subject_to(X(:,k+1) == forward_dt(X(:,k),U(:,k+1)));
        end
    case 'ExpEuler'
        for k = 1:N
            opti.subject_to(X(:,k+1) == forward_dt(X(:,k),U(:,k+1),dt));
        end
end


% endpoint constraints

% start and end at rest (kinematically)
opti.subject_to(X([1,2,4,5],1) == Start_ldl);
opti.subject_to(X([2,5],N+1) == [0;0]); % don't impose end position yet
if constrainEndTrack
    opti.subject_to(X([1,4],N+1) == [l1_track(N+1);l2_track(N+1)]);
end

opti.subject_to(U(:,1) == 0); %#ok<UNRCH>
opti.subject_to(act(:,1) == Start_act); 

if isnumeric(setEndAct) 
    % end at set activation levels
    % first excitation doesn't matter (only for other case)
    
    iZero = setEndAct == 0;
    if any(iZero)
        opti.subject_to(act(find(iZero),N+1) < End_act_tol(iZero)); %#ok<FNDSB> %CasADi doesn't handle logical indexing.
    end
    iOne = setEndAct == 1;
    if any(iOne)
        opti.subject_to(act(find(iOne),N+1) > 1 - End_act_tol(iOne)); %#ok<FNDSB>
    end
    iMid = ~iZero & ~iOne;
    if any(iMid)
        opti.subject_to(act(find(iMid),N+1) == setEndAct(iMid));  %#ok<FNDSB>
    end
else
    % end activation is steady (no rate of change) and determined by
    % optimizer
    
    opti.subject_to(-End_dact_tol < o1ActDyn(act(1,end),U(1,end),Psim.alp) < End_dact_tol)
    opti.subject_to(-End_dact_tol < o1ActDyn(act(2,end),U(2,end),Psim.alp) < End_dact_tol)
    
    if EndForceDiffTol > 0
        opti.subject_to(-EndForceDiffTol < (F1(end) - F2(end))/Psim.Fmax < EndForceDiffTol);
    elseif EndForceDiffTol == 0
        opti.subject_to((F1(end) - F2(end))/Psim.Fmax == 0);
    else
        error('EndForceDiffTol must be > 0')
    end
end

% path constraints
opti.subject_to( 0 <= U(:) <= 1);
opti.subject_to([l1(:);l2(:)] > 0 )
opti.subject_to(0 <= act(:) <= 1)

F1_active = F_funActive([l1;dl1],a1);
F2_active = F_funActive([l2;dl2],a2);

absSmooth = absSmoothFun(smoothType);
Power1 = dl1.*F1_active;
Power2 = dl2.*F2_active;
absPower1 = absSmooth(Power1,smoothPower);
absPower2 = absSmooth(Power2,smoothPower);

% Objective: regularize controls
J_reg = W.act*sumsqr([a1(:);a2(:)]) +...
    W.exc*sumsqr(U(:)) +...
    W.pow*(sum(absPower1(:)+absPower2(:)));
J_trackErr = sumsqr(l1-l1_track)+sumsqr(l2-l2_track);
J_trackErrEnd = W.endTrack*((l1(end)-l1_track(end))^2+(l2(end)-l2_track(end))^2);
J = J_reg+J_trackErr+J_trackErrEnd;
opti.minimize(J);

% initialize guess
opti.set_initial(X,X_guess);
opti.set_initial(U,U_guess);

% solve optimization problem
p_opts = struct();
s_opts = struct('max_iter',1000,'linear_solver',linSolver);
opti.solver('ipopt',p_opts, s_opts);
%opti.solver('sqpmethod',struct('convexify_strategy','eigen-clip','qpsol','qrqp'));

success = false;
try
    sol = opti.solve();
    % Get Result
    x_t = sol.value(X);
    u_inp = sol.value(U);
    J_sol = sol.value(J);
    J_reg_sol = sol.value(J_reg);
    J_trackErr_sol = sol.value(J_trackErr);
    J_trackErrEnd_sol = sol.value(J_trackErrEnd);
    success = true;
catch
    warning('Optimization failed, going into debug mode')
    x_t = opti.debug.value(X);
    u_inp = opti.debug.value(U);
    J_sol = opti.debug.value(J);
end

%%
F_t1 = full(F_fun(x_t(1:2,:),x_t(3,:))); % get net force for muscle 1
F_t2 = full(F_fun(x_t(4:5,:),x_t(6,:))); % get net force for muscle 2

F_tActive1 = full(F_funActive(x_t(1:2,:),x_t(3,:))); % get active force for muscle 1
F_tActive2 = full(F_funActive(x_t(4:5,:),x_t(6,:))); % get active force for muscle 2

F_tPassive1 = F_t1 - F_tActive1;
F_tPassive2 = F_t2 - F_tActive2;

if showFigs
% make an x_t for plotting. just need to add some extra activations
x_tplot = NaN(10,N+1);
x_tplot([1,2,5:7,10],:) = x_t;

plotMuscleAntagonistResults(t_track,x_tplot,u_inp,[[NaN;NaN],[l1_track;l2_track]],...
    [F_tActive1;F_tPassive1;F_tActive2;F_tPassive2],Psim,true);

% check smooth approximation of power
if W.pow > 0
    figure;
    subplot(2,1,1)
    plot(t,[abs(Psim.Vmax*sol.value(Power1));Psim.Vmax*sol.value(absPower1)],'-')
    legend('abs Power','smooth abs Power')
    ylabel('Power (W)')
    title('Antagonist')
    
    subplot(2,1,2)
    plot(t,[abs(Psim.Vmax*sol.value(Power2));Psim.Vmax*sol.value(absPower2)],'-')
    legend('abs Power','smooth abs Power')
    title('Agonist')
    ylabel('Power (W)')
    xlabel('Time (s)')
end
end



disp(['Cost is: ',num2str(J_sol)])

if success && ~isempty(saveName)
    save(saveName,'x_t','u_inp','*_track','Start*','End*','X_guess','U_guess','Psim','D','N','T','*_sol','F_t*','user_input');
end