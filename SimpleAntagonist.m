function [forward_dt,F_fun,F_fun_mscl] =  SimpleAntagonist(Psim)

import casadi.*

% state, inputs, parameter vectors
l1 = MX.sym('l1',2); a1 = MX.sym('a1',1); a1_12 = MX.sym('a1_12',2);
u1 = MX.sym('u1',1);

l2 = MX.sym('l2',2); a2 = MX.sym('a2',1); a2_12 = MX.sym('a2_12',2);
u2 = MX.sym('u2',1);

% extract Parameters

% simulation (physical) parameters
M = Psim.M;
Fmax = Psim.Fmax;
Vmax = Psim.Vmax;
lO = Psim.lO;

% simulation time parameters
dt = Psim.dt; % time increment

% Muscle model parameters
b = Psim.b;
c = Psim.c;
d = Psim.d;
alp = Psim.alp;

% Specify hill-type model based on Murtola 2023

if mod(b(1),2) == 0
    fl = Otten_fl_even(l1(1)/lO,b);
else
    fl = Otten_fl_smooth(l1(1)/lO,b,200);
end
fv = Otten_fv_smooth(-l1(2)/lO/Vmax,d,200);
fp = Murtola_fp_smooth(l1(1)/lO,c,200);

F1mscl = Fmax*(a1(1)*(fl.*fv));
F_fun_mscl = Function('F',{l1;a1},{F1mscl});

F1 = F1mscl+Fmax*fp;
F_fun = Function('F',{l1;a1},{F1});

F2 = F_fun(l2,a2);


% set dynamics based on activation type

if ~isfield(Psim,'actOrder')
    Psim.actOrder = 3;
end
switch Psim.actOrder
    case 1
        s = Psim.s; % smoothing parameter for discontinuity in o1 dynamics

        % RHS
        x_dot1 = [l1(2); (-F1+F2)/M; o1ActDynSmooth(a1,u1,alp,s)];
        x_dot2 = [l2(2); (F1-F2)/M; o1ActDynSmooth(a2,u2,alp,s)];
        x_dot = [x_dot1;x_dot2];
        
        %simulation
        ode = struct('x',[l1;a1;l2;a2],'p',[u1;u2],'ode',x_dot);
        opts = struct('tf',dt,'number_of_finite_elements',1);
        intg = integrator('intg','rk', ode, opts);
        
        % create function to run a short simulation and output the result
        X0 = MX.sym('X0',6);
        P = MX.sym('P',2);
        res = intg('x0',X0,'p',P); % symbolic evaluation of integral        
    case 3
        bet = Psim.bet;
        
        % RHS
        x_dot1 = [l1(2); (-F1+F2)/M; o3ActDyn([a1_12;a1],u1,alp,bet)];
        x_dot2 = [l2(2); (F1-F2)/M; o3ActDyn([a2_12;a2],u2,alp,bet)];
        x_dot = [x_dot1;x_dot2];
        % Muscles pull against eachother and connecting mass
        % O3 dynamics according to Murtola and Richards 2023
        
        %simulation
        ode = struct('x',[l1;a1_12;a1;l2;a2_12;a2],'p',[u1;u2],'ode',x_dot);
        opts = struct('tf',dt,'number_of_finite_elements',1);
        intg = integrator('intg','rk', ode, opts);
        
        % create function to run a short simulation and output the result
        X0 = MX.sym('X0',10);
        P = MX.sym('P',2);
        res = intg('x0',X0,'p',P); % symbolic evaluation of integral
    otherwise
        error('Unidentified order of activation dynamics')
end






forward_dt = Function('FD_dt',{X0;P},{res.xf});
