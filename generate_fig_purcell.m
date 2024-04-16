function out=generate_fig_purcell()

clear all
% -- Filament parameters --
% N: number of links
% gamma: ratio between the hydrodynamic drag coefficients 
%        (1/2 for Stokes flow, 1/70 for agar gel)
%
% Sp: Sperm number (typical range 1-10)
% kd: Bending stiffness (default = 1)
% L: Lengths of the segments (default = 1)

N = 3; 
gamma = 1/2; 
Sp = 4; 
kd = 1; 
L = ones(1,N); 

% -- CPG model parameters --
% omega: % (default = 2*pi, T=1 unit)
% tau: Strength of the activity
% coup: Oscillators coupling strength (scaling by N^2?)
% coupR: Differential coupling
% sigma_amp: sensitivity strength

omega = 2*pi; 
tau = 8; 
coup = 1; 
coupR = coup; 
sigma_amp = 4;

% --- Simulation options --
% t_switch: sigma switching time interval
% switch_width: time period for switching (0=step func)
% T: simulation final time
% tpnum: time steps in total

t_switch = 100;
switch_width = 0;
T=6; 
tpnum=1000*T;

sigma = @(t) sigma_amp;
psi = @(t) psi_custom(t,t_switch,switch_width);

Iload=0; %1: load init_guess

% -- pack the parameters in params and opts ---
params = struct();
params.N = N;params.gamma = gamma;params.Sp=Sp;params.kd=kd;params.L=L;
params.omega=omega;params.tau=tau;params.coup=coup;params.coupR=coupR;params.sigma=sigma;params.psi=psi;

opts = struct();
opts.t_switch=t_switch; opts.switch_width=switch_width; opts.T=T; opts.tpnum=tpnum;

% -- inits ----%

Cinit=0.42;
Kinit=1.0;
Pinit=-0.99;
%z0 = [0;0;2*pi*rand;0*ones(N-1,1)]; % initial geometry
%z0 = [0;0;0;0*pi*ones(N-1,1)]; % initial geometry
z0 = [0;0;0;Cinit*cos((2*pi*Kinit/(N-1))*(1:N-1)')]; % initial geometry
p0 = Pinit*(2*pi/(N-1))*(1:N-1)'; % initial CPG state

init = [z0;p0];

if Iload==1
    load('init_guess.mat');
    Ng=(length(init)-3)/2+1;
    Pertb=0.0*rand(1,N-1); % pertubation in phase variables
    ThgS=interp1(linspace(0,1,Ng),[0;cumsum(init(4:Ng+2))],linspace(0,1,N),'linear');
    Thg=[diff(ThgS)];
    Phg=interp1(linspace(0,1,Ng-2),init(Ng+4:2*Ng+1),linspace(0,1,N-2),'linear');
    pp0=[Thg';Phg'];
    init=[0;0;0;Thg'+Pertb';0;Phg'];
end

inits=init;

% -- generate figures ----%

% fig_sims_swimL
[tps, traj]=generate_fig_main(params, inits, opts);
generate_fig_purcell_draw(tps, traj, params);

end
