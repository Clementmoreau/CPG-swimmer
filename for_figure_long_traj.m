% Motion of a swimmer driven by CPG
% Updated 24 Jan 2024 with omega turn term

addpath(genpath('./helpers'))

%% Set the parameters.

% *** Filament parameters. ***

N = 10; % number of links
gamma = 1/70; % ratio between the hydrodynamic drag coefficients (1/2 for Stokes flow, 1/70 for agar gel)
Sp = 4; % Sperm number (typical range 1-10)
kd = 1; % Bending stiffness (default = 1)
L = ones(1,N); % Lengths of the segments (default = 1)

% *** CPG model parameters ***

% for core model
omega = 2*pi; % T=1 unit for intrinsic oscillation
tau = 8; % Strength of the active torque
coup = 0.1; % Oscillators coupling strength
coupR = coup; % Differential coupling (not in use currently)

% for variable proprioception
t_switch = 50; % how often proprioception sign switches
sigma_1 = @(t) 12*(cos(2*pi*t/t_switch)^1);

% for omega-turns
sigma_2 = @(t) 3*sin(2*pi*t/t_switch/sqrt(2))^25;
alpha_omega = 0.99*tau/kd;%2*(rand-1/2) * tau/kd; % target angle at each junction, maximum absolute value is tau/kd
K_omega = 2*omega; % how fast the target angle will be reached

% Pack the physical parameters in the params structure.
params = struct();
params.N = N;params.gamma = gamma;params.Sp=Sp;params.kd=kd;params.L=L;
params.omega=omega;params.tau=tau;params.coup=coup;params.coupR=coupR;
params.sigma_1 = sigma_1;
params.sigma_2 = sigma_2;
params.K_omega = K_omega; params.alpha_omega = alpha_omega;

% Simulation options
T=400; % final time
tpnum=T*50; % number of time steps
tps=linspace(0,T,tpnum); % time step vector
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Plot the sigma_1-sigma_2 function. 
sig1 =zeros(1,length(tps));
sig2 =zeros(1,length(tps));
for i=1:length(tps)
    sig1(i) = sigma_1(tps(i));
    sig2(i) = sigma_2(tps(i));
end
figure(2);clf;plot(sig1,sig2);
%axis equal

% Output options.
Itraj=1; % 1: to visualise trajecotry
IMovie=0; % 1: to generate a movie

%Initial condition. (x,y,theta,alpha_1...alpha_N-1)

% z0 = [0;0;2*pi*rand;0*ones(N-1,1)]; % initial geometry
z0 = [0;0;0;0*ones(N-1,1)]; 

p0 = zeros(N-1,1); % initial CPG state

init = [z0;p0];

%% Solve the equation.

dZ=@(t,z) cpg_sym(t,z,params);
tic
[tps,traj]=ode15s(dZ,tps,init,opts);
toc

%% post-processing

Xc=[];Yc=[];
X1=[];Y1=[];
for i = 1:length(tps)
     [X,Y,TH]=coordinates_filament(traj(i,:),params);
     Xc=[Xc sum(X)/(N+1)];Yc=[Yc sum(Y)/(N+1)];
     X1=[X1 X(1)];
     Y1=[Y1 Y(1)];
end

figsize = 400;

iframe=15;

%-------- draw trajectories
fig7=figure(7);clf;
set(gcf, 'Position',  [1648,-322,946,894])
    hold on
    plot(X1(1:i),Y1(1:i),'k')    
    axis equal
    set(gca,'FontSize',20)
    box on

snap = 350;

set(gcf, 'Position',  [1648,-322,434,434])
    hold on

    col = colormap(othercolor('Paired6',length(tps)));

for i=1:snap:length(tps)%floor(0.1*length(tps)):17:length(tps)
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'LineWidth',3,'Color',col(i,:))    
end
    xlabel('$x$','FontSize',14,'Interpreter','latex');
    ylabel('$y$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex') 
    axis equal
    grid on    
    set(gca,'FontSize',20)
    hold off 

    c = colorbar;
    c.Ticks=[];

    %exportgraphics(gcf,'fig_turn_long.pdf','ContentType','vector')

