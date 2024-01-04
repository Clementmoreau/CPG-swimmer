% Motion of a swimmer with 'locking' (anisotropic stiffness).
% Figure for small amplitude regime.

close all;
addpath(genpath('./helpers'))

% Number of links.
N = 20; 

% Physical parameters.

Sp = 10; % Sperm number
gamma = 0.5; % Ratio between RFT drag coefficients
L = ones(1,N); % Lengths of the links (default=1)
% Angular actuation at one end.
amp = 1e-3;
% Torque actuation at one end.
% ampT = 15;
% Actuation frequency.
omega = 1;
% Anisotropic stiffness function.
a1 = 1e-7;
a2 = 0;
x0 = 0;
stiffness = @(x) 1 + a2*(1+tanh(-1*(x-x0)/a1))/2 ; % tanh-based
% stiffness = @(x) 1 + ((x-x0)>0).*((x-x0)/a1).^p; % polynomial-based
% stiffness = @(x) 1; % constant

% Integration time and ODE precision.
nb_per = 10;
T=nb_per*2*pi;
tpnum = nb_per*200;
tps=linspace(0,T,tpnum); %time step
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

%initial condition (x,y,theta,alpha_1...alpha_N-1)
z0=[0;0;0;0*ones(N-1,1)];
% z0=[0;0;0;(rand(N-1,1)-1/2)];

% Solve the equation.
% dZ=@(t,z) elasticity_forced_torque(t,z,N,L,gamma,Sp,ampT,omega,stiffness);
dZ=@(t,z) elasticity(t,z,N,L,gamma,Sp,amp,omega,stiffness);
tic
[tps,traj]=ode15s(dZ,tps,z0,opts);
toc

%% Post-processing ang plotting

% Full coordinates at each time step.
Xc=zeros(1,length(tps));Yc=zeros(1,length(tps));
XX=zeros(N+1,length(tps));YY=zeros(N+1,length(tps));
for j = 1:length(tps)
    [X,Y,TH]=coordinates_filament(traj(j,:),N,L);
    Xc(j) = sum(X)/(N+1);
    Yc(j) = sum(Y)/(N+1);
    XX(:,j) = X;
    YY(:,j) = Y;
end

% Prepare the figure.
fig = figure(1);clf;
set(gcf, 'Position',  [1, 100, 1200, 1100])
set(gcf,'color','w');

tl = tiledlayout(3,1);

nexttile(1);
hold on
col = colormap(hsv(200));
for i=1:4:200
    plot(XX(:,end-200+i),YY(:,end-200+i),'LineWidth',2,'Color',col(i,:))
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
%axis equal
box on
grid on
set(gca,'FontSize',20)
 title(['Beating pattern'], 'Interpreter','latex')
