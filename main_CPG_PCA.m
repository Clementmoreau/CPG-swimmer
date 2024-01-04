% Motion of a swimmer driven by CPG
% with PCA analysis

addpath(genpath('./helpers'))

%% Set the parameters.

% Filament parameters.
N = 10; % number of links
gamma = 1/70; % ratio between the hydrodynamic drag coefficients (1/2 for Stokes flow, 1/70 for agar gel)
Sp = 4; % Sperm number (typical range 1-10)
kd = 1; % Bending stiffness (default = 1)
L = ones(1,N); % Lengths of the segments (default = 1)

% CPG model parameters.
omega = 2*pi; % T=1 unit
tau = 8;%*2^(1/4); % Strength of the activity
coup = 1; % Oscillators coupling strength
coupR = coup; % Differential coupling
sigma_amp = 1;
sigma = @(t) sigma_amp; %* 2*(mod(floor(t/8),2)-3/4); % Proprioception strength
%psi = @(t) pi*(mod(floor(t/10),2));

t_switch = 20;
switch_width = 2;
psi = @(t) psi_custom(t,t_switch,switch_width);
% psi = @(t) pi*t/10;
%psi = @(t) 0;

% sig_test = sigma*2*(mod(floor(t/5),2)-1/2);

% Pack the physical parameters in the params structure.
params = struct();
params.N = N;params.gamma = gamma;params.Sp=Sp;params.kd=kd;params.L=L;
params.omega=omega;params.tau=tau;params.coup=coup;params.coupR=coupR;params.sigma=sigma;params.psi=psi;

% Simulation options
T=400; % final time
tpnum=T*100; % number of time steps
tps=linspace(0,T,tpnum); % time step vector
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

PSI =[];
for i=1:length(tps)
    PSI(i) = psi(tps(i));
end
figure(2);clf;plot(tps,PSI);

% Output options.
Itraj=1; % 1: to visualise trajecotry
IMovie=0; % 1: to generate a movie

%Initial condition. (x,y,theta,alpha_1...alpha_N-1)

% z0 = [0;0;2*pi*rand;0*ones(N-1,1)]; % initial geometry
z0 = [0;0;0;0*ones(N-1,1)]; 

p0 = (1.5)*(2*pi/(N-1))*(1:N-1)'; % initial CPG state

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

iframe=20;

%% -------- draw trajectories

fig7=figure(7);clf;
set(gcf, 'Position',  [1, 1, 3*figsize, 3*figsize])
for i=length(tps):length(tps)%i=length(tps)-1000:10:length(tps)
    plot(Xc(1:i),Yc(1:i),'k','LineWidth',2)
    hold on
    plot(X1(1:i),Y1(1:i),'k')
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'b','LineWidth',6)
    hold off     
    axis equal
    set(gca,'FontSize',20)
    title(['$T$ =',num2str(tps(i))],'Interpreter','latex')
    drawnow
end

%fig6=figure(6);clf;
%set(gcf, 'Position',  [figsize, 640, figsize*2, figsize])
    Larc=cumsum(L)/sum(L);
    Tarc=tps;
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc(1:end-1));
    
    Phdat=mod(traj(:,N+3:end),2*pi);

fig5=figure(5);clf;
set(gcf, 'Position',  [figsize, 400, figsize*2, figsize*0.5])    
    sA=surf(Tmesh,Lmesh,traj(:,4:N+2)');
    sA.EdgeColor = 'none';
    %image(PhData','CDataMapping','scaled');
    cbA=colorbar;
    colormap jet
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    title('curvature','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on    
    set(gca,'FontSize',20)

    Phdat=mod(traj(:,N+3:end),2*pi);
fig6=figure(6);clf;
set(gcf, 'Position',  [figsize, 100, figsize*2, figsize*0.5])    
    sPh=surf(Tmesh,Lmesh,Phdat');
    sPh.EdgeColor = 'none';
    cbPh=colorbar;
    colormap hsv
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    title('phase','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')    
    set(cbPh,'TickLabelInterpreter','latex'); 
    set(cbPh,'Limits',[0 2*pi]); 
    set(cbPh,'Ticks',[0 pi 2*pi]);
    set(cbPh,'TickLabels',{'0', '$\pi$', '$2\pi$'});       
    axis tight
    view(0,90)
    grid on    
    set(gca,'FontSize',20)

shapedata = traj(:,4:N+2);
[shapedataPCA,score,latent] = pca(shapedata);
figure(1);clf;
tiledlayout(9,1);
for i = 1:9
    nexttile(i);
    plot(tps,score(:,i));
end

figure(2);clf
plot3(score(:,1),score(:,2),score(:,3))

figure(3);clf
plot(tps,sum(score(:,5:end).^2,2))


% smoothen with averaging



