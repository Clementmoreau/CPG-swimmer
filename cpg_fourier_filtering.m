% Motion of a swimmer driven by CPG

addpath(genpath('./helpers'))

%% Set the parameters.

% Filament parameters.
N = 31; % number of links
gamma = 1/70; % ratio between the hydrodynamic drag coefficients (1/2 for Stokes flow, 1/70 for agar gel)
Sp = 1.5; % Sperm number (typical range 1-10)
kd = 1; % Bending stiffness (default = 1)
L = ones(1,N); % Lengths of the segments (default = 1)

% CPG model parameters.
omega = 2*pi; % T=1 unit
tau = 10; % Strength of the activity
coup = 0.2; % Oscillators coupling strength
coupR = coup; % Differential coupling
sigma_amp = -1;
sigma = @(t) sigma_amp; %* 2*(mod(floor(t/8),2)-3/4); % Proprioception strength
%psi = @(t) pi*(mod(floor(t/10),2));

t_switch = 8;
switch_width = 0.5;
psi = @(t) psi_custom(t,t_switch,switch_width);
% psi = @(t) pi*t/10;
%psi = @(t) 0;

% sig_test = sigma*2*(mod(floor(t/5),2)-1/2);

% Pack the physical parameters in the params structure.
params = struct();
params.N = N;params.gamma = gamma;params.Sp=Sp;params.kd=kd;params.L=L;
params.omega=omega;params.tau=tau;params.coup=coup;params.coupR=coupR;params.sigma=sigma;params.psi=psi;

% Simulation options
T=6; % final time
tpnum=T*250; % number of time steps
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

% Using the symmetric formulation.
dZ=@(t,z) cpg_sym(t,z,params);
tic
[tps,traj]=ode15s(dZ,tps,init,opts);
toc


%% visualisation of the effect of fourier filtering

% Maximum number of modes we will keep. 
    nb_modes = min(ceil((N-1)/2),10);

    % prepare the figures.
    figure(5);clf;tiledlayout(4,1)
    set(gcf,'Position',[1 1 600 1500]);
    figure(6);clf;tiledlayout(4,1)
    set(gcf,'Position',[600 1 600 1500]);

% Loop over number of modes we keep, from 1 to the maximal number.
   for i_filt = 12:14

    % careful of taking the transpose of traj because fft is performed on the
    % colums by default.

    % Fourier transform of the shape.
    shape_filtered = filter_high_mode(traj(:,4:N+2)',i_filt);

    % Should it be taken over phi with mod 2*pi? probably not but...
    % phi_filtered = filter_high_mode(mod(traj(:,N+3:end)',2*pi),i_filt);

    % Fourier transform of the phase.
    phi_filtered = filter_high_mode(mod(traj(:,N+3:end)',2*pi),i_filt);

    Larc=cumsum(L)/sum(L);
    Tarc=tps;
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc(1:end-1));

    % Plot the whole filtered shape evolution.
    figure(5);nexttile(i_filt-11);
    sA=surf(Tmesh,Lmesh,shape_filtered);
    sA.EdgeColor = 'none';
    cbA=colorbar;
    colormap jet
    %xlabel('time','FontSize',14,'Interpreter','latex');
    %ylabel('length','FontSize',14,'Interpreter','latex');
    title(['curvature, ',num2str(i_filt),' modes kept'],'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on    

    % Plot the whole filtered shape evolution.
    figure(6);nexttile(i_filt-11);  
    sPh=surf(Tmesh,Lmesh,mod(phi_filtered,2*pi));
    sPh.EdgeColor = 'none';
    cbPh=colorbar;
    colormap hsv
    %xlabel('time','FontSize',14,'Interpreter','latex');
    %ylabel('length','FontSize',14,'Interpreter','latex');
    title(['phase, ',num2str(i_filt),' modes'],'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    set(cbPh,'TickLabelInterpreter','latex'); 
    set(cbPh,'Limits',[0 2*pi]); 
    set(cbPh,'Ticks',[0 pi 2*pi]);
    set(cbPh,'TickLabels',{'0', '$\pi$', '$2\pi$'});       
    axis tight
    view(0,90)
    grid on    

    drawnow
   end

   % plot the original one for comparison
   figure(5);nexttile;
   sA=surf(Tmesh,Lmesh,traj(:,4:N+2)');
    sA.EdgeColor = 'none';
    cbA=colorbar;
    colormap jet
    %xlabel('time','FontSize',14,'Interpreter','latex');
    %ylabel('length','FontSize',14,'Interpreter','latex');
    title(['curvature without filtering'],'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on 

    figure(6);nexttile;
    sA=surf(Tmesh,Lmesh,mod(traj(:,N+3:end)',2*pi));
    sA.EdgeColor = 'none';
    cbA=colorbar;
    colormap hsv
    %xlabel('time','FontSize',14,'Interpreter','latex');
    %ylabel('length','FontSize',14,'Interpreter','latex');
    title(['phase without filtering'],'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on 

%% usual figures

% usual calculations of swimmer coordinates
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

%-------- draw trajectories
fig7=figure(7);clf;
set(gcf, 'Position',  [1, 1, 3*figsize, 3*figsize])
for i=1:iframe:length(tps)%i=length(tps)-1000:10:length(tps)
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

% draw snapshots of swimmer
fig8=figure(8);clf;
set(gcf, 'Position',  [figsize, 640, figsize, figsize])
    hold on
for i=1:100:length(tps)%floor(0.1*length(tps)):17:length(tps)
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'LineWidth',3)    
end
    xlabel('$x$','FontSize',14,'Interpreter','latex');
    ylabel('$y$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex') 
    axis equal
    grid on    
    set(gca,'FontSize',20)
    hold off 

% time series of each oscillator
fig9=figure(9);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    phi1=traj(:,N+3);
    phi2=traj(:,N+4);
    phiD=mod(diff(traj(:,N+3:end),1,2)+pi,2*pi)-pi;
    phiav=mean(phiD,2);
    plot(tps,phiD(:,1)/pi,tps, phiav/pi,'LineWidth',3)  
    hold on
    for i=1:N-2
        plot(tps,phiD(:,i)/pi,'LineWidth',1) 
    end
    hold off
    xlabel('$t$','FontSize',14,'Interpreter','latex');
    ylabel('$\phi_i/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on    
    set(gca,'FontSize',20)
    hold off 




