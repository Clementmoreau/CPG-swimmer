% Motion of a swimmer driven by CPG
% Updated 24 Jan 2024 with omega turn term

addpath(genpath('./helpers'))

%% Set the parameters.

% *** Filament parameters. ***

N = 10; % number of links
gamma = 1/2; % ratio between the hydrodynamic drag coefficients (1/2 for Stokes flow, 1/70 for agar gel)
Sp = 6; % Sperm number (typical range 1-10)
kd = 1; % Bending stiffness (default = 1)
L = ones(1,N); % Lengths of the segments (default = 1)

% *** CPG model parameters ***

% for core model
omega = 2*pi; % T=1 unit for intrinsic oscillation
tau = 10; % Strength of the active torque
coup = 2; % Oscillators coupling strength
coupR = coup; % Differential coupling (not in use currently)

% for variable proprioception
sigma_amp = 0.5; % strength of proprioception
t_switch =50; % how often proprioception sign switches
switch_width = 3; % how fast the switch occurs
sigma = @(t) 2*sigma_amp*sigma_custom(t,t_switch,switch_width) - sigma_amp; % propioception function
sigma_1 = @(t) 15*cos(t/t_switch/sqrt(2))^3;
psi = @(t) pi; % phase in proproceptive term %psi_custom(t,t_switch,switch_width);

% for omega-turns
cutoff = @(x) 2*(sigma_amp - abs(x)); % Cutoff function activating the omega turn when sigma is close to zero
sigma_2 = @(t) 0.2*sin(t/t_switch)^3;
alpha_omega = 0.9*tau/kd;%2*(rand-1/2) * tau/kd; % target angle at each junction, maximum absolute value is tau/kd
K_omega = 5*omega; % how fast the target angle will be reached

% Pack the physical parameters in the params structure.
params = struct();
params.N = N;params.gamma = gamma;params.Sp=Sp;params.kd=kd;params.L=L;
params.omega=omega;params.tau=tau;params.coup=coup;params.coupR=coupR;params.sigma=sigma;params.psi=psi;
params.sigma_amp = sigma_amp;
params.sigma_1 = sigma_1;
params.sigma_2 = sigma_2;
params.cutoff = cutoff; params.K_omega = K_omega; params.alpha_omega = alpha_omega;

% Simulation options
T=650; % final time
tpnum=T*50; % number of time steps
tps=linspace(0,T,tpnum); % time step vector
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Plot the sigma function. 
SIG =[];
for i=1:length(tps)
    SIG(i) = sigma(tps(i));
end
figure(2);clf;plot(tps,SIG);

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

iframe=15;

%-------- draw trajectories
if Itraj==1
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

     

dpi = '-r400';
%saveas(fig5,'fig5','epsc')
%saveas(fig6,'fig6','epsc')
%saveas(fig7,'fig7','epsc') 
%saveas(fig8,'fig8','epsc') 
%saveas(fig9,'fig9','epsc') 

end


% ----------------------------
% Swimming Movie [Imovie=1]
% ----------------------------
if IMovie==1
figure(1);clf

ltot=1;
xmin=min(traj(:,1))-ltot;
xmax=max(traj(:,1))+ltot;
ymin=min(traj(:,2))-ltot;
ymax=max(traj(:,2))+ltot;

ifig=1;
for i=1:3:length(tps)
    plot(traj(1:i,1),traj(1:i,2),'LineWidth',1.3)%trajectory of (x,y)
    hold on
    %title('Trajectory','FontSize',14);
    xlabel('$x$','FontSize',14,'Interpreter','latex');
    ylabel('$y$','FontSize',14,'Interpreter','latex');
    [Xr,Yr,TH]=coordinates_filament(traj(i,:),params);
    %[Xr,Yr]=coordinates_swimmer(traj(i,:));
    plot(Xr,Yr,'k','LineWidth',2)
    axis equal
    axis([xmin xmax ymin ymax])
    set(gca,'FontSize',24)
    set(gca,'TickLabelInterpreter','latex')
    hold off
    drawnow;
    if IMovie==1
    M1(ifig)=getframe(gcf);
%    close
    clf
    ifig=ifig+1;
    end
end

disp('Now.....converting to a MP4 file....')
fps=60;
animtitle=strcat('anim');
vid1aj = VideoWriter(animtitle,'MPEG-4');
vid1aj.FrameRate=fps;
vid1aj.Quality=100;
open(vid1aj)
for frame = 1:size(M1,2)
   currFrame = M1( frame );   
   writeVideo( vid1aj, currFrame );  
end
close( vid1aj );
disp('.....Completed!')
close;
end
