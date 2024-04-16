function [tps, traj]=generate_fig_main(params, inits, opts)


% Motion of a swimmer driven by CPG

addpath(genpath('./helpers'))


% Unpack the physical parameters in the params structure.
N= params.N;
gamma=params.gamma;
Sp=params.Sp;
kd=params.kd;
L=params.L;
omega=params.omega;
tau=params.tau;
coup=params.coup;
coupR=params.coupR;
sigma=params.sigma;
psi=params.psi;

% Unpack the options from the opts structure.

t_switch=opts.t_switch;
switch_width=opts.switch_width;
T=opts.T;
tpnum=opts.tpnum;

% load initial data.
init=inits;

%% Set the parameters.

%sigma_amp=sigma;


% Simulation options
tps=linspace(0,T,tpnum); % time step vector
option = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Input Optiions.
Isym=1; % 0: asym, 1:sym for the proprioceptice control


%Initial condition. (x,y,theta,alpha_1...alpha_N-1)


%% Solve the equation.

if Isym==1
    dZ=@(t,z) cpg_sym(t,z,params); % symmetrized proprioceptive control
else
    dZ=@(t,z) cpg(t,z,params);   % asymmetric forward proprioceptive control
end

tic
[tps,traj]=ode15s(dZ,tps,init,option);
toc

%{

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

%-------- draw trajectories
if Itraj==1
fig7=figure(7);clf;
set(gcf, 'Position',  [1, 1, 3*figsize, 3*figsize])
for i=1:4:length(tps)%i=length(tps)-1000:10:length(tps)
    plot(Xc(1:i),Yc(1:i),'k','LineWidth',2)
    hold on
    plot(X1(1:i),Y1(1:i),'k')
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'b','LineWidth',6)
    hold off 
    box on
    axis equal
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex') 
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
    box on
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
    box on
    set(gca,'FontSize',20)
    hold off 


    phi1=traj(:,N+3);
    phi2=traj(:,N+4);
    the1=traj(:,4);
    the2=traj(:,5);    
fig10=figure(10);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    hold on
    plot(phi1/pi,phi2/pi,'LineWidth',2)  
    hold off
    xlabel('$\phi_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\phi_2/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on  
    box on
    set(gca,'FontSize',20)
    hold off 

fig11=figure(11);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    hold on
    scatter3(the1/pi,the2/pi, mod(phi2/pi,2),20,tps,'filled')  
    hold off
    xlabel('$\alpha_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\alpha_2/\pi$','FontSize',14,'Interpreter','latex');
    zlabel('$\phi_2/\pi$','FontSize',14,'Interpreter','latex');    
    set(gca,'TickLabelInterpreter','latex')    
    axis equal
    cb=colorbar;
    colormap jet;
    set(cb, "TickLabelInterpreter",'latex');
    grid on    
    box on
    set(gca,'FontSize',20)
    hold off  


%fig6=figure(6);clf;
%set(gcf, 'Position',  [figsize, 640, figsize*2, figsize])
    Larc=cumsum(L)/sum(L);
    Tarc=tps;
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc(1:end-1));
    Adat=traj(:,4:N+2);
    Phdat=mod(traj(:,N+3:end),2*pi);

    Larc1=linspace(0,1,N);
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc1);
    Adat=cat(2,traj(:,4:N+2),traj(:,N+2));
    Phdat=cat(2,mod(traj(:,N+3:end),2*pi),mod(traj(:,end),2*pi));

fig5=figure(5);clf;
set(gcf, 'Position',  [figsize, 640, figsize*2, figsize*0.5])    
    sA=surf(Tmesh,Lmesh,Adat');
    sA.EdgeColor = 'none';
    %image(PhData','CDataMapping','scaled');
    cbA=colorbar;
    colormap jet
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    ylim([0,1])
    title('curvature','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on    
    box on
    set(gca,'FontSize',20)

fig6=figure(6);clf;
set(gcf, 'Position',  [figsize, 300, figsize*2, figsize*0.5])    
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
    box on
    set(gca,'FontSize',20)

     

dpi = '-r400';
saveas(fig5,'fig5','epsc')
saveas(fig6,'fig6','epsc')
saveas(fig7,'fig7','epsc') 
saveas(fig8,'fig8','epsc') 
saveas(fig9,'fig9','epsc') 


saveas(fig10,'fig10','epsc') 
saveas(fig11,'fig11','epsc') 

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

%}

end
