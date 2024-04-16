function out=generate_fig_purcellSPO(tps, traj, params)


% Motion of a swimmer driven by CPG

addpath(genpath('./helpers'))

% Output options.
Itraj=1; % 1: to visualise trajecotry
IMovie=0; % 1: to generate a movie

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
    axis equal
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex') 
    title(['$T$ =',num2str(tps(i))],'Interpreter','latex')
    drawnow
end


fig8=figure(8);clf;
set(gcf, 'Position',  [figsize, 640, figsize, figsize])
    hold on
for i=1:10:length(tps)%floor(0.1*length(tps)):17:length(tps)
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
    ylabel('$\delta\phi_i/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on    
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
    set(gca,'FontSize',20)
    hold off 

fig11=figure(11);clf;
set(gcf, 'Position',  [2*figsize, 120, figsize, figsize])
    hold on
    plot(the1/pi,the2/pi,'LineWidth',2)  
    hold off
    xlabel('$\alpha_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\alpha_2/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on    
    set(gca,'FontSize',20)
    hold off

fig12=figure(12);clf;
set(gcf, 'Position',  [2*figsize, 320, figsize, figsize])
    hold on
    plot(tps,the1/pi,tps,the2/pi, tps, mod(phi2/pi,2),'LineWidth',2)  
    hold off
    xlabel('$t$','FontSize',14,'Interpreter','latex');
    ylabel('$\alpha_1/\pi,\alpha_2/\pi,\phi_2/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    legend('$\alpha_1/\pi$','$\alpha_2/\pi$','$\phi_2/\pi$','Interpreter','latex')
    axis tight
    grid on    
    set(gca,'FontSize',20)
    hold off    

%fig6=figure(6);clf;
%set(gcf, 'Position',  [figsize, 640, figsize*2, figsize])
    Larc=cumsum(L)/sum(L);
    Tarc=tps/Tppe; % normalized by time period of the orbit
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc(1:end-1));
    Adat=traj(:,4:N+2);
    Phdat=mod(traj(:,N+3:end),2*pi);

fig5=figure(5);clf;
set(gcf, 'Position',  [figsize, 640, figsize, figsize])    
    sA=surf(Tmesh,Lmesh,Adat');
    sA.EdgeColor = 'none';
    %image(PhData','CDataMapping','scaled');
    cbA=colorbar;
    colormap jet
    xlabel('$t/T$','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    title('curvature','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on    
    set(gca,'FontSize',20)

fig6=figure(6);clf;
set(gcf, 'Position',  [figsize, 300, figsize, figsize])    
    sPh=surf(Tmesh,Lmesh,Phdat');
    sPh.EdgeColor = 'none';
    cbPh=colorbar;
    colormap hsv
    xlabel('$t/T$','FontSize',14,'Interpreter','latex');
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
saveas(fig4,'fig_purcellSP0','epsc')
saveas(fig5,'fig_purcellSP1','epsc')
saveas(fig6,'fig_purcellSP2','epsc')
saveas(fig7,'fig_purcellSP3','epsc') 
saveas(fig8,'fig_purcellSP4','epsc') 
saveas(fig9,'fig_purcellSP5','epsc') 

saveas(fig10,'fig_purcellSP6','epsc') 
saveas(fig11,'fig_purcellSP7','epsc') 
saveas(fig12,'fig_purcellSP8','epsc') 

end



end
